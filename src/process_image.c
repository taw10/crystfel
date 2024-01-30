/*
 * process_image.c
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2022 Thomas White <taw@physics.org>
 *   2014-2017 Valerio Mariani <valerio.mariani@desy.de>
 *   2017      Stijn de Graaf
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <unistd.h>
#include <sys/stat.h>

#include <utils.h>
#include <index.h>
#include <peaks.h>
#include <filters.h>
#include <thread-pool.h>
#include <geometry.h>
#include <stream.h>
#include <reflist-utils.h>
#include <integration.h>
#include <detgeom.h>
#include <image-msgpack.h>
#include <profile.h>

#include "process_image.h"
#include "predict-refine.h"
#include "im-sandbox.h"
#include "im-zmq.h"
#include "im-asapo.h"
#include "peaks.h"
#include "peakfinder8.h"

static float **backup_image_data(float **dp, struct detgeom *det)
{
	float **bu;
	int i;

	bu = malloc(det->n_panels * sizeof(float *));
	if ( bu == NULL ) return NULL;

	for ( i=0; i<det->n_panels; i++ ) {

		size_t data_size;

		data_size = det->panels[i].w * det->panels[i].h * sizeof(float);
		bu[i] = malloc(data_size);
		if ( bu[i] == NULL ) {
			free(bu);
			ERROR("Failed to allocate pre-filter backup.\n");
			return NULL;
		}

		memcpy(bu[i], dp[i], data_size);

	}

	return bu;
}


static void restore_image_data(float **dp, struct detgeom *det, float **bu)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {
		size_t data_size;
		data_size = det->panels[i].w * det->panels[i].h * sizeof(float);
		memcpy(dp[i], bu[i], data_size);
		free(bu[i]);
	}
	free(bu);
}


static struct image *file_wait_open_read(const char *filename,
                                         const char *event,
                                         DataTemplate *dtempl,
                                         struct sb_shm *sb_shared,
                                         char *last_task,
                                         signed int wait_for_file,
                                         int cookie,
                                         int no_image_data,
                                         int no_mask_data,
                                         ImageDataArrays *ida)
{
	signed int file_wait_time = wait_for_file;
	int wait_message_done = 0;
	int read_retry_done = 0;
	int r;
	struct image *image;

	set_last_task(last_task, "wait for file");

	do {

		struct stat statbuf;

		sb_shared->pings[cookie]++;
		r = stat(filename, &statbuf);
		if ( r ) {

			if ( (wait_for_file != 0) && (file_wait_time != 0) ) {

				if ( !wait_message_done ) {
					STATUS("Waiting for '%s'\n", filename);
					wait_message_done = 1;
				}

				sleep(1);
				if ( wait_for_file != -1 ) {
					file_wait_time--;
				}
				continue;

			}

			ERROR("File not found: %s (process_image)\n", filename);
			return NULL;
		}

	} while ( r );

	do {

		set_last_task(last_task, "read file");
		sb_shared->pings[cookie]++;

		profile_start("image-read");
		image = image_read(dtempl, filename, event,
		                   no_image_data, no_mask_data, ida);
		profile_end("image-read");
		if ( image == NULL ) {
			if ( wait_for_file && !read_retry_done ) {
				read_retry_done = 1;
				STATUS("File '%s' exists but could not be read."
				       "  Trying again after 10 seconds.\n",
				       filename);
				sleep(10);
				continue;
			}
			ERROR("Couldn't read image: %s\n", filename);
			return NULL;
		}

	} while ( image == NULL );

	return image;
}


void process_image(const struct index_args *iargs, struct pattern_args *pargs,
                   Stream *st, int cookie, const char *tmpdir,
                   int serial, struct sb_shm *sb_shared,
                   char *last_task, struct im_asapo *asapostuff,
                   Mille *mille, ImageDataArrays *ida)
{
	struct image *image;
	int i;
	int r;
	int ret;
	char *rn;
	float **prefilter;
	int any_crystals;

	if ( pargs->zmq_data != NULL ) {

		set_last_task(last_task, "unpacking ZMQ data");
		profile_start("read-zmq-data");
		image = image_read_data_block(iargs->dtempl,
		                              pargs->zmq_data,
		                              pargs->zmq_data_size,
		                              NULL,
		                              iargs->data_format,
		                              serial,
		                              iargs->no_image_data,
		                              iargs->no_mask_data,
		                              ida);
		profile_end("read-zmq-data");
		if ( image == NULL ) return;

		/* image_read_data_block() will leave the filename/event as
		 * NULL, because there's no file (duh).  Fill them in now with
		 * the values passed down to us. For ZMQ, these values are just
		 * placeholders. */
		image->filename = strdup(pargs->filename);
		image->ev = strdup(pargs->event);

	} else if ( pargs->asapo_data != NULL ) {

		set_last_task(last_task, "unpacking ASAP::O data");
		profile_start("read-asapo-data");
		image = image_read_data_block(iargs->dtempl,
		                              pargs->asapo_data,
		                              pargs->asapo_data_size,
		                              pargs->asapo_meta,
		                              iargs->data_format,
		                              serial,
		                              iargs->no_image_data,
		                              iargs->no_mask_data,
		                              ida);
		profile_end("read-asapo-data");
		if ( image == NULL ) return;

		/* image_read_data_block() will leave the filename/event as
		 * NULL, because there's no file (duh).  Fill them in now with
		 * the values passed down to us from ASAP::O. */
		image->filename = strdup(pargs->filename);
		image->ev = strdup(pargs->event);

	} else {
		profile_start("file-wait-open-read");
		image = file_wait_open_read(pargs->filename, pargs->event,
		                            iargs->dtempl,
		                            sb_shared, last_task,
		                            iargs->wait_for_file,
		                            cookie,
		                            iargs->no_image_data,
		                            iargs->no_mask_data,
		                            ida);
		profile_end("file-wait-open-read");
		if ( image == NULL ) {
			if ( iargs->wait_for_file != 0 ) {
				pthread_mutex_lock(&sb_shared->totals_lock);
				sb_shared->should_shutdown = 1;
				pthread_mutex_unlock(&sb_shared->totals_lock);
			}
			return;
		}
	}

	image->serial = serial;

	/* Take snapshot of image before applying horrible noise filters */
	set_last_task(last_task, "image filter");
	profile_start("image-filter");
	sb_shared->pings[cookie]++;

	if ( (iargs->peak_search.median_filter > 0) || iargs->peak_search.noisefilter ) {
		profile_start("data-backup");
		prefilter = backup_image_data(image->dp, image->detgeom);
		profile_end("data-backup");
	} else {
		prefilter = NULL;
	}

	if ( iargs->peak_search.median_filter > 0 ) {
		profile_start("median-filter");
		filter_median(image, iargs->peak_search.median_filter);
		profile_end("median-filter");
	}

	if ( iargs->peak_search.noisefilter ) {
		profile_start("median-filter");
		filter_noise(image);
		profile_end("noise-filter");
	}
	profile_end("image-filter");

	set_last_task(last_task, "resolution range");
	sb_shared->pings[cookie]++;
	mark_resolution_range_as_bad(image, iargs->highres, +INFINITY);

	sb_shared->pings[cookie]++;
	profile_start("peak-search");
	switch ( iargs->peak_search.method ) {

		ImageFeatureList *peaks;

		case PEAK_HDF5:
		case PEAK_CXI:
		set_last_task(last_task, "peaksearch:hdf5orcxi");
		peaks = image_read_peaks(iargs->dtempl,
		                         pargs->filename,
		                         pargs->event,
		                         iargs->peak_search.half_pixel_shift);
		if ( iargs->peak_search.revalidate ) {
			ImageFeatureList *npeaks = validate_peaks(image, peaks,
			                                          iargs->peak_search.min_snr,
			                                          iargs->peak_search.pk_inn,
			                                          iargs->peak_search.pk_mid,
			                                          iargs->peak_search.pk_out,
			                                          iargs->peak_search.use_saturated,
			                                          iargs->peak_search.check_hdf5_snr);
			image_feature_list_free(peaks);
			image->features = npeaks;
		}
		break;

		case PEAK_ZAEF:
		set_last_task(last_task, "peaksearch:zaef");
		image->features = search_peaks(image, iargs->peak_search.threshold,
		                               iargs->peak_search.min_sq_gradient,
		                               iargs->peak_search.min_snr,
		                               iargs->peak_search.pk_inn,
		                               iargs->peak_search.pk_mid,
		                               iargs->peak_search.pk_out,
		                               iargs->peak_search.use_saturated);
		break;

		case PEAK_PEAKFINDER8:
		set_last_task(last_task, "peaksearch:pf8");
		image->features = peakfinder8(image, 2048,
		                              iargs->peak_search.threshold,
		                              iargs->peak_search.min_snr,
		                              iargs->peak_search.min_pix_count,
		                              iargs->peak_search.max_pix_count,
		                              iargs->peak_search.local_bg_radius,
		                              iargs->peak_search.min_res,
		                              iargs->peak_search.max_res,
		                              iargs->peak_search.use_saturated,
		                              iargs->peak_search.peakfinder8_fast,
		                              iargs->pf_private);
		break;

		case PEAK_PEAKFINDER9:
		set_last_task(last_task, "peaksearch:pf9");
		image->features = search_peaks_peakfinder9(image,
		                                           iargs->peak_search.min_snr_biggest_pix,
		                                           iargs->peak_search.min_snr_peak_pix,
		                                           iargs->peak_search.min_snr,
		                                           iargs->peak_search.min_sig,
		                                           iargs->peak_search.min_peak_over_neighbour,
		                                           iargs->peak_search.local_bg_radius);
		break;

		case PEAK_MSGPACK:
		peaks = image_msgpack_read_peaks(iargs->dtempl,
		                                 pargs->zmq_data,
		                                 pargs->zmq_data_size,
		                                 iargs->peak_search.half_pixel_shift);
		if ( iargs->peak_search.revalidate ) {
			ImageFeatureList *npeaks = validate_peaks(image, peaks,
			                                          iargs->peak_search.min_snr,
			                                          iargs->peak_search.pk_inn,
			                                          iargs->peak_search.pk_mid,
			                                          iargs->peak_search.pk_out,
			                                          iargs->peak_search.use_saturated,
			                                          iargs->peak_search.check_hdf5_snr);
			image_feature_list_free(peaks);
			image->features = npeaks;
		}
		break;

		case PEAK_NONE:
		case PEAK_ERROR:
		break;

	}
	if ( image->features == NULL ) {
		ERROR("Peaksearch failed for image %s" "(event %s).\n",
		      image->filename, image->ev);
	}
	profile_end("peak-search");

	image->peak_resolution = estimate_peak_resolution(image->features,
	                                                  image->lambda,
	                                                  image->detgeom);

	if ( prefilter != NULL ) {
		profile_start("restore-filter-backup");
		restore_image_data(image->dp, image->detgeom, prefilter);
		profile_end("restore-filter-backup");
	}

	rn = getcwd(NULL, 0);

	r = chdir(tmpdir);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		return;
	}

	/* Set beam parameters */
	if ( iargs->fix_divergence >= 0.0 ) {
		image->div = iargs->fix_divergence;
	} else {
		image->div = 0.0;
	}

	if ( image_feature_count(image->features) < iargs->min_peaks ) {
		r = chdir(rn);
		if ( r ) {
			ERROR("Failed to chdir: %s\n", strerror(errno));
			return;
		}
		free(rn);
		image->hit = 0;

		if ( iargs->stream_nonhits ) {
			goto streamwrite;
		} else {
			goto out;
		}
	}
	image->hit = 1;

	/* Index the pattern */
	set_last_task(last_task, "indexing");
	profile_start("index");
	index_pattern_4(image, iargs->ipriv, &sb_shared->pings[cookie],
	                last_task, mille, iargs->max_mille_level);
	profile_end("index");

	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		return;
	}
	free(rn);

	/* Set beam/crystal parameters */
	set_last_task(last_task, "prediction params");
	if ( iargs->fix_profile_r >= 0.0 ) {
		for ( i=0; i<image->n_crystals; i++ ) {
			crystal_set_profile_radius(image->crystals[i].cr,
			                           iargs->fix_profile_r);
			crystal_set_mosaicity(image->crystals[i].cr, 0.0);
		}
	} else {
		for ( i=0; i<image->n_crystals; i++ ) {
			crystal_set_profile_radius(image->crystals[i].cr, 0.02e9);
			crystal_set_mosaicity(image->crystals[i].cr, 0.0);
		}
	}

	if ( iargs->fix_profile_r < 0.0 ) {
		for ( i=0; i<image->n_crystals; i++ ) {
			if ( refine_radius(image->crystals[i].cr, image) ) {
				ERROR("WARNING: Radius determination failed\n");
			}
		}
	}

	/* Integrate! */
	if ( !iargs->cell_params_only ) {
		set_last_task(last_task, "integration");
		profile_start("integration");
		sb_shared->pings[cookie]++;
		integrate_all_5(image, iargs->int_meth, PMODEL_XSPHERE,
		                iargs->push_res,
		                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
		                iargs->int_diag, iargs->int_diag_h,
		                iargs->int_diag_k, iargs->int_diag_l,
		                &sb_shared->term_lock, iargs->overpredict);
		profile_end("integration");
	}

streamwrite:
	set_last_task(last_task, "stream write");
	profile_start("write-stream");
	sb_shared->pings[cookie]++;
	ret = stream_write_chunk(st, image, iargs->stream_flags);
	if ( ret != 0 ) {
		ERROR("Error writing stream file.\n");
	}
	profile_end("write-stream");

	int n = 0;
	for ( i=0; i<image->n_crystals; i++ ) {
		n += crystal_get_num_implausible_reflections(image->crystals[i].cr);
	}
	if ( n > 0 ) {
		STATUS("WARNING: %i implausibly negative reflection%s in %s "
		       "%s\n", n, n>1?"s":"", image->filename, image->ev);
	}

	im_asapo_send(asapostuff, image, image->hit);

out:
	/* Count crystals which are still good */
	set_last_task(last_task, "process_image finalisation");
	sb_shared->pings[cookie]++;
	pthread_mutex_lock(&sb_shared->totals_lock);
	any_crystals = 0;
	for ( i=0; i<image->n_crystals; i++ ) {
		if ( crystal_get_user_flag(image->crystals[i].cr) == 0 ) {
			sb_shared->n_crystals++;
			any_crystals = 1;
		}
	}
	sb_shared->n_processed++;
	sb_shared->n_hits += image->hit;
	sb_shared->n_hadcrystals += any_crystals;
	pthread_mutex_unlock(&sb_shared->totals_lock);

	/* Free image (including detgeom) */
	profile_start("free-image");
	image_free(image);
	profile_end("free-image");

	set_last_task(last_task, "sandbox");
}
