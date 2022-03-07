/*
 * process_image.c
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
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
#include <time-accounts.h>

#include "process_image.h"
#include "predict-refine.h"
#include "im-sandbox.h"
#include "im-zmq.h"

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
                                         TimeAccounts *taccs,
                                         char *last_task,
                                         signed int wait_for_file,
                                         int cookie,
                                         int no_image_data,
                                         int no_mask_data)
{
	signed int file_wait_time = wait_for_file;
	int wait_message_done = 0;
	int read_retry_done = 0;
	int r;
	struct image *image;

	time_accounts_set(taccs, TACC_WAITFILE);
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

		time_accounts_set(taccs, TACC_IMAGE_DATA);
		set_last_task(last_task, "read file");
		sb_shared->pings[cookie]++;

		image = image_read_with_time_accounting(dtempl, filename, event,
		                                        no_image_data, no_mask_data,
		                                        taccs);
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
                   int serial, struct sb_shm *sb_shared, TimeAccounts *taccs,
                   char *last_task)
{
	struct image *image;
	int i;
	int r;
	int ret;
	char *rn;
	float **prefilter;
	int any_crystals;

	if ( pargs->zmq_data != NULL ) {
		time_accounts_set(taccs, TACC_IMAGE_DATA);
		set_last_task(last_task, "unpacking messagepack object");
		image = image_read_data_block(iargs->dtempl,
		                              pargs->zmq_data,
		                              pargs->zmq_data_size,
		                              iargs->data_format,
		                              serial,
		                              iargs->no_image_data,
		                              iargs->no_mask_data);
		if ( image == NULL ) return;
	} else {
		image = file_wait_open_read(pargs->filename, pargs->event,
		                            iargs->dtempl,
		                            sb_shared, taccs, last_task,
		                            iargs->wait_for_file,
		                            cookie,
		                            iargs->no_image_data,
		                            iargs->no_mask_data);
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
	time_accounts_set(taccs, TACC_FILTER);
	set_last_task(last_task, "image filter");
	sb_shared->pings[cookie]++;
	prefilter = backup_image_data(image->dp, image->detgeom);

	if ( iargs->median_filter > 0 ) {
		filter_median(image, iargs->median_filter);
	}

	if ( iargs->noisefilter ) {
		filter_noise(image);
	}

	time_accounts_set(taccs, TACC_RESRANGE);
	set_last_task(last_task, "resolution range");
	sb_shared->pings[cookie]++;
	mark_resolution_range_as_bad(image, iargs->highres, +INFINITY);

	time_accounts_set(taccs, TACC_PEAKSEARCH);
	sb_shared->pings[cookie]++;
	switch ( iargs->peaks ) {

		case PEAK_HDF5:
		case PEAK_CXI:
		set_last_task(last_task, "peaksearch:hdf5orcxi");
		image->features = image_read_peaks(iargs->dtempl,
		                                   pargs->filename,
		                                   pargs->event,
		                                   iargs->half_pixel_shift);
		if ( image->features == NULL ) {
			ERROR("Failed to get peaks from HDF5 file.\n");
		}
		if ( !iargs->no_revalidate ) {
			validate_peaks(image, iargs->min_snr,
				       iargs->pk_inn, iargs->pk_mid,
		                       iargs->pk_out, iargs->use_saturated,
				       iargs->check_hdf5_snr);
		}
		break;

		case PEAK_ZAEF:
		set_last_task(last_task, "peaksearch:zaef");
		search_peaks(image, iargs->threshold,
		             iargs->min_sq_gradient, iargs->min_snr,
		             iargs->pk_inn, iargs->pk_mid, iargs->pk_out,
		             iargs->use_saturated);
		break;

		case PEAK_PEAKFINDER8:
		set_last_task(last_task, "peaksearch:pf8");
		if ( search_peaks_peakfinder8(image, 2048,
		                              iargs->threshold,
		                              iargs->min_snr,
		                              iargs->min_pix_count,
		                              iargs->max_pix_count,
		                              iargs->local_bg_radius,
		                              iargs->min_res,
		                              iargs->max_res,
		                              iargs->use_saturated) ) {
			ERROR("Failed to find peaks in image %s"
			      "(event %s).\n",
			      image->filename, image->ev);
		}
		break;

		case PEAK_PEAKFINDER9:
		set_last_task(last_task, "peaksearch:pf9");
		if ( search_peaks_peakfinder9(image,
		                              iargs->min_snr_biggest_pix,
		                              iargs->min_snr_peak_pix,
		                              iargs->min_snr,
		                              iargs->min_sig,
		                              iargs->min_peak_over_neighbour,
		                              iargs->local_bg_radius) )
		{
			ERROR("Failed to find peaks in image %s"
			      "(event %s).\n",
			      image->filename, image->ev);
		}
		break;

		case PEAK_MSGPACK:
		image->features = image_msgpack_read_peaks(iargs->dtempl,
		                                           pargs->zmq_data,
		                                           pargs->zmq_data_size,
		                                           iargs->half_pixel_shift);
		break;

		case PEAK_NONE:
		case PEAK_ERROR:
		break;

	}

	image->peak_resolution = estimate_peak_resolution(image->features,
	                                                  image->lambda,
	                                                  image->detgeom);

	restore_image_data(image->dp, image->detgeom, prefilter);

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
	time_accounts_set(taccs, TACC_INDEXING);
	set_last_task(last_task, "indexing");
	index_pattern_3(image, iargs->ipriv, &sb_shared->pings[cookie],
	                last_task);

	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		return;
	}
	free(rn);

	/* Set beam/crystal parameters */
	time_accounts_set(taccs, TACC_PREDPARAMS);
	set_last_task(last_task, "prediction params");
	if ( iargs->fix_profile_r >= 0.0 ) {
		for ( i=0; i<image->n_crystals; i++ ) {
			crystal_set_profile_radius(image->crystals[i],
			                           iargs->fix_profile_r);
			crystal_set_mosaicity(image->crystals[i], 0.0);
		}
	} else {
		for ( i=0; i<image->n_crystals; i++ ) {
			crystal_set_profile_radius(image->crystals[i], 0.02e9);
			crystal_set_mosaicity(image->crystals[i], 0.0);
		}
	}

	if ( iargs->fix_profile_r < 0.0 ) {
		for ( i=0; i<image->n_crystals; i++ ) {
			if ( refine_radius(image->crystals[i], image) ) {
				ERROR("WARNING: Radius determination failed\n");
			}
		}
	}

	/* Integrate! */
	time_accounts_set(taccs, TACC_INTEGRATION);
	set_last_task(last_task, "integration");
	sb_shared->pings[cookie]++;
	integrate_all_5(image, iargs->int_meth, PMODEL_XSPHERE,
	                iargs->push_res,
	                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
	                iargs->int_diag, iargs->int_diag_h,
	                iargs->int_diag_k, iargs->int_diag_l,
	                &sb_shared->term_lock, iargs->overpredict);

streamwrite:
	time_accounts_set(taccs, TACC_WRITESTREAM);
	set_last_task(last_task, "stream write");
	sb_shared->pings[cookie]++;
	ret = stream_write_chunk(st, image, iargs->stream_flags);
	if ( ret != 0 ) {
		ERROR("Error writing stream file.\n");
	}

	int n = 0;
	for ( i=0; i<image->n_crystals; i++ ) {
		n += crystal_get_num_implausible_reflections(image->crystals[i]);
	}
	if ( n > 0 ) {
		STATUS("WARNING: %i implausibly negative reflection%s in %s "
		       "%s\n", n, n>1?"s":"", image->filename, image->ev);
	}

out:
	/* Count crystals which are still good */
	time_accounts_set(taccs, TACC_TOTALS);
	set_last_task(last_task, "process_image finalisation");
	sb_shared->pings[cookie]++;
	pthread_mutex_lock(&sb_shared->totals_lock);
	any_crystals = 0;
	for ( i=0; i<image->n_crystals; i++ ) {
		if ( crystal_get_user_flag(image->crystals[i]) == 0 ) {
			sb_shared->n_crystals++;
			any_crystals = 1;
		}
	}
	sb_shared->n_processed++;
	sb_shared->n_hits += image->hit;
	sb_shared->n_hadcrystals += any_crystals;
	pthread_mutex_unlock(&sb_shared->totals_lock);

	/* Free image (including detgeom) */
	image_free(image);

	set_last_task(last_task, "sandbox");
}
