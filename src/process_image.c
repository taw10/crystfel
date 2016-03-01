/*
 * process_image.c
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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
#include <hdf5.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <unistd.h>

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "thread-pool.h"
#include "geometry.h"
#include "stream.h"
#include "reflist-utils.h"
#include "process_image.h"
#include "integration.h"
#include "predict-refine.h"
#include "im-sandbox.h"
#include "time-accounts.h"


static float **backup_image_data(float **dp, struct detector *det)
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


static void restore_image_data(float **dp, struct detector *det, float **bu)
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


void process_image(const struct index_args *iargs, struct pattern_args *pargs,
                   Stream *st, int cookie, const char *tmpdir,
                   int serial, struct sb_shm *sb_shared, TimeAccounts *taccs)
{
	int check;
	struct hdfile *hdfile;
	struct image image;
	int i;
	int r;
	int ret;
	char *rn;
	float **prefilter;
	int any_crystals;

	image.features = NULL;
	image.copyme = iargs->copyme;
	image.id = cookie;
	image.filename = pargs->filename_p_e->filename;
	image.event = pargs->filename_p_e->ev;
	image.beam = iargs->beam;
	image.det = copy_geom(iargs->det);
	image.crystals = NULL;
	image.n_crystals = 0;
	image.serial = serial;
	image.indexed_by = INDEXING_NONE;

	time_accounts_set(taccs, TACC_HDF5OPEN);
	hdfile = hdfile_open(image.filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open file: %s\n", image.filename);
		return;
	}

	time_accounts_set(taccs, TACC_HDF5READ);
	check = hdf5_read2(hdfile, &image, image.event, 0);
	if ( check ) {
		return;
	}

	/* Take snapshot of image before applying horrible noise filters */
	time_accounts_set(taccs, TACC_FILTER);
	prefilter = backup_image_data(image.dp, image.det);

	if ( iargs->median_filter > 0 ) {
		filter_median(&image, iargs->median_filter);
	}

	if ( iargs->noisefilter ) {
		filter_noise(&image);
	}

	time_accounts_set(taccs, TACC_RESRANGE);
	mark_resolution_range_as_bad(&image, iargs->highres, +INFINITY);

	time_accounts_set(taccs, TACC_PEAKSEARCH);
	switch ( iargs->peaks ) {

		case PEAK_HDF5:
		if ( get_peaks(&image, hdfile, iargs->hdf5_peak_path) ) {
			ERROR("Failed to get peaks from HDF5 file.\n");
		}
		if ( !iargs->no_revalidate ) {
			validate_peaks(&image, iargs->min_snr,
				       iargs->pk_inn, iargs->pk_mid,
		                       iargs->pk_out, iargs->use_saturated,
				       iargs->check_hdf5_snr);
		}
		break;

		case PEAK_CXI:
		if ( get_peaks_cxi(&image, hdfile, iargs->hdf5_peak_path,
		                   pargs->filename_p_e) ) {
			ERROR("Failed to get peaks from CXI file.\n");
		}
		if ( !iargs->no_revalidate ) {
			validate_peaks(&image, iargs->min_snr,
				       iargs->pk_inn, iargs->pk_mid,
		                       iargs->pk_out, iargs->use_saturated,
				       iargs->check_hdf5_snr);
		}
		break;

		case PEAK_ZAEF:
		search_peaks(&image, iargs->threshold,
		             iargs->min_gradient, iargs->min_snr,
		             iargs->pk_inn, iargs->pk_mid,iargs->pk_out,
		             iargs->use_saturated);
		break;

	}

	restore_image_data(image.dp, image.det, prefilter);

	rn = getcwd(NULL, 0);

	r = chdir(tmpdir);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		hdfile_close(hdfile);
		return;
	}

	/* Set beam parameters */
	if ( iargs->fix_divergence >= 0.0 ) {
		image.div = iargs->fix_divergence;
	} else {
		image.div = 0.0;
	}
	if ( iargs->fix_bandwidth >= 0.0 ) {
		image.bw = iargs->fix_bandwidth;
	} else {
		image.bw = 0.00000001;
	}

	/* Index the pattern */
	time_accounts_set(taccs, TACC_INDEXING);
	index_pattern(&image, iargs->indm, iargs->ipriv);

	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		hdfile_close(hdfile);
		return;
	}
	free(rn);

	/* Set beam/crystal parameters */
	time_accounts_set(taccs, TACC_PREDPARAMS);
	if ( iargs->fix_profile_r >= 0.0 ) {
		for ( i=0; i<image.n_crystals; i++ ) {
			crystal_set_profile_radius(image.crystals[i],
			                           iargs->fix_profile_r);
			crystal_set_mosaicity(image.crystals[i], 0.0);
		}
	} else {
		for ( i=0; i<image.n_crystals; i++ ) {
			crystal_set_profile_radius(image.crystals[i], 0.02e9);
			crystal_set_mosaicity(image.crystals[i], 0.0);
		}
	}

	if ( iargs->fix_profile_r < 0.0 ) {
		for ( i=0; i<image.n_crystals; i++ ) {
			refine_radius(image.crystals[i], &image);
		}
	}

	/* Integrate! */
	time_accounts_set(taccs, TACC_INTEGRATION);
	integrate_all_4(&image, iargs->int_meth, PMODEL_SCSPHERE,
	                iargs->push_res,
	                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
	                iargs->int_diag, iargs->int_diag_h,
	                iargs->int_diag_k, iargs->int_diag_l,
	                &sb_shared->term_lock);

	time_accounts_set(taccs, TACC_WRITESTREAM);
	ret = write_chunk(st, &image, hdfile,
	                  iargs->stream_peaks, iargs->stream_refls,
	                  pargs->filename_p_e->ev);
	if ( ret != 0 ) {
		ERROR("Error writing stream file.\n");
	}

	int n = 0;
	for ( i=0; i<image.n_crystals; i++ ) {
		n += crystal_get_num_implausible_reflections(image.crystals[i]);
	}
	if ( n > 0 ) {
		STATUS("WARNING: %i implausibly negative reflection%s in %s "
		       "%s\n", n, n>1?"s":"", image.filename,
		       get_event_string(image.event));
	}

	/* Count crystals which are still good */
	time_accounts_set(taccs, TACC_TOTALS);
	pthread_mutex_lock(&sb_shared->totals_lock);
	any_crystals = 0;
	for ( i=0; i<image.n_crystals; i++ ) {
		if ( crystal_get_user_flag(image.crystals[i]) == 0 ) {
			sb_shared->n_crystals++;
			any_crystals = 1;
		}
	}
	sb_shared->n_processed++;
	sb_shared->n_hadcrystals += any_crystals;
	pthread_mutex_unlock(&sb_shared->totals_lock);

	for ( i=0; i<image.n_crystals; i++ ) {
		cell_free(crystal_get_cell(image.crystals[i]));
		reflist_free(crystal_get_reflections(image.crystals[i]));
		crystal_free(image.crystals[i]);
	}
	free(image.crystals);

	for ( i=0; i<image.det->n_panels; i++ ) {
		free(image.dp[i]);
		free(image.bad[i]);
		free(image.sat[i]);
	}
	free(image.dp);
	free(image.bad);
	free(image.sat);

	image_feature_list_free(image.features);
	free_detector_geometry(image.det);
	hdfile_close(hdfile);
}
