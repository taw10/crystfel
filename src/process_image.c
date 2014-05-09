/*
 * process_image.c
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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
#include <unistd.h>

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "thread-pool.h"
#include "beam-parameters.h"
#include "geometry.h"
#include "stream.h"
#include "reflist-utils.h"
#include "process_image.h"
#include "integration.h"


void process_image(const struct index_args *iargs, struct pattern_args *pargs,
                   Stream *st, int cookie, const char *tmpdir, int results_pipe)
{
	float *data_for_measurement;
	size_t data_size;
	int check;
	struct hdfile *hdfile;
	struct image image;
	int i;
	int r;
	char *rn;

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.copyme = iargs->copyme;
	image.id = cookie;
	image.filename = pargs->filename_p_e->filename;
	image.event = pargs->filename_p_e->ev;
	image.beam = iargs->beam;
	image.det = iargs->det;
	image.crystals = NULL;
	image.n_crystals = 0;

	hdfile = hdfile_open(image.filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open file: %s\n", image.filename);
		return;
	}

	check = hdf5_read2(hdfile, &image, image.event, 0);
	if ( check ) {
		return;
	}

	/* Take snapshot of image after CM subtraction but before applying
	 * horrible noise filters to it */
	data_size = image.width * image.height * sizeof(float);
	data_for_measurement = malloc(data_size);
	memcpy(data_for_measurement, image.data, data_size);

	if ( iargs->median_filter > 0 ) {
		filter_median(&image, iargs->median_filter);
	}

	if ( iargs->noisefilter ) {
		filter_noise(&image);
	}

	mark_resolution_range_as_bad(&image, iargs->highres, +INFINITY);

	hdfile = hdfile_open(image.filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open file %s.\n", image.filename);
		return;
	}

	switch ( iargs->peaks ) {

		case PEAK_HDF5:
		/* Get peaks from HDF5 */

		if ( !single_panel_data_source(iargs->det, iargs->element) ) {
			ERROR("Peaks from HDF5 file not supported with multiple panel data sources.\n");
		}

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

		case PEAK_ZAEF:
		search_peaks(&image, iargs->threshold,
		             iargs->min_gradient, iargs->min_snr,
		             iargs->pk_inn, iargs->pk_mid,iargs->pk_out,
		             iargs->use_saturated);
		break;

	}

	/* Get rid of noise-filtered version at this point
	 * - it was strictly for the purposes of peak detection. */
	free(image.data);
	image.data = data_for_measurement;

	rn = getcwd(NULL, 0);

	r = chdir(tmpdir);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		hdfile_close(hdfile);
		return;
	}

	/* Index the pattern */
	index_pattern(&image, iargs->indm, iargs->ipriv);

	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		hdfile_close(hdfile);
		return;
	}

	pargs->n_crystals = image.n_crystals;
	for ( i=0; i<image.n_crystals; i++ ) {
		crystal_set_image(image.crystals[i], &image);
	}

	/* Default parameters */
	image.div = image.beam->divergence;
	image.bw = image.beam->bandwidth;
	for ( i=0; i<image.n_crystals; i++ ) {
		crystal_set_profile_radius(image.crystals[i],
		                           image.beam->profile_radius);
		crystal_set_mosaicity(image.crystals[i], 0.0);  /* radians */
	}

	/* Integrate all the crystals at once - need all the crystals so that
	 * overlaps can be detected. */
	integrate_all_4(&image, iargs->int_meth, PMODEL_SPHERE, iargs->push_res,
	                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
	                iargs->int_diag, iargs->int_diag_h,
	                iargs->int_diag_k, iargs->int_diag_l, results_pipe);

	write_chunk(st, &image, hdfile,
	            iargs->stream_peaks, iargs->stream_refls,
	            pargs->filename_p_e->ev);

	for ( i=0; i<image.n_crystals; i++ ) {
		cell_free(crystal_get_cell(image.crystals[i]));
		reflist_free(crystal_get_reflections(image.crystals[i]));
		crystal_free(image.crystals[i]);
	}
	free(image.crystals);

	for ( i=0; i<image.det->n_panels; i++ ) {
		free(image.dp[i]);
		free(image.bad[i]);
	}
	free(image.dp);
	free(image.bad);

	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	image_feature_list_free(image.features);
	hdfile_close(hdfile);
}
