/*
 * process_image.c
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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


void process_image(const struct index_args *iargs, struct pattern_args *pargs,
                   Stream *st, int cookie)
{
	float *data_for_measurement;
	size_t data_size;
	int check;
	struct hdfile *hdfile;
	struct image image;
	int i;
	char filename[1024];

	/* Prefix to jump out of temporary folder */
	if ( pargs->filename[0] != '/' ) {
		snprintf(filename, 1023, "../../%s", pargs->filename);
	} else {
		snprintf(filename, 1023, "%s", pargs->filename);
	}

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.copyme = iargs->copyme;
	image.id = cookie;
	image.filename = pargs->filename;  /* Relative to top level */
	image.beam = iargs->beam;
	image.det = iargs->det;
	image.crystals = NULL;
	image.n_crystals = 0;

	hdfile = hdfile_open(filename);  /* Relative to temporary folder */
	if ( hdfile == NULL ) return;

	if ( iargs->element != NULL ) {

		int r;
		r = hdfile_set_image(hdfile, iargs->element);
		if ( r ) {
			ERROR("Couldn't select path '%s'\n", iargs->element);
			hdfile_close(hdfile);
			return;
		}

	} else {

		int r;
		r = hdfile_set_first_image(hdfile, "/");
		if ( r ) {
			ERROR("Couldn't select first path\n");
			hdfile_close(hdfile);
			return;
		}

	}

	check = hdf5_read(hdfile, &image, 1);
	if ( check ) {
		hdfile_close(hdfile);
		return;
	}

	if ( (image.width != image.det->max_fs + 1 )
	  || (image.height != image.det->max_ss + 1))
	{
		ERROR("Image size doesn't match geometry size"
			" - rejecting image.\n");
		ERROR("Image size: %i,%i.  Geometry size: %i,%i\n",
		      image.width, image.height,
		      image.det->max_fs + 1, image.det->max_ss + 1);
		hdfile_close(hdfile);
		return;
	}

	fill_in_values(image.det, hdfile);
	fill_in_beam_parameters(image.beam, hdfile);

	image.lambda = ph_en_to_lambda(eV_to_J(image.beam->photon_energy));

	if ( (image.beam->photon_energy < 0.0) || (image.lambda > 1000) ) {
		/* Error message covers a silly value in the beam file or in
		 * the HDF5 file. */
		ERROR("Nonsensical wavelength (%e m or %e eV) value for %s.\n",
		      image.lambda, image.beam->photon_energy, image.filename);
		hdfile_close(hdfile);
		return;
	}

	if ( iargs->cmfilter ) filter_cm(&image);

	/* Take snapshot of image after CM subtraction but before
	 * the aggressive noise filter. */
	data_size = image.width * image.height * sizeof(float);
	data_for_measurement = malloc(data_size);

	if ( iargs->noisefilter ) {
		filter_noise(&image, data_for_measurement);
	} else {
		memcpy(data_for_measurement, image.data, data_size);
	}

	switch ( iargs->peaks ) {

		case PEAK_HDF5:
		// Get peaks from HDF5
		if (get_peaks(&image, hdfile,
			iargs->hdf5_peak_path)) {
			ERROR("Failed to get peaks from HDF5 file.\n");
		}
		if ( !iargs->no_revalidate ) {
			validate_peaks(&image, iargs->min_int_snr,
				       iargs->ir_inn, iargs->ir_mid,
				       iargs->ir_out, iargs->use_saturated);
		}
		break;

		case PEAK_ZAEF:
		search_peaks(&image, iargs->threshold,
		             iargs->min_gradient, iargs->min_snr,
		             iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
		             iargs->use_saturated);
		break;

	}

	/* Get rid of noise-filtered version at this point
	 * - it was strictly for the purposes of peak detection. */
	free(image.data);
	image.data = data_for_measurement;

	/* Index the pattern */
	index_pattern(&image, iargs->indm, iargs->ipriv);

	pargs->n_crystals = image.n_crystals;

	/* Default beam parameters */
	image.div = image.beam->divergence;
	image.bw = image.beam->bandwidth;

	/* Integrate each crystal's diffraction spots */
	for ( i=0; i<image.n_crystals; i++ ) {

		RefList *reflections;

		/* Set default crystal parameter(s) */
		crystal_set_profile_radius(image.crystals[i],
		                           image.beam->profile_radius);

		if ( iargs->integrate_found ) {
			reflections = select_intersections(&image,
			                                   image.crystals[i]);
		} else {
			reflections = find_intersections(&image,
			                                 image.crystals[i]);
		}

		crystal_set_reflections(image.crystals[i], reflections);

	}

	/* Integrate all the crystals at once - need all the crystals so that
	 * overlaps can be detected. */
	integrate_reflections(&image, iargs->closer,
	                              iargs->bgsub,
	                              iargs->min_int_snr,
	                              iargs->ir_inn,
	                              iargs->ir_mid,
	                              iargs->ir_out,
	                              iargs->integrate_saturated,
	                              iargs->res_cutoff);

	write_chunk(st, &image, hdfile,
	            iargs->stream_peaks, iargs->stream_refls);

	for ( i=0; i<image.n_crystals; i++ ) {
		crystal_free(image.crystals[i]);
	}

	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	image_feature_list_free(image.features);
	hdfile_close(hdfile);
}
