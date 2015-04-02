/*
 * process_image.h
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

#ifndef PROCESS_IMAGE_H
#define PROCESS_IMAGE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "integration.h"


enum {
	PEAK_ZAEF,
	PEAK_HDF5,
	PEAK_CXI,
};


/* Information about the indexing process which is common to all patterns */
struct index_args
{
	UnitCell *cell;
	int cmfilter;
	int noisefilter;
	int median_filter;
	int satcorr;
	float threshold;
	float min_gradient;
	float min_snr;
	int check_hdf5_snr;
	struct detector *det;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	int peaks;                /* Peak detection method */
	float tols[4];
	struct beam_params *beam;
	char *hdf5_peak_path;
	float pk_inn;
	float pk_mid;
	float pk_out;
	float ir_inn;
	float ir_mid;
	float ir_out;
	struct copy_hdf5_field *copyme;
	int integrate_saturated;
	int use_saturated;
	int no_revalidate;
	int stream_peaks;
	int stream_refls;
	IntegrationMethod int_meth;
	IntDiag int_diag;
	signed int int_diag_h;
	signed int int_diag_k;
	signed int int_diag_l;
	float push_res;
	float highres;
	float fix_profile_r;
	float fix_bandwidth;
	float fix_divergence;
	int predict_refine;
};


/* Information about the indexing process for one pattern */
struct pattern_args
{
	/* "Input" */
	struct filename_plus_event *filename_p_e;

	/* "Output" */
	int n_crystals;
};


extern void process_image(const struct index_args *iargs,
                          struct pattern_args *pargs, Stream *st,
                          int cookie, const char *tmpdir, int results_pipe,
                          int serial);


#endif	/* PROCESS_IMAGEs_H */
