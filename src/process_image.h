/*
 * process_image.h
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2016 Thomas White <taw@physics.org>
 *   2014-2017 Valerio Mariani <valerio.mariani@desy.de>
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

struct index_args;

#include "integration.h"
#include "im-sandbox.h"
#include "time-accounts.h"
#include "taketwo.h"


enum {
        PEAK_PEAKFINDER8,
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
	IndexingPrivate *ipriv;
	int peaks;                /* Peak detection method */
	float tols[4];
	struct beam_params *beam;
	char *hdf5_peak_path;
	int half_pixel_shift;
	float pk_inn;
	float pk_mid;
	float pk_out;
	float ir_inn;
	float ir_mid;
	float ir_out;
	int min_res;
	int max_res;
	int max_n_peaks;
	int min_pix_count;
	int max_pix_count;
	int local_bg_radius;
	int min_peaks;
	struct imagefile_field_list *copyme;
	int integrate_saturated;
	int use_saturated;
	int no_revalidate;
	int stream_peaks;
	int stream_refls;
	int stream_nonhits;
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
	char *felix_options;
	int profile;  /* Whether or not to do wall clock profiling */
	struct taketwo_options taketwo_opts;
};


/* Information about the indexing process for one pattern */
struct pattern_args
{
	/* "Input" */
	struct filename_plus_event *filename_p_e;
};


extern void process_image(const struct index_args *iargs,
                          struct pattern_args *pargs, Stream *st,
                          int cookie, const char *tmpdir, int serial,
                          struct sb_shm *sb_shared, TimeAccounts *taccs);


#endif	/* PROCESS_IMAGE_H */
