/*
 * indexamajig.c
 *
 * Index patterns, output hkl+intensity etc.
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
 *   2017      Valerio Mariani <valerio.mariani@desy.de>
 *   2017-2018 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <argp.h>
#include <gsl/gsl_errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <utils.h>
#include <index.h>
#include <peaks.h>
#include <filters.h>
#include <thread-pool.h>
#include <geometry.h>
#include <stream.h>
#include <reflist-utils.h>
#include <cell-utils.h>
#include <integration.h>
#include <image.h>
#include <datatemplate.h>

#include "im-sandbox.h"
#include "version.h"
#include "json-utils.h"


struct indexamajig_arguments
{
	struct index_args iargs;  /* These are the options that will be
	                           * given to process_image */
	char *filename;
	char *geom_filename;
	char *outfile;
	char *prefix;
	int check_prefix;
	int n_proc;
	char *cellfile;
	char *indm_str;
	int basename;
	char *zmq_addr;
	int serial_start;
	char *temp_location;
	int if_refine;
	int if_checkcell;
	int if_peaks;
	int if_multi;
	int if_retry;
	int profile;  /* Whether to do wall-clock time profiling */
	char **copy_headers;
	int n_copy_headers;
	char *harvest_file;

	TakeTwoOptions **taketwo_opts_ptr;
	FelixOptions **felix_opts_ptr;
	XGandalfOptions **xgandalf_opts_ptr;
	PinkIndexerOptions **pinkindexer_opts_ptr;
	FromFileOptions **fromfile_opts_ptr;
};


static double nan_if_neg(double n)
{
	if ( n < 0.0 ) {
		return NAN;
	} else {
		return n;
	}
}


static void write_json_cell(FILE *fh, const char *name, UnitCell *cell)
{
	double a, b, c, al, be, ga;
	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	if ( cell == NULL ) {
		fprintf(fh, "    \"%s\": null,\n", name);
	} else {
		fprintf(fh, "    \"%s\": {\n", name);
		fprintf(fh, "      \"lattice_type\": \"%s\",\n",
		        str_lattice(cell_get_lattice_type(cell)));
		fprintf(fh, "      \"centering\": \"%c\",\n",
		        cell_get_centering(cell));
		fprintf(fh, "      \"unique_axis\": \"%c\",\n",
		        cell_get_unique_axis(cell));
		fprintf(fh, "      \"a_Angstrom\": %f,\n", a*1e10);
		fprintf(fh, "      \"b_Angstrom\": %f,\n", b*1e10);
		fprintf(fh, "      \"c_Angstrom\": %f,\n", c*1e10);
		fprintf(fh, "      \"alpha_deg\": %f,\n", rad2deg(al));
		fprintf(fh, "      \"beta_deg\": %f,\n", rad2deg(be));
		fprintf(fh, "      \"gamma_deg\": %f\n", rad2deg(ga));  /* NB No comma */
		fprintf(fh, "    },\n");
	}
}


static void write_methods(FILE *fh, const char *name, IndexingPrivate *ipriv)
{
	int i, n;
	const IndexingMethod *methods;

	fprintf(fh, "    \"%s\": [", name);

	if ( ipriv != NULL ) {
		methods = indexing_methods(ipriv, &n);
		for ( i=0; i<n; i++ ) {
			fprintf(fh, " \"%s\"", indexer_str(methods[i]));
			if ( i < n-1 ) {
				fprintf(fh, ",");
			}
		}
	}
	fprintf(fh, " ],\n");
}


static void write_harvest_file(struct index_args *args,
                               const char *filename,
                               int if_multi, int if_refine, int if_retry,
                               int if_peaks, int if_checkcell)
{
	FILE *fh;
	char *tmp;

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Unable to write parameter harvesting file.\n");
		return;
	}

	fprintf(fh, "{\n");
	fprintf(fh, "  \"input\": {\n");
	write_float(fh, 0, "highres", args->highres);
	fprintf(fh, "  },\n");

	fprintf(fh, "  \"peaksearch\": {\n");
	write_str(fh, 1, "method", str_peaksearch(args->peaks));
	write_float(fh, 1, "radius_inner_px", args->pk_inn);
	write_float(fh, 1, "radius_middle_px", args->pk_mid);
	write_float(fh, 1, "radius_outer_px", args->pk_out);
	write_bool(fh, 1, "noise_filter", args->noisefilter);
	write_int(fh, 1, "median_filter", args->median_filter);
	write_float(fh, 1, "threshold_adu", args->threshold);
	write_float(fh, 1, "min_squared_gradient_adu2", args->min_sq_gradient);
	write_float(fh, 1, "min_snr", args->min_snr);
	write_bool(fh, 1, "check_hdf5_snr", args->check_hdf5_snr);
	write_bool(fh, 1, "half_pixel_shift", args->half_pixel_shift);
	write_int(fh, 1, "min_res_px", args->min_res);
	write_int(fh, 1, "max_res_px", args->max_res);
	write_int(fh, 1, "min_pixel_count", args->min_pix_count);
	write_int(fh, 1, "max_pixel_count", args->max_pix_count);
	write_int(fh, 1, "local_bg_radius_px", args->local_bg_radius);
	write_bool(fh, 1, "use_saturated", args->use_saturated);
	write_bool(fh, 1, "revalidate_hdf5", 1-args->no_revalidate);
	write_float(fh, 1, "min_snr_of_biggest_pixel", args->min_snr_biggest_pix);
	write_float(fh, 1, "min_snr_of_peak_pixel", args->min_snr_peak_pix);
	write_float(fh, 1, "min_sig_adu", args->min_sig);
	write_float(fh, 0, "min_peak_over_neighbour_adu", args->min_peak_over_neighbour);
	fprintf(fh, "  },\n");

	fprintf(fh, "  \"hitfinding\": {\n");
	write_int(fh, 0, "min_num_peaks", args->min_peaks);
	fprintf(fh, "  },\n");

	if ( args->ipriv == NULL ) {
		fprintf(fh, "  \"indexing\": null\n");
		fprintf(fh, "  \"integration\": null\n");
	} else {
		fprintf(fh, "  \"indexing\": {\n");
		write_methods(fh, "methods", args->ipriv);
		write_json_cell(fh, "target_cell", args->cell);
		write_float(fh, 1, "tolerance_a_percent", 100.0*args->tols[0]);
		write_float(fh, 1, "tolerance_b_percent", 100.0*args->tols[1]);
		write_float(fh, 1, "tolerance_c_percent", 100.0*args->tols[2]);
		write_float(fh, 1, "tolerance_alpha_deg", rad2deg(args->tols[3]));
		write_float(fh, 1, "tolerance_beta_deg", rad2deg(args->tols[4]));
		write_float(fh, 1, "tolerance_gamma_deg", rad2deg(args->tols[5]));
		write_bool(fh, 1, "multi_lattice", if_multi);
		write_bool(fh, 1, "refine", if_refine);
		write_bool(fh, 1, "retry", if_retry);
		write_bool(fh, 1, "check_peaks", if_peaks);
		write_bool(fh, 1, "check_cell", if_checkcell);
		write_float(fh, 1, "wavelength_estimate_m", args->wavelength_estimate);
		write_float(fh, 0, "clen_estimate_m", args->clen_estimate);
		fprintf(fh, "  },\n");

		fprintf(fh, "  \"integration\": {\n");
		tmp = str_integration_method(args->int_meth);
		write_str(fh, 1, "method", tmp);
		free(tmp);
		write_float(fh, 1, "radius_inner_px", args->ir_inn);
		write_float(fh, 1, "radius_middle_px", args->ir_mid);
		write_float(fh, 1, "radius_outer_px", args->ir_out);
		write_float(fh, 1, "push_res_invm", args->push_res);
		write_float(fh, 1, "fix_profile_radius_invm", nan_if_neg(args->fix_profile_r));
		write_float(fh, 1, "fix_divergence_rad", nan_if_neg(args->fix_divergence));
		write_bool(fh, 0, "overpredict", args->overpredict);
		fprintf(fh, "  }\n"); /* NB No comma */
	}

	fprintf(fh, "}\n");

	fclose(fh);
}


static void show_version(FILE *fh, struct argp_state *state)
{
	printf("CrystFEL: %s\n", crystfel_version_string());
	printf("%s\n", crystfel_licence_string());
}


static void add_copy_header(struct indexamajig_arguments *args,
                            const char *header)
{
	char **new_copy_headers = realloc(args->copy_headers,
	                                  (args->n_copy_headers+1)*sizeof(char *));
	if ( new_copy_headers == NULL ) {
		ERROR("Failed to add copy header '%s'\n", header);
		return;
	}

	args->copy_headers = new_copy_headers;
	args->copy_headers[args->n_copy_headers++] = strdup(header);
}


static error_t parse_arg(int key, char *arg, struct argp_state *state)
{
	float tmp;
	int r;
	struct indexamajig_arguments *args = state->input;

	switch ( key ) {

		case ARGP_KEY_INIT :
		state->child_inputs[0] = args->taketwo_opts_ptr;
		state->child_inputs[1] = args->felix_opts_ptr;
		state->child_inputs[2] = args->xgandalf_opts_ptr;
		state->child_inputs[3] = args->pinkindexer_opts_ptr;
		state->child_inputs[4] = args->fromfile_opts_ptr;
		break;

		case 'h' :
		argp_state_help(state, stdout, ARGP_HELP_STD_HELP);
		break;  /* argp_state_help doesn't return */

		case 'v' :
		show_version(stdout, state);
		exit(0);

		case 'i' :
		args->filename = strdup(arg);
		break;

		case 'o' :
		args->outfile = strdup(arg);
		break;

		case 'x' :
		args->prefix = strdup(arg);
		break;

		case 'j' :
		args->n_proc = atoi(arg);
		break;

		case 'g' :
		args->geom_filename = arg;
		break;

		case 201 :
		args->basename = 1;
		break;

		case 202 :
		args->check_prefix = 0;
		break;

		case 203 :
		if ( sscanf(arg, "%f", &tmp) != 1 ) {
			ERROR("Invalid value for --highres\n");
			return EINVAL;
		}
		args->iargs.highres = 1.0 / (tmp/1e10); /* A -> m^-1 */
		break;

		case 204 :
		args->profile = 1;
		break;

		case 205 :
		args->temp_location = strdup(arg);
		break;

		case 206 :
		if (sscanf(arg, "%d", &args->iargs.wait_for_file) != 1)
		{
			ERROR("Invalid value for --wait-for-file\n");
			return EINVAL;
		}
		break;

		case 207 :
		args->zmq_addr = strdup(arg);
		break;

		case 208 :
		args->iargs.no_image_data = 1;
		break;

		case 209 :
		ERROR("--spectrum-filename is no longer used.\n");
		return 1;

		case 210 :
		args->iargs.no_mask_data = 1;
		break;

		/* ---------- Peak search ---------- */

		case 't' :
		args->iargs.threshold = strtof(arg, NULL);
		break;

		case 301 :
		args->iargs.peaks = parse_peaksearch(arg);
		if ( args->iargs.peaks == PEAK_ERROR ) {
			ERROR("Unrecognised peak detection method '%s'\n", arg);
			return EINVAL;
		}
		break;

		case 302 :
		r = sscanf(arg, "%f,%f,%f", &args->iargs.pk_inn,
		           &args->iargs.pk_mid, &args->iargs.pk_out);
		if ( (r != 3) || (args->iargs.pk_inn < 0) ) {
			ERROR("Invalid parameters for '--peak-radius'\n");
			return EINVAL;
		}
		break;

		case 303 :
		if (sscanf(arg, "%d", &args->iargs.min_peaks) != 1)
		{
			ERROR("Invalid value for --min-peaks\n");
			return EINVAL;
		}
		break;

		case 304 :
		ERROR("The option --hdf5-peak-path is no longer used.\n");
		ERROR("Set the peak path in the geometry file.\n");
		break;

		case 305 :
		if (sscanf(arg, "%d", &args->iargs.median_filter) != 1)
		{
			ERROR("Invalid value for --median-filter\n");
			return EINVAL;
		}
		break;

		case 306 :
		args->iargs.noisefilter = 1;
		break;

		case 307 :
		if (sscanf(arg, "%f", &args->iargs.min_sq_gradient) != 1)
		{
			ERROR("Invalid value for --min-squared-gradient\n");
			return EINVAL;
		}
		break;

		case 308 :
		if (sscanf(arg, "%f", &args->iargs.min_snr) != 1)
		{
			ERROR("Invalid value for --min-snr\n");
			return EINVAL;
		}
		break;

		case 309 :
		if (sscanf(arg, "%d", &args->iargs.min_pix_count) != 1)
		{
			ERROR("Invalid value for --min-pix-count\n");
			return EINVAL;
		}
		break;

		case 310 :
		if (sscanf(arg, "%d", &args->iargs.max_pix_count) != 1)
		{
			ERROR("Invalid value for --max-pix-count\n");
			return EINVAL;
		}
		break;

		case 311 :
		if (sscanf(arg, "%d", &args->iargs.local_bg_radius) != 1)
		{
			ERROR("Invalid value for --local-bg-radius\n");
			return EINVAL;
		}
		break;

		case 312 :
		if (sscanf(arg, "%d", &args->iargs.min_res) != 1)
		{
			ERROR("Invalid value for --min-res\n");
			return EINVAL;
		}
		break;

		case 313 :
		if (sscanf(arg, "%d", &args->iargs.max_res) != 1)
		{
			ERROR("Invalid value for --max-res\n");
			return EINVAL;
		}
		break;

		case 314 :
		if (sscanf(arg, "%f", &args->iargs.min_snr_biggest_pix) != 1)
		{
			ERROR("Invalid value for --max-snr-biggest-pix\n");
			return EINVAL;
		}
		break;

		case 315 :
		if (sscanf(arg, "%f", &args->iargs.min_snr_peak_pix) != 1)
		{
			ERROR("Invalid value for --max-snr-peak-pix\n");
			return EINVAL;
		}
		break;

		case 316 :
		if (sscanf(arg, "%f", &args->iargs.min_sig) != 1)
		{
			ERROR("Invalid value for --max-ssig\n");
			return EINVAL;
		}
		break;

		case 317 :
		if (sscanf(arg, "%f", &args->iargs.min_peak_over_neighbour) != 1)
		{
			ERROR("Invalid value for --max-peak-over-neighbour\n");
			return EINVAL;
		}
		break;

		case 318 :
		args->iargs.use_saturated = 0;
		break;

		case 319 :
		args->iargs.no_revalidate = 1;
		break;

		case 320 :
		args->iargs.half_pixel_shift = 0;
		break;

		case 321 :
		args->iargs.check_hdf5_snr = 1;
		break;

		/* ---------- Indexing ---------- */

		case 400 :
		case 'z' :
		args->indm_str = strdup(arg);
		break;

		case 'p' :
		args->cellfile = strdup(arg);
		break;

		case 401 :
		/* Values in 'tols' are in frac (not %) and rad
		 * Conversion happens a few lines below */
		r = sscanf(arg, "%f,%f,%f,%f,%f,%f",
		           &args->iargs.tols[0], &args->iargs.tols[1], &args->iargs.tols[2],
		           &args->iargs.tols[3], &args->iargs.tols[4], &args->iargs.tols[5]);
		if ( r != 6 ) {
			/* Try old format */
			r = sscanf(arg, "%f,%f,%f,%f",
			           &args->iargs.tols[0], &args->iargs.tols[1],
			           &args->iargs.tols[2], &args->iargs.tols[3]);
			if ( r != 4 ) {
				ERROR("Invalid parameters for '--tolerance'\n");
				return EINVAL;
			}
			args->iargs.tols[4] = args->iargs.tols[3];
			args->iargs.tols[5] = args->iargs.tols[3];
		}

		/* Percent to fraction */
		args->iargs.tols[0] /= 100.0;
		args->iargs.tols[1] /= 100.0;
		args->iargs.tols[2] /= 100.0;
		args->iargs.tols[3] = deg2rad(args->iargs.tols[3]);
		args->iargs.tols[4] = deg2rad(args->iargs.tols[4]);
		args->iargs.tols[5] = deg2rad(args->iargs.tols[5]);
		break;

		case 402 :
		args->if_checkcell = 0;
		break;

		case 403 :
		args->if_checkcell = 1;  /* This is the default */
		break;

		case 404 :
		args->if_multi = 1;
		break;

		case 405 :
		args->if_multi = 0;  /* This is the default */
		break;

		case 406 :
		args->if_retry = 0;
		break;

		case 407 :
		args->if_retry = 1;  /* This is the default */
		break;

		case 408 :
		args->if_refine = 0;
		break;

		case 409 :
		args->if_refine = 1;  /* This is the default */
		break;

		case 410 :
		args->if_peaks = 0;
		break;

		case 411 :
		args->if_peaks = 1;  /* This is the default */
		break;

		case 412 :
		ERROR("The option --no-cell-combinations is no longer used.\n");
		/* .. but we can still carry on.  Results will probably be
		 *  better than the user expected. */
		break;

		case 413 :
		if (sscanf(arg, "%f", &args->iargs.wavelength_estimate) != 1)
		{
			ERROR("Invalid value for --wavelength-estimate\n");
			return EINVAL;
		}
		break;

		case 414 :
		if (sscanf(arg, "%d", &args->iargs.n_threads) != 1)
		{
			ERROR("Invalid value for --max-indexer-threads\n");
			return EINVAL;
		}
		break;

		case 415 :
		if (sscanf(arg, "%f", &args->iargs.clen_estimate) != 1)
		{
			ERROR("Invalid value for --camera-length-estimate\n");
			return EINVAL;
		}
		break;

		/* ---------- Integration ---------- */

		case 501 :
		args->iargs.int_meth = integration_method(arg, &r);
		if ( r ) {
			ERROR("Invalid integration method '%s'\n", arg);
			return EINVAL;
		}
		break;

		case 502 :
		if ( sscanf(arg, "%f", &args->iargs.fix_profile_r) != 1 ) {
			ERROR("Invalid value for --fix-profile-radius\n");
			return EINVAL;
		}
		break;

		case 503 :
		ERROR("The option --fix-bandwidth is no longer used.\n");
		ERROR("Set the bandwidth in the geometry file instead.\n");
		break;

		case 504 :
		if ( sscanf(arg, "%f", &args->iargs.fix_divergence) != 1 ) {
			ERROR("Invalid value for --fix-divergence\n");
			return EINVAL;
		}
		break;

		case 505 :
		r = sscanf(arg, "%f,%f,%f", &args->iargs.ir_inn,
		           &args->iargs.ir_mid, &args->iargs.ir_out);
		if ( r != 3 ) {
			ERROR("Invalid parameters for '--int-radius'\n");
			return EINVAL;
		}
		break;

		case 506 :
		if ( strcmp(arg, "random") == 0 ) {
			args->iargs.int_diag = INTDIAG_RANDOM;
		}

		if ( strcmp(arg, "all") == 0 ) {
			args->iargs.int_diag = INTDIAG_ALL;
		}

		if ( strcmp(arg, "negative") == 0 ) {
			args->iargs.int_diag = INTDIAG_NEGATIVE;
		}

		if ( strcmp(arg, "implausible") == 0 ) {
			args->iargs.int_diag = INTDIAG_IMPLAUSIBLE;
		}

		if ( strcmp(arg, "strong") == 0 ) {
			args->iargs.int_diag = INTDIAG_STRONG;
		}

		r = sscanf(arg, "%i,%i,%i", &args->iargs.int_diag_h,
		           &args->iargs.int_diag_k, &args->iargs.int_diag_l);
		if ( r == 3 ) {
			args->iargs.int_diag = INTDIAG_INDICES;
		}

		if ( (args->iargs.int_diag == INTDIAG_NONE)
		  && (strcmp(arg, "none") != 0) )
		{
			ERROR("Invalid value for --int-diag.\n");
			return EINVAL;
		}

		break;

		case 507 :
		if ( sscanf(arg, "%f", &args->iargs.push_res) != 1 ) {
			ERROR("Invalid value for --push-res\n");
			return EINVAL;
		}
		args->iargs.push_res *= 1e9;  /* nm^-1 -> m^-1 */
		break;

		case 508 :
		args->iargs.overpredict = 1;
		break;

		/* ---------- Output ---------- */

		case 601 :
		args->iargs.stream_nonhits = 0;
		break;

		case 602 :
		add_copy_header(args, arg);
		break;

		case 603 :
		args->iargs.stream_flags = CLEAR_BIT(args->iargs.stream_flags,
		                                     STREAM_PEAKS);
		break;

		case 604 :
		args->iargs.stream_flags = CLEAR_BIT(args->iargs.stream_flags,
		                                     STREAM_REFLECTIONS);
		break;

		case 605 :
		if ( sscanf(arg, "%d", &args->serial_start) != 1 ) {
			ERROR("Invalid value for --serial-start\n");
			return EINVAL;
		}
		break;

		case 606 :
		args->harvest_file = strdup(arg);
		break;

		default :
		return ARGP_ERR_UNKNOWN;

	}

	return 0;
}


int main(int argc, char *argv[])
{
	FILE *fh = NULL;
	Stream *st;
	struct indexamajig_arguments args;
	char *tmpdir;  /* e.g. /tmp/indexamajig.12345 */
	char *rn;  /* e.g. /home/taw/indexing */
	int r;
	int timeout = 240;
	TakeTwoOptions *taketwo_opts = NULL;
	FelixOptions *felix_opts = NULL;
	XGandalfOptions *xgandalf_opts = NULL;
	PinkIndexerOptions *pinkindexer_opts = NULL;
	FromFileOptions *fromfile_opts = NULL;
	double wl_from_dt;

	/* Defaults for "top level" arguments */
	args.filename = NULL;
	args.geom_filename = NULL;
	args.outfile = NULL;
	args.temp_location = strdup(".");
	args.prefix = strdup("");
	args.check_prefix = 1;
	args.n_proc = 1;
	args.cellfile = NULL;
	args.indm_str = NULL;
	args.basename = 0;
	args.zmq_addr = NULL;
	args.serial_start = 1;
	args.if_peaks = 1;
	args.if_multi = 0;
	args.if_retry = 1;
	args.if_refine = 1;
	args.if_checkcell = 1;
	args.profile = 0;
	args.copy_headers = NULL;
	args.n_copy_headers = 0;
	args.harvest_file = NULL;
	args.taketwo_opts_ptr = &taketwo_opts;
	args.felix_opts_ptr = &felix_opts;
	args.xgandalf_opts_ptr = &xgandalf_opts;
	args.pinkindexer_opts_ptr = &pinkindexer_opts;
	args.fromfile_opts_ptr = &fromfile_opts;

	/* Defaults for process_image arguments */
	args.iargs.cell = NULL;
	args.iargs.noisefilter = 0;
	args.iargs.median_filter = 0;
	args.iargs.tols[0] = 0.05;  /* frac (not %) */
	args.iargs.tols[1] = 0.05;  /* frac (not %) */
	args.iargs.tols[2] = 0.05;  /* frac (not %) */
	args.iargs.tols[3] = deg2rad(1.5); /* radians */
	args.iargs.tols[4] = deg2rad(1.5); /* radians */
	args.iargs.tols[5] = deg2rad(1.5); /* radians */
	args.iargs.threshold = 800.0;
	args.iargs.min_sq_gradient = 100000.0;
	args.iargs.min_snr = 5.0;
	args.iargs.min_pix_count = 2;
	args.iargs.max_pix_count = 200;
	args.iargs.min_res = 0;
	args.iargs.max_res = 1200;
	args.iargs.local_bg_radius = 3;
	args.iargs.min_snr_biggest_pix = 7.0;    /* peak finder 9  */
	args.iargs.min_snr_peak_pix = 6.0;
	args.iargs.min_sig = 11.0;
	args.iargs.min_peak_over_neighbour = -INFINITY;
	args.iargs.check_hdf5_snr = 0;
	args.iargs.dtempl = NULL;
	args.iargs.peaks = PEAK_ZAEF;
	args.iargs.half_pixel_shift = 1;
	args.iargs.pk_inn = -1.0;
	args.iargs.pk_mid = -1.0;
	args.iargs.pk_out = -1.0;
	args.iargs.ir_inn = -1.0;
	args.iargs.ir_mid = -1.0;
	args.iargs.ir_out = -1.0;
	args.iargs.use_saturated = 1;
	args.iargs.no_revalidate = 0;
	args.iargs.stream_flags = STREAM_PEAKS | STREAM_REFLECTIONS;
	args.iargs.stream_nonhits = 1;
	args.iargs.int_diag = INTDIAG_NONE;
	args.iargs.min_peaks = 0;
	args.iargs.overpredict = 0;
	args.iargs.wait_for_file = 0;
	args.iargs.ipriv = NULL;  /* No default */
	args.iargs.int_meth = integration_method("rings-nocen-nosat-nograd", NULL);
	args.iargs.push_res = +INFINITY;
	args.iargs.highres = +INFINITY;
	args.iargs.fix_profile_r = -1.0;
	args.iargs.fix_divergence = -1.0;
	args.iargs.no_image_data = 0;
	args.iargs.no_mask_data = 0;
	args.iargs.wavelength_estimate = NAN;
	args.iargs.clen_estimate = NAN;
	args.iargs.n_threads = 1;

	argp_program_version_hook = show_version;

	static char doc[] = "Index and integrate snapshot diffraction images.\v"
	                    "For more information including a tutorial, visit "
	                    "https://www.desy.de/~twhite/crystfel";

	static struct argp_option options[] = {

		{NULL, 0, 0, OPTION_DOC, "Basic options:", 2},

		{NULL, 'h', NULL, OPTION_HIDDEN, NULL},
		{NULL, 'v', NULL, OPTION_HIDDEN, NULL},

		{"input", 'i', "infile", 0, "List of input image filenames"},
		{"output", 'o', "filename.stream", 0, "Output stream filename"},
		{"geometry",'g', "experiment.geom", 0, "Detector geometry filename"},
		{"prefix", 'x', "/path/to/images/", OPTION_NO_USAGE, "Prefix filenames from input "
		        "file"},
		{NULL, 'j', "nproc", 0, "Run this many analyses in parallel, default 1"},
		{"basename", 201, NULL, OPTION_NO_USAGE, "Remove director parts from the "
		        "filenames"},
		{"no-check-prefix", 202, NULL, OPTION_NO_USAGE, "Don't attempt to correct the "
		        "--prefix"},
		{"highres", 203, "res", OPTION_NO_USAGE, "Absolute resolution cutoff in Angstroms"},
		{"profile", 204, NULL, OPTION_NO_USAGE, "Show timing data for performance "
		        "monitoring"},
		{"temp-dir", 205, "path", OPTION_NO_USAGE, "Location for temporary folder"},
		{"wait-for-file", 206, "seconds", OPTION_NO_USAGE, "Wait for each file before "
		        "processing"},
		{"zmq-input", 207, "addr", OPTION_NO_USAGE, "Receive data over ZeroMQ from "
			"this location"},
		{"no-image-data", 208, NULL, OPTION_NO_USAGE, "Do not load image data"},
		{"spectrum-file", 209, "fn", OPTION_NO_USAGE | OPTION_HIDDEN,
		       "File containing radiation spectrum"},
		{"no-mask-data", 210, NULL, OPTION_NO_USAGE, "Do not load mask data"},

		{NULL, 0, 0, OPTION_DOC, "Peak search options:", 3},
		{"peaks", 301, "method", 0, "Peak search method.  Default: zaef"},
		{"peak-radius", 302, "r1,r2,r3", OPTION_NO_USAGE, "Radii for peak search"},
		{"min-peaks", 303, "n", OPTION_NO_USAGE, "Minimum number of peaks for indexing"},
		{"hdf5-peaks", 304, "p", OPTION_HIDDEN, "Location of peak table in HDF5 file"},
		{"median-filter", 305, "n", OPTION_NO_USAGE, "Apply median filter to image data"},
		{"filter-noise", 306, NULL, OPTION_NO_USAGE, "Apply noise filter to image data"},
		{"threshold", 't', "adu", OPTION_NO_USAGE, "Threshold for peak detection "
		        "(zaef only, default 800)"},
		{"min-squared-gradient", 307, "n", OPTION_NO_USAGE, "Minimum squared gradient "
		        "(zaef only, default 100000)"},
		{"min-gradient", 307, "n", OPTION_ALIAS | OPTION_HIDDEN, NULL},
		{"min-snr", 308, "n", OPTION_NO_USAGE, "Minimum signal/noise ratio for peaks "
		        "(zaef,peakfinder8,peakfinder9 only, default 5)"},
		{"min-pix-count", 309, "n", OPTION_NO_USAGE, "Minimum number of pixels per peak "
		        "(peakfinder8 only, default 2)"},
		{"max-pix-count", 310, "n", OPTION_NO_USAGE, "Maximum number of pixels per peak "
		        "(peakfinder8 only, default 2)"},
		{"local-bg-radius", 311, "n", OPTION_NO_USAGE, "Radius (pixels) for local "
		        "background estimation (peakfinder8/9 only, default 3)"},
		{"min-res", 312, "n", OPTION_NO_USAGE, "Minimum resoultion (pixels) for peak "
		        "search (peakfinder8 only, default 0)"},
		{"max-res", 313, "n", OPTION_NO_USAGE, "Maximum resoultion (pixels) for peak "
		        "search (peakfinder8 only, default 1200)"},
		{"min-snr-biggest-pix", 314, "n", OPTION_NO_USAGE, "Minimum SNR of the biggest "
		        "pixel in the peak (peakfinder9 only)"},
		{"min-snr-peak-pix", 315, "n", OPTION_NO_USAGE, "Minimum SNR of peak pixel "
		        "(peakfinder9 only)"},
		{"min-sig", 316, "n", OPTION_NO_USAGE, "Minimum standard deviation of the "
		        "background (peakfinder9 only)"},
		{"min-peak-over-neighbour", 317, "n", OPTION_NO_USAGE, "Minimum difference between "
		        "highest pixel and neighbours (peakfinder9 only, just for speed)"},
		{"no-use-saturated", 318, NULL, OPTION_NO_USAGE, "Reject saturated peaks"},
		{"no-revalidate", 319, NULL, OPTION_NO_USAGE, "Don't re-integrate and check HDF5 "
		        "peaks"},
		{"no-half-pixel-shift", 320, NULL, OPTION_NO_USAGE, "Don't offset HDF5 peak "
		        "locations by 0.5 pixels"},
		{"check-hdf5-snr", 321, NULL, OPTION_NO_USAGE, "Check SNR for peaks from HDF5 or "
		        "CXI (see --min-snr)"},

		{NULL, 0, 0, OPTION_DOC, "Indexing options:", 4},
		{"indexing", 400, "method", 0, "List of indexing methods"},
		{NULL, 'z', "method", OPTION_HIDDEN | OPTION_ALIAS, NULL},
		{"pdb", 'p', "parameters.cell", 0, "PDB or CrystFEL Unit Cell File"},
		{"tolerance", 401, "a,b,c,al,be,ga", OPTION_NO_USAGE, "Tolerances for cell "
		        "comparison in percent and degrees, default 5,5,5,1.5,1.5,1.5"},
		{"no-check-cell", 402, NULL, OPTION_NO_USAGE, "Don't check cell parameters "
		        "against target cell"},
		{"check-cell", 403, NULL, OPTION_HIDDEN, NULL},
		{"multi", 404, NULL, OPTION_NO_USAGE, "Repeat indexing to index multiple hits"},
		{"no-multi", 405, NULL, OPTION_HIDDEN, NULL},
		{"no-retry", 406, NULL, OPTION_NO_USAGE, "Don't repeat indexing to increase "
		        "indexing rate"},
		{"retry", 407, NULL, OPTION_HIDDEN, NULL},
		{"no-refine", 408, NULL, OPTION_NO_USAGE, "Skip prediction refinement"},
		{"refine", 409, NULL, OPTION_HIDDEN, NULL},
		{"no-check-peaks", 410, NULL, OPTION_NO_USAGE, "Don't check that most peaks can be "
		        "accounted for by the indexing solution"},
		{"check-peaks", 411, NULL, OPTION_HIDDEN, NULL},
		{"no-cell-combinations", 412, NULL, OPTION_HIDDEN, NULL},
		{"wavelength-estimate", 413, "metres", 0,
		        "Estimate of the incident radiation wavelength, in metres."},
		{"max-indexer-threads", 414, "n", 0,
		        "Maximum number of threads allowed for indexing engines."},
		{"camera-length-estimate", 415, "metres", 0,
		        "Estimate of the camera length, in metres."},

		{NULL, 0, 0, OPTION_DOC, "Integration options:", 5},
		{"integration", 501, "method", OPTION_NO_USAGE, "Integration method"},
		{"fix-profile-radius", 502, "r", OPTION_NO_USAGE, "Fix profile radius for spot "
		        "prediction, instead of automatically determining"},
		{"fix-bandwidth", 503, "bw", OPTION_NO_USAGE, "Set the bandwidth for spot "
		        "prediction"},
		{"fix-divergence", 504, "deg", OPTION_NO_USAGE, "Set the divergence (full angle) "
		        "for spot prediction"},
		{"int-radius", 505, "r1,r2,r3", 0, "Set the integration radii (inner,mid,outer)"},
		{"int-diag", 506, "condition", 0, "Show debugging information about reflections"},
		{"push-res", 507, "dist", 0, "Integrate higher than apparent resolution cutoff (m^-1)"},
		{"overpredict", 508, NULL, 0, "Over-predict reflections"},

		{NULL, 0, 0, OPTION_DOC, "Output options:", 6},
		{"no-non-hits-in-stream", 601, NULL, OPTION_NO_USAGE, "Don't include non-hits in "
		        "stream (see --min-peaks)"},
		{"copy-hdf5-field", 602, "f", OPTION_HIDDEN, NULL},
		{"copy-header", 602, "f", OPTION_NO_USAGE, "Put the value of this image header "
		        "field into the stream"},
		{"no-peaks-in-stream", 603, NULL, OPTION_NO_USAGE, "Don't put peak search results "
		        "in stream"},
		{"no-refls-in-stream", 604, NULL, OPTION_NO_USAGE, "Don't put integration results "
		        "in stream"},
		{"serial-start", 605, "n", OPTION_NO_USAGE, "Start the serial numbers in the stream "
		        "here"},
		{"harvest-file", 606, "filename", OPTION_NO_USAGE, "Write the actual parameters "
			"used in JSON format"},

		{NULL, 0, 0, OPTION_DOC, "More information:", 99},

		{0}

	};

	static struct argp_child argp_children[] = {
		{&taketwo_argp, 0, NULL, -2},
		{&felix_argp, 0, NULL, -2},
		{&xgandalf_argp, 0, NULL, -2},
		{&pinkIndexer_argp, 0, NULL, -2},
		{&fromfile_argp, 0, NULL, -2},
		{0}
	};

	static struct argp argp = { options, parse_arg, NULL, doc,
	                            argp_children, NULL, NULL };
	if ( argp_parse(&argp, argc, argv, 0, NULL, &args) ) return 1;

	/* Check for minimal information */
	if ( (args.filename == NULL) && (args.zmq_addr == NULL) ) {
		ERROR("You need to provide the input filename (use -i)\n");
		return 1;
	}
	if ( args.geom_filename == NULL ) {
		ERROR("You need to specify the geometry filename (use -g)\n");
		return 1;
	}
	if ( args.outfile == NULL ) {
		ERROR("You need to specify the output filename (use -o)\n");
		return 1;
	}

	if ( (args.filename != NULL) && (args.zmq_addr != NULL) ) {
		ERROR("You must only specify one of --input and --zmq-input.\n");
		return 1;
	}

	/* Open input */
	if ( args.filename != NULL ) {
		if ( strcmp(args.filename, "-") == 0 ) {
			fh = stdin;
		} else {
			fh = fopen(args.filename, "r");
		}
		if ( fh == NULL ) {
			ERROR("Failed to open input file '%s'\n", args.filename);
			return 1;
		}
	}

	/* Check prefix (if given) */
	if ( args.check_prefix ) {
		args.prefix = check_prefix(args.prefix);
	}

	/* Check number of processes */
	if ( args.n_proc == 0 ) {
		ERROR("Invalid number of processes.\n");
		return 1;
	}

	/* Load data template (new API) */
	args.iargs.dtempl = data_template_new_from_file(args.geom_filename);
	if ( args.iargs.dtempl == NULL ) {
		ERROR("Failed to read detector geometry from '%s'"
		      " (for new API)\n", args.geom_filename);
		return 1;
	}

	/* Add any headers we need to copy */
	for ( r=0; r<args.n_copy_headers; r++ ) {
		data_template_add_copy_header(args.iargs.dtempl,
		                              args.copy_headers[r]);
	}

	/* If no integration radii were given, apply the defaults */
	if ( args.iargs.ir_inn < 0 ) {
		STATUS("WARNING: You did not specify --int-radius.\n");
		STATUS("WARNING: I will use the default values, which are"
		       " probably not appropriate for your patterns.\n");
		args.iargs.ir_inn = 4.0;
		args.iargs.ir_mid = 5.0;
		args.iargs.ir_out = 7.0;
	}

	/* If no peak radii were given, copy the integration radii */
	if ( args.iargs.pk_inn < 0.0 ) {
		args.iargs.pk_inn = args.iargs.ir_inn;
		args.iargs.pk_mid = args.iargs.ir_mid;
		args.iargs.pk_out = args.iargs.ir_out;
	}

	/* Load unit cell (if given) */
	if ( args.cellfile != NULL ) {
		args.iargs.cell = load_cell_from_file(args.cellfile);
		if ( args.iargs.cell == NULL ) {
			ERROR("Couldn't read unit cell (from %s)\n", args.cellfile);
			return 1;
		}
		free(args.cellfile);
	} else {
		args.iargs.cell = NULL;
	}

	tmpdir = create_tempdir(args.temp_location);
	if ( tmpdir == NULL ) return 1;

	/* Change into temporary folder, temporarily, to control the crap
	 * dropped by indexing programs during setup */
	rn = getcwd(NULL, 0);
	r = chdir(tmpdir);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		return 1;
	}

	/* Auto-detect indexing methods if 'requested' */
	if ( args.indm_str == NULL ) {

		STATUS("No indexing methods specified.  I will try to ");
		STATUS("automatically detect the available methods.\n");
		STATUS("To disable auto-detection of indexing methods, specify ");
		STATUS("which methods to use with --indexing=<methods>.\n");
		STATUS("Use --indexing=none to disable indexing and integration.\n");

		args.indm_str = detect_indexing_methods(args.iargs.cell);

	}

	wl_from_dt = data_template_get_wavelength_if_possible(args.iargs.dtempl);
	if ( !isnan(wl_from_dt) ) {
		if ( !isnan(args.iargs.wavelength_estimate) ) {
			ERROR("WARNING: Ignoring your value for "
			      "--wavelength-estimate because the geometry file "
			      "already contains a static value.\n");
		}
		args.iargs.wavelength_estimate = wl_from_dt;
	}

	/* Prepare the indexing system */
	if ( args.indm_str == NULL ) {

		ERROR("No indexing method specified, and no usable indexing ");
		ERROR("methods auto-detected.\n");
		ERROR("Install some indexing programs (mosflm,dirax etc), or ");
		ERROR("try again with --indexing=none.\n");
		return 1;

	} else if ( strcmp(args.indm_str, "none") == 0 ) {

		STATUS("Indexing/integration disabled.\n");
		if ( args.iargs.cell != NULL ) {
			STATUS("Ignoring your unit cell.\n");
		}
		args.iargs.ipriv = NULL;

	} else {

		int i, n;
		const IndexingMethod *methods;
		IndexingFlags flags = 0;

		if ( args.iargs.cell != NULL ) {
			STATUS("This is what I understood your unit cell to be:\n");
			cell_print(args.iargs.cell);
		} else {
			STATUS("No reference unit cell provided.\n");
		}

		if ( args.if_checkcell ) {
			flags |= INDEXING_CHECK_CELL;
		}
		if ( args.if_refine ) {
			flags |= INDEXING_REFINE;
		}
		if ( args.if_peaks ) {
			flags |= INDEXING_CHECK_PEAKS;
		}
		if ( args.if_multi ) {
			flags |= INDEXING_MULTI;
		}
		if ( args.if_retry ) {
			flags |= INDEXING_RETRY;
		}

		args.iargs.ipriv = setup_indexing(args.indm_str,
		                                  args.iargs.cell,
		                                  args.iargs.tols,
		                                  flags,
		                                  args.iargs.wavelength_estimate,
		                                  args.iargs.clen_estimate,
		                                  args.iargs.n_threads,
		                                  taketwo_opts,
		                                  xgandalf_opts,
		                                  pinkindexer_opts,
		                                  felix_opts,
		                                  fromfile_opts);

		free(args.filename);

		if ( args.iargs.ipriv == NULL ) {
			ERROR("Failed to set up indexing system\n");
			return 1;
		}

		methods = indexing_methods(args.iargs.ipriv, &n);
		for ( i=0; i<n; i++ ) {
			if ( (methods[i] & INDEXING_METHOD_MASK) == INDEXING_PINKINDEXER ) {
				/* Extend timeout if using pinkIndexer */
				timeout = 3000;
				break;
			}
		}

	}

	/* Change back to where we were before.  Sandbox code will create
	 * worker subdirs inside the temporary folder, and process_image will
	 * change into them. */
	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		return 1;
	}
	free(rn);

	/* Open output stream */
	st = stream_open_for_write(args.outfile, args.iargs.dtempl);
	if ( st == NULL ) {
		ERROR("Failed to open stream '%s'\n", args.outfile);
		return 1;
	}

	/* Write audit info */
	stream_write_commandline_args(st, argc, argv);
	stream_write_geometry_file(st, args.geom_filename);
	stream_write_target_cell(st, args.iargs.cell);
	stream_write_indexing_methods(st, args.indm_str);

	if ( (args.harvest_file != NULL) && (args.serial_start <= 1) ) {
		write_harvest_file(&args.iargs, args.harvest_file,
		                   args.if_multi, args.if_refine, args.if_retry,
		                   args.if_peaks, args.if_checkcell);
	}

	free(args.outfile);
	free(args.indm_str);

	gsl_set_error_handler_off();

	r = create_sandbox(&args.iargs, args.n_proc, args.prefix, args.basename,
	                   fh, st, tmpdir, args.serial_start, args.zmq_addr,
	                   timeout, args.profile);

	cell_free(args.iargs.cell);
	free(args.prefix);
	free(args.temp_location);
	free(tmpdir);
	data_template_free(args.iargs.dtempl);
	stream_close(st);
	cleanup_indexing(args.iargs.ipriv);

	return r;
}
