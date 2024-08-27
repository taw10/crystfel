/*
 * im-argparse.c
 *
 * Command line argument parsing for indexamajig
 *
 * Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2023 Thomas White <taw@physics.org>
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

#include "version.h"


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
	struct im_zmq_params zmq_params;
	struct im_asapo_params asapo_params;
	int serial_start;
	char *temp_location;
	int if_refine;
	int if_checkcell;
	int if_peaks;
	int if_multi;
	int if_retry;
	int profile;  /* Whether to do wall-clock time profiling */
	int no_data_timeout;
	char **copy_headers;
	int n_copy_headers;
	char *harvest_file;
	int cpu_pin;

	struct taketwo_options **taketwo_opts_ptr;
	struct felix_options **felix_opts_ptr;
	struct xgandalf_options **xgandalf_opts_ptr;
	struct pinkindexer_options **pinkindexer_opts_ptr;
	struct fromfile_options **fromfile_opts_ptr;
	struct asdf_options **asdf_opts_ptr;
};


static void show_version(FILE *fh, struct argp_state *state)
{
	printf("CrystFEL: %s\n", crystfel_version_string());
	printf("%s\n", crystfel_licence_string());
}


static DataSourceType parse_data_format(const char *str)
{
	if ( strcmp(str, "hdf5") == 0 ) return DATA_SOURCE_TYPE_HDF5;
	if ( strcmp(str, "msgpack") == 0 ) return DATA_SOURCE_TYPE_MSGPACK;
	if ( strcmp(str, "seedee") == 0 ) return DATA_SOURCE_TYPE_SEEDEE;
	/* CBF and CBFGZ should be added here once image-cbf.c supports
	 * in-memory access */
	return DATA_SOURCE_TYPE_UNKNOWN;
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
		state->child_inputs[5] = args->asdf_opts_ptr;
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
		free(args->prefix);
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
		args->zmq_params.addr = strdup(arg);
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

		case 211 :
		if ( args->zmq_params.n_subscriptions == 256 ) {
			ERROR("Too many ZMQ subscriptions.\n");
			return 1;
		}
		args->zmq_params.subscriptions[args->zmq_params.n_subscriptions++] = strdup(arg);
		break;

		case 212 :
		args->zmq_params.request = strdup(arg);
		break;

		case 213 :
		args->asapo_params.endpoint = strdup(arg);
		break;

		case 214 :
		args->asapo_params.token = strdup(arg);
		break;

		case 215 :
		args->asapo_params.beamtime = strdup(arg);
		break;

		case 217 :
		args->asapo_params.group_id = strdup(arg);
		break;

		case 218 :
		args->asapo_params.source = strdup(arg);
		break;

		case 219 :
		args->iargs.data_format = parse_data_format(arg);
		if ( args->iargs.data_format == DATA_SOURCE_TYPE_UNKNOWN ) {
			ERROR("Unrecognised data format '%s'\n", arg);
			return EINVAL;
		}
		break;

		case 220 :
		args->asapo_params.stream = strdup(arg);
		break;

		case 221 :
		args->asapo_params.wait_for_stream = 1;
		break;

		case 222 :
		args->asapo_params.write_output_stream = 1;
		break;

		case 223 :
		args->cpu_pin = 1;
		break;

		case 224 :
		if (sscanf(arg, "%d", &args->no_data_timeout) != 1)
		{
			ERROR("Invalid value for --no-data-timeout\n");
			return EINVAL;
		}
		break;

		case 225 :
		if (sscanf(arg, "%d", &args->asapo_params.consumer_timeout_ms) != 1)
		{
			ERROR("Invalid value for --asapo-consumer-timeout\n");
			return EINVAL;
		}
		break;

		case 226 :
		args->asapo_params.use_ack = 1;
		break;

		/* ---------- Peak search ---------- */

		case 't' :
		args->iargs.peak_search.threshold = strtof(arg, NULL);
		break;

		case 301 :
		args->iargs.peak_search.method = parse_peaksearch(arg);
		if ( args->iargs.peak_search.method == PEAK_ERROR ) {
			ERROR("Unrecognised peak detection method '%s'\n", arg);
			return EINVAL;
		}
		break;

		case 302 :
		r = sscanf(arg, "%f,%f,%f", &args->iargs.peak_search.pk_inn,
		           &args->iargs.peak_search.pk_mid, &args->iargs.peak_search.pk_out);
		if ( (r != 3) || (args->iargs.peak_search.pk_inn < 0) ) {
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
		ERROR("The option --hdf5-peaks is no longer used.\n");
		ERROR("Set the peak path in the geometry file.\n");
		break;

		case 305 :
		if (sscanf(arg, "%d", &args->iargs.peak_search.median_filter) != 1)
		{
			ERROR("Invalid value for --median-filter\n");
			return EINVAL;
		}
		break;

		case 306 :
		args->iargs.peak_search.noisefilter = 1;
		break;

		case 307 :
		if (sscanf(arg, "%f", &args->iargs.peak_search.min_sq_gradient) != 1)
		{
			ERROR("Invalid value for --min-squared-gradient\n");
			return EINVAL;
		}
		break;

		case 308 :
		if (sscanf(arg, "%f", &args->iargs.peak_search.min_snr) != 1)
		{
			ERROR("Invalid value for --min-snr\n");
			return EINVAL;
		}
		break;

		case 309 :
		if (sscanf(arg, "%d", &args->iargs.peak_search.min_pix_count) != 1)
		{
			ERROR("Invalid value for --min-pix-count\n");
			return EINVAL;
		}
		break;

		case 310 :
		if (sscanf(arg, "%d", &args->iargs.peak_search.max_pix_count) != 1)
		{
			ERROR("Invalid value for --max-pix-count\n");
			return EINVAL;
		}
		break;

		case 311 :
		if (sscanf(arg, "%d", &args->iargs.peak_search.local_bg_radius) != 1)
		{
			ERROR("Invalid value for --local-bg-radius\n");
			return EINVAL;
		}
		break;

		case 312 :
		if (sscanf(arg, "%d", &args->iargs.peak_search.min_res) != 1)
		{
			ERROR("Invalid value for --min-res\n");
			return EINVAL;
		}
		break;

		case 313 :
		if (sscanf(arg, "%d", &args->iargs.peak_search.max_res) != 1)
		{
			ERROR("Invalid value for --max-res\n");
			return EINVAL;
		}
		break;

		case 314 :
		if (sscanf(arg, "%f", &args->iargs.peak_search.min_snr_biggest_pix) != 1)
		{
			ERROR("Invalid value for --max-snr-biggest-pix\n");
			return EINVAL;
		}
		break;

		case 315 :
		if (sscanf(arg, "%f", &args->iargs.peak_search.min_snr_peak_pix) != 1)
		{
			ERROR("Invalid value for --max-snr-peak-pix\n");
			return EINVAL;
		}
		break;

		case 316 :
		if (sscanf(arg, "%f", &args->iargs.peak_search.min_sig) != 1)
		{
			ERROR("Invalid value for --max-ssig\n");
			return EINVAL;
		}
		break;

		case 317 :
		if (sscanf(arg, "%f", &args->iargs.peak_search.min_peak_over_neighbour) != 1)
		{
			ERROR("Invalid value for --max-peak-over-neighbour\n");
			return EINVAL;
		}
		break;

		case 318 :
		args->iargs.peak_search.use_saturated = 0;
		break;

		case 319 :
		args->iargs.peak_search.revalidate = 0;
		break;

		case 320 :
		args->iargs.peak_search.half_pixel_shift = 0;
		break;

		case 321 :
		args->iargs.peak_search.check_hdf5_snr = 1;
		break;

		case 322:
		args->iargs.peak_search.peakfinder8_fast = 1;
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

		case 416 :
		args->iargs.mille = 1;
		break;

		case 417 :
		free(args->milledir);
		args->milledir = strdup(arg);
		break;

		case 418 :
		if (sscanf(arg, "%d", &args->iargs.max_mille_level) != 1)
		{
			ERROR("Invalid value for --max-mille-level\n");
			return EINVAL;
		}
		if ( args->iargs.max_mille_level < 0 ) {
			ERROR("Invalid value for --max-mille-level\n");
			return EINVAL;
		}
		break;

		case 419 :
		free(args->millefile);
		args->millefile = strdup(arg);
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

		case 509 :
		args->iargs.cell_params_only = 1;
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


struct indexamajig_arguments *parse_args(int argc, char *argv[])
{
	struct indexamajig_arguments args;
	int r;
	struct taketwo_options *taketwo_opts = NULL;
	struct felix_options *felix_opts = NULL;
	struct xgandalf_options *xgandalf_opts = NULL;
	struct pinkindexer_options *pinkindexer_opts = NULL;
	struct fromfile_options *fromfile_opts = NULL;
	struct asdf_options *asdf_opts = NULL;
	int err = 0;

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
	args.zmq_params.addr = NULL;
	args.zmq_params.request = NULL;
	args.zmq_params.n_subscriptions = 0;
	args.asapo_params.endpoint = NULL;
	args.asapo_params.token = NULL;
	args.asapo_params.beamtime = NULL;
	args.asapo_params.group_id = NULL;
	args.asapo_params.source = NULL;
	args.asapo_params.stream = NULL;
	args.asapo_params.wait_for_stream = 0;
	args.asapo_params.write_output_stream = 0;
	args.asapo_params.consumer_timeout_ms = 500;
	args.asapo_params.use_ack = 0;
	args.cpu_pin = 0;
	args.serial_start = 1;
	args.if_peaks = 1;
	args.if_multi = 0;
	args.if_retry = 1;
	args.if_refine = 1;
	args.if_checkcell = 1;
	args.profile = 0;
	args.no_data_timeout = 60;
	args.copy_headers = NULL;
	args.n_copy_headers = 0;
	args.harvest_file = NULL;
	args.taketwo_opts_ptr = &taketwo_opts;
	args.felix_opts_ptr = &felix_opts;
	args.xgandalf_opts_ptr = &xgandalf_opts;
	args.pinkindexer_opts_ptr = &pinkindexer_opts;
	args.fromfile_opts_ptr = &fromfile_opts;
	args.asdf_opts_ptr = &asdf_opts;
	args.worker = 0;
	args.fd_stream = 0;
	args.fd_mille = 0;

	/* Defaults for process_image arguments */
	args.iargs.cell = NULL;
	args.iargs.peak_search.noisefilter = 0;
	args.iargs.peak_search.median_filter = 0;
	args.iargs.tols[0] = 0.05;  /* frac (not %) */
	args.iargs.tols[1] = 0.05;  /* frac (not %) */
	args.iargs.tols[2] = 0.05;  /* frac (not %) */
	args.iargs.tols[3] = deg2rad(1.5); /* radians */
	args.iargs.tols[4] = deg2rad(1.5); /* radians */
	args.iargs.tols[5] = deg2rad(1.5); /* radians */
	args.iargs.peak_search.threshold = 800.0;
	args.iargs.peak_search.min_sq_gradient = 100000.0;
	args.iargs.peak_search.min_snr = 5.0;
	args.iargs.peak_search.min_pix_count = 2;
	args.iargs.peak_search.max_pix_count = 200;
	args.iargs.peak_search.min_res = 0;
	args.iargs.peak_search.max_res = 1200;
	args.iargs.peak_search.local_bg_radius = 3;
	args.iargs.peak_search.min_snr_biggest_pix = 7.0;    /* peak finder 9  */
	args.iargs.peak_search.min_snr_peak_pix = 6.0;
	args.iargs.peak_search.min_sig = 11.0;
	args.iargs.peak_search.min_peak_over_neighbour = -INFINITY;
	args.iargs.peak_search.check_hdf5_snr = 0;
	args.iargs.peak_search.peakfinder8_fast = 0;
	args.iargs.pf_private = NULL;
	args.iargs.dtempl = NULL;
	args.iargs.peak_search.method = PEAK_ZAEF;
	args.iargs.peak_search.half_pixel_shift = 1;
	args.iargs.peak_search.pk_inn = -1.0;
	args.iargs.peak_search.pk_mid = -1.0;
	args.iargs.peak_search.pk_out = -1.0;
	args.iargs.ir_inn = -1.0;
	args.iargs.ir_mid = -1.0;
	args.iargs.ir_out = -1.0;
	args.iargs.peak_search.use_saturated = 1;
	args.iargs.peak_search.revalidate = 1;
	args.iargs.stream_flags = STREAM_PEAKS | STREAM_REFLECTIONS;
	args.iargs.stream_nonhits = 1;
	args.iargs.int_diag = INTDIAG_NONE;
	args.iargs.min_peaks = 0;
	args.iargs.overpredict = 0;
	args.iargs.cell_params_only = 0;
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
	args.iargs.data_format = DATA_SOURCE_TYPE_UNKNOWN;
	args.iargs.mille = 0;
	args.iargs.milledir = strdup(".");
	args.iargs.max_mille_level = 99;

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
		{"zmq-subscribe", 211, "tag", OPTION_NO_USAGE, "Subscribe to ZMQ message"
			"type"},
		{"zmq-request", 212, "str", OPTION_NO_USAGE, "Request messages using"
			"this string."},
		{"asapo-endpoint", 213, "str", OPTION_NO_USAGE, "ASAP::O endpoint"},
		{"asapo-token", 214, "str", OPTION_NO_USAGE, "ASAP::O token"},
		{"asapo-beamtime", 215, "str", OPTION_NO_USAGE, "ASAP::O beamtime ID"},
		{"asapo-group", 217, "str", OPTION_NO_USAGE, "ASAP::O group ID"},
		{"asapo-source", 218, "str", OPTION_NO_USAGE, "ASAP::O data source"},
		{"data-format", 219, "str", OPTION_NO_USAGE, "Streamed data format"},
		{"asapo-stream", 220, "str", OPTION_NO_USAGE, "ASAP::O stream name"},
		{"asapo-wait-for-stream", 221, NULL, OPTION_NO_USAGE,
		        "Wait for ASAP::O stream to appear"},
		{"asapo-output-stream", 222, NULL, OPTION_NO_USAGE,
			"Create an ASAP::O hits-only stream"},
		{"cpu-pin", 223, NULL, OPTION_NO_USAGE, "Pin worker processes to CPUs"},
		{"no-data-timeout", 224, "s", OPTION_NO_USAGE,
			"Shut down after this many seconds without ASAP::O data"},
		{"asapo-consumer-timeout", 225, "ms", OPTION_NO_USAGE,
			"ASAP::O get_next timeout for one frame (milliseconds)"},
		{"asapo-acks", 226, NULL, OPTION_NO_USAGE, "Use ASAP::O acknowledgements"},

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
		        "or MsgPack peaks"},
		{"no-half-pixel-shift", 320, NULL, OPTION_NO_USAGE, "Don't offset HDF5 peak "
		        "locations by 0.5 pixels"},
		{"check-hdf5-snr", 321, NULL, OPTION_NO_USAGE, "Check SNR for peaks from HDF5, "
		        "CXI or MsgPack (see --min-snr)"},
		{"peakfinder8-fast", 322, NULL, OPTION_NO_USAGE, "peakfinder8 fast execution"},

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
		{"mille", 416, NULL, 0,
		        "Generate data for detector geometry refinement using Millepede"},
		{"mille-dir", 417, "dirname", 0, "Save Millepede data in folder"},
		{"max-mille-level", 418, "n", 0, "Maximum geometry refinement level"},

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
		{"cell-parameters-only", 509, NULL, 0, "Don't predict reflections at all"},

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

		{"worker", 701, NULL, OPTION_HIDDEN, "Be a worker process"},
		{"fd-input", 702, NULL, OPTION_HIDDEN, "File descriptor for stream"},
		{"fd-stream", 703, NULL, OPTION_HIDDEN, "File descriptor for input"},
		{"shm-queue", 704, NULL, OPTION_HIDDEN, "SHM handle for queue structure"},
		{"fd-mille", 705, NULL, OPTION_HIDDEN, "File descriptor for Mille data"},

		{NULL, 0, 0, OPTION_DOC, "More information:", 99},

		{0}

	};

	static struct argp_child argp_children[] = {
		{&taketwo_argp, 0, NULL, -2},
		{&felix_argp, 0, NULL, -2},
		{&xgandalf_argp, 0, NULL, -2},
		{&pinkIndexer_argp, 0, NULL, -2},
		{&fromfile_argp, 0, NULL, -2},
		{&asdf_argp, 0, NULL, -2},
		{0}
	};

	#if defined(__APPLE__)
	int q;
	for ( q=0; q<argc; q++ ) {
		if ( strcmp(argv[q], "--help") == 0 ) {
			fprintf(stderr, "\n"
			        "WARNING: 'indexamajig --help' crashes on some Mac OS versions.\n"
			        "This is a known problem, and is due to a bug in external library "
			        "code.\n\n");
			fflush(stderr);
			break;
		}
	}
	#endif

	static struct argp argp = { options, parse_arg, NULL, doc,
	                            argp_children, NULL, NULL };
	if ( argp_parse(&argp, argc, argv, 0, NULL, &args) ) return NULL;

	return args;
}
