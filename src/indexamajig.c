/*
 * indexamajig.c
 *
 * Index patterns, output hkl+intensity etc.
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
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
#include <getopt.h>
#include <hdf5.h>
#include <gsl/gsl_errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "version.h"
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
#include "cell-utils.h"
#include "integration.h"

#include "im-sandbox.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Process and index FEL diffraction images.\n"
"\n"
" -h, --help               Display this help message.\n"
"     --version            Print CrystFEL version number and exit.\n"
"\n"
" -i, --input=<filename>   Specify file containing list of images to process.\n"
"                           '-' means stdin, which is the default.\n"
" -o, --output=<filename>  Write output stream to this file. '-' for stdout.\n"
"                           Default: indexamajig.stream\n"
"\n"
"     --indexing=<methods> Use 'methods' for indexing.  Provide one or more\n"
"                           methods separated by commas.\n"
"                           See 'man indexamajig' for details.\n"
" -g. --geometry=<file>    Get detector geometry from file.\n"
" -p, --pdb=<file>         File (PDB or CrystFEL unit cell format) from which\n"
"                           to get the unit cell.  Default: 'molecule.pdb'.\n"
"     --basename           Remove the directory parts of the filenames.\n"
" -x, --prefix=<p>         Prefix filenames from input file with <p>.\n"
"     --peaks=<method>     Use 'method' for finding peaks.  Choose from:\n"
"                           zaef  : Use Zaefferer (2000) gradient detection.\n"
"                                    This is the default method.\n"
"                           hdf5  : Get from a table in HDF5 file.\n"
"                           cxi   : Get from CXI format HDF5 file.\n"
"     --hdf5-peaks=<p>     Find peaks table in HDF5 file here.\n"
"                           Default: /processing/hitfinder/peakinfo\n"
"     --integration=<meth> Perform final pattern integration using <meth>.\n"
"\n\n"
"For more control over the process, you might need:\n\n"
"    --tolerance=<tol>   Set the tolerances for cell comparison.\n"
"                          Default: 5,5,5,1.5.\n"
"    --filter-noise      Apply an aggressive noise filter which sets all\n"
"                         pixels in each 3x3 region to zero if any of them\n"
"                         have negative values.  Intensity measurement will\n"
"                         be performed on the image as it was before this.\n"
"    --median-filter=<n> Apply a median filter to the image data.  Intensity\n"
"                         measurement will be performed on the image as it\n"
"                         was before this.  The side length of the median\n"
"                         filter box will be 2<n>+1.  Default: 0 (no filter).\n"
"    --no-sat-corr       Don't correct values of saturated peaks using a\n"
"                         table included in the HDF5 file.\n"
"    --threshold=<n>     Only accept peaks above <n> ADU.  Default: 800.\n"
"    --min-gradient=<n>  Minimum squared gradient for Zaefferer peak search.\n"
"                         Default: 100,000.\n"
"    --min-snr=<n>       Minimum signal-to-noise ratio for peaks.\n"
"                         Default: 5.\n"
"    --check-hdf5-snr    Check SNR for peaks from --peaks=hdf5.\n"
"    --peak-radius=<r>   Integration radii for peak search.\n"
"    --int-radius=<r>    Set the integration radii.  Default: 4,5,7.\n"
"    --push-res=<n>      Integrate higher than apparent resolution cutoff.\n"
"    --highres=<n>       Absolute resolution cutoff in Angstroms.\n"
"    --fix-profile-radius Fix the reciprocal space profile radius for spot\n"
"                         prediction (default: automatically determine.\n"
"    --fix-bandwidth     Set the bandwidth for spot prediction.\n"
"    --fix-divergence    Set the divergence (full angle) for spot prediction.\n"
"\n"
"\nFor time-resolved stuff, you might want to use:\n\n"
"     --copy-hdf5-field <f>  Copy the value of field <f> into the stream. You\n"
"                             can use this option as many times as you need.\n"
"\n"
"\nOptions for greater performance:\n\n"
" -j <n>                   Run <n> analyses in parallel.  Default 1.\n"
" --temp-dir=<path>        Put the temporary folder under <path>.\n"
"\n"
"\nOptions you probably won't need:\n\n"
"     --no-check-prefix    Don't attempt to correct the --prefix.\n"
"     --no-use-saturated   During the initial peak search, reject\n"
"                           peaks which contain pixels above max_adu.\n"
"     --no-revalidate      Don't re-integrate and check HDF5 peaks for\n"
"                           validity.\n"
"     --no-peaks-in-stream Do not record peak search results in the stream.\n"
"     --no-refls-in-stream Do not record integrated reflections in the stream.\n"
"     --int-diag=<cond>    Show debugging information about reflections.\n"
"     --no-refine          Skip the prediction refinement step.\n"
);
}


static void add_geom_beam_stuff_to_copy_hdf5(struct copy_hdf5_field *copyme,
                                             struct detector *det,
                                             struct beam_params *beam)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p = &det->panels[i];

		if ( p->clen_from != NULL ) {
			add_copy_hdf5_field(copyme, p->clen_from);
		}
	}

	if ( beam->photon_energy_from != NULL ) {
		add_copy_hdf5_field(copyme, beam->photon_energy_from);
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *outfile = NULL;
	FILE *fh;
	Stream *st;
	int config_checkprefix = 1;
	int config_basename = 0;
	int integrate_saturated = 0;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	char *indm_str = NULL;
	char *cellfile = NULL;
	char *prefix = NULL;
	char *speaks = NULL;
	char *toler = NULL;
	int n_proc = 1;
	struct index_args iargs;
	char *intrad = NULL;
	char *pkrad = NULL;
	char *int_str = NULL;
	char *tempdir = NULL;
	char *int_diag = NULL;
	char *geom_filename = NULL;
	struct beam_params beam;
	int have_push_res = 0;

	/* Defaults */
	iargs.cell = NULL;
	iargs.noisefilter = 0;
	iargs.median_filter = 0;
	iargs.satcorr = 1;
	iargs.tols[0] = 5.0;
	iargs.tols[1] = 5.0;
	iargs.tols[2] = 5.0;
	iargs.tols[3] = 1.5;
	iargs.threshold = 800.0;
	iargs.min_gradient = 100000.0;
	iargs.min_snr = 5.0;
	iargs.check_hdf5_snr = 0;
	iargs.det = NULL;
	iargs.peaks = PEAK_ZAEF;
	iargs.beam = &beam;
	iargs.hdf5_peak_path = NULL;
	iargs.copyme = NULL;
	iargs.pk_inn = -1.0;
	iargs.pk_mid = -1.0;
	iargs.pk_out = -1.0;
	iargs.ir_inn = 4.0;
	iargs.ir_mid = 5.0;
	iargs.ir_out = 7.0;
	iargs.use_saturated = 1;
	iargs.no_revalidate = 0;
	iargs.stream_peaks = 1;
	iargs.stream_refls = 1;
	iargs.int_diag = INTDIAG_NONE;
	iargs.copyme = new_copy_hdf5_field_list();
	if ( iargs.copyme == NULL ) {
		ERROR("Couldn't allocate HDF5 field list.\n");
		return 1;
	}
	iargs.indm = NULL;  /* No default */
	iargs.ipriv = NULL;  /* No default */
	iargs.int_meth = integration_method("rings-nocen", NULL);
	iargs.push_res = 0.0;
	iargs.highres = +INFINITY;
	iargs.fix_profile_r = -1.0;
	iargs.fix_bandwidth = -1.0;
	iargs.fix_divergence = -1.0;
	iargs.predict_refine = 1;

	/* Long options */
	const struct option longopts[] = {

		/* Options with long and short versions */
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,               'v'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"indexing",           1, NULL,               'z'},
		{"geometry",           1, NULL,               'g'},
		{"pdb",                1, NULL,               'p'},
		{"prefix",             1, NULL,               'x'},
		{"threshold",          1, NULL,               't'},
		{"beam",               1, NULL,               'b'},

		/* Long-only options with no arguments */
		{"filter-noise",       0, &iargs.noisefilter,        1},
		{"no-check-prefix",    0, &config_checkprefix,       0},
		{"basename",           0, &config_basename,          1},
		{"no-peaks-in-stream", 0, &iargs.stream_peaks,       0},
		{"no-refls-in-stream", 0, &iargs.stream_refls,       0},
		{"integrate-saturated",0, &integrate_saturated,      1},
		{"no-use-saturated",   0, &iargs.use_saturated,      0},
		{"no-revalidate",      0, &iargs.no_revalidate,      1},
		{"check-hdf5-snr",     0, &iargs.check_hdf5_snr,     1},
		{"no-refine",          0, &iargs.predict_refine,     0},

		/* Long-only options which don't actually do anything */
		{"no-sat-corr",        0, &iargs.satcorr,            0},
		{"sat-corr",           0, &iargs.satcorr,            1},
		{"no-check-hdf5-snr",  0, &iargs.check_hdf5_snr,     0},
		{"use-saturated",      0, &iargs.use_saturated,      1},

		/* Long-only options with arguments */
		{"peaks",              1, NULL,                2},
		{"cell-reduction",     1, NULL,                3},
		{"min-gradient",       1, NULL,                4},
		{"record",             1, NULL,                5},
		{"cpus",               1, NULL,                6},
		{"cpugroup",           1, NULL,                7},
		{"cpuoffset",          1, NULL,                8},
		{"hdf5-peaks",         1, NULL,                9},
		{"copy-hdf5-field",    1, NULL,               10},
		{"min-snr",            1, NULL,               11},
		{"tolerance",          1, NULL,               13},
		{"int-radius",         1, NULL,               14},
		{"median-filter",      1, NULL,               15},
		{"integration",        1, NULL,               16},
		{"temp-dir",           1, NULL,               17},
		{"int-diag",           1, NULL,               18},
		{"push-res",           1, NULL,               19},
		{"res-push",           1, NULL,               19}, /* compat */
		{"peak-radius",        1, NULL,               20},
		{"highres",            1, NULL,               21},
		{"fix-profile-radius", 1, NULL,               22},
		{"fix-bandwidth",      1, NULL,               23},
		{"fix-divergence",     1, NULL,               24},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:z:p:x:j:g:t:vb:",
	                        longopts, NULL)) != -1)
	{
		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'v' :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'b' :
			ERROR("WARNING: This version of CrystFEL no longer "
			      "uses beam files.  Please remove the beam file "
			      "from your indexamajig command line.\n");
			return 1;

			case 'i' :
			filename = strdup(optarg);
			break;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 'z' :
			indm_str = strdup(optarg);
			break;

			case 'p' :
			cellfile = strdup(optarg);
			break;

			case 'x' :
			prefix = strdup(optarg);
			break;

			case 'j' :
			n_proc = atoi(optarg);
			break;

			case 'g' :
			geom_filename = optarg;
			break;

			case 't' :
			iargs.threshold = strtof(optarg, NULL);
			break;

			case 2 :
			speaks = strdup(optarg);
			break;

			case 3 :
			ERROR("The option '--cell-reduction' is no longer "
			      "used.\n"
			      "The complete indexing behaviour is now "
			      "controlled using '--indexing'.\n"
			      "See 'man indexamajig' for details of the "
			      "available methods.\n");
			return 1;

			case 4 :
			iargs.min_gradient = strtof(optarg, NULL);
			break;

			case 5 :
			ERROR("The option '--record' is no longer used.\n"
			      "Use '--no-peaks-in-stream' and"
			      "'--no-refls-in-stream' if you need to control"
			      "the contents of the stream.\n");
			return 1;

			case 6 :
			case 7 :
			case 8 :
			ERROR("The options --cpus, --cpugroup and --cpuoffset"
			      " are no longer used by indexamajig.\n");
			break;

			case 9 :
			free(iargs.hdf5_peak_path);
			iargs.hdf5_peak_path = strdup(optarg);
			break;

			case 10 :
			add_copy_hdf5_field(iargs.copyme, optarg);
			break;

			case 11 :
			iargs.min_snr = strtof(optarg, NULL);
			break;

			case 13 :
			toler = strdup(optarg);
			break;

			case 14 :
			intrad = strdup(optarg);
			break;

			case 15 :
			iargs.median_filter = atoi(optarg);
			break;

			case 16 :
			int_str = strdup(optarg);
			break;

			case 17 :
			tempdir = strdup(optarg);
			break;

			case 18 :
			int_diag = strdup(optarg);
			break;

			case 19 :
			if ( sscanf(optarg, "%f", &iargs.push_res) != 1 ) {
				ERROR("Invalid value for --push-res\n");
				return 1;
			}
			iargs.push_res *= 1e9;  /* nm^-1 -> m^-1 */
			have_push_res = 1;
			break;

			case 20 :
			pkrad = strdup(optarg);
			break;

			case 21 :
			if ( sscanf(optarg, "%f", &iargs.highres) != 1 ) {
				ERROR("Invalid value for --highres\n");
				return 1;
			}
			/* A -> m^-1 */
			iargs.highres = 1.0 / (iargs.highres/1e10);
			break;

			case 22 :
			if ( sscanf(optarg, "%f", &iargs.fix_profile_r) != 1 ) {
				ERROR("Invalid value for "
				      "--fix-profile-radius\n");
				return 1;
			}
			break;

			case 23 :
			if ( sscanf(optarg, "%f", &iargs.fix_bandwidth) != 1 ) {
				ERROR("Invalid value for --fix-bandwidth\n");
				return 1;
			}
			break;

			case 24 :
			if ( sscanf(optarg, "%f", &iargs.fix_divergence) != 1 ) {
				ERROR("Invalid value for --fix-divergence\n");
				return 1;
			}
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( tempdir == NULL ) {
		tempdir = strdup(".");
	}

	if ( filename == NULL ) {
		filename = strdup("-");
	}
	if ( strcmp(filename, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(filename, "r");
	}
	if ( fh == NULL ) {
		ERROR("Failed to open input file '%s'\n", filename);
		return 1;
	}
	free(filename);

	if ( speaks == NULL ) {
		speaks = strdup("zaef");
		STATUS("You didn't specify a peak detection method.\n");
		STATUS("I'm using 'zaef' for you.\n");
	}
	if ( strcmp(speaks, "zaef") == 0 ) {
		iargs.peaks = PEAK_ZAEF;
	} else if ( strcmp(speaks, "hdf5") == 0 ) {
		iargs.peaks = PEAK_HDF5;
	} else if ( strcmp(speaks, "cxi") == 0 ) {
		iargs.peaks = PEAK_CXI;
	} else {
		ERROR("Unrecognised peak detection method '%s'\n", speaks);
		return 1;
	}
	free(speaks);

	/* Set default path for peaks, if appropriate */
	if ( iargs.hdf5_peak_path == NULL ) {
		if ( iargs.peaks == PEAK_HDF5 ) {
			iargs.hdf5_peak_path = strdup("/processing/hitfinder/peakinfo");
		} else if ( iargs.peaks == PEAK_CXI ) {
			iargs.hdf5_peak_path = strdup("/entry_1/result_1");
		}
	}

	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
	}

	if ( n_proc == 0 ) {
		ERROR("Invalid number of processes.\n");
		return 1;
	}

	iargs.det = get_detector_geometry(geom_filename, iargs.beam);
	if ( iargs.det == NULL ) {
		ERROR("Failed to read detector geometry from  '%s'\n",
		      geom_filename);
		return 1;
	}

	if ( indm_str == NULL ) {

		STATUS("You didn't specify an indexing method, so I  won't try "
		       " to index anything.\n"
		       "If that isn't what you wanted, re-run with"
		       " --indexing=<methods>.\n");
		indm = NULL;

	} else {

		indm = build_indexer_list(indm_str);
		if ( indm == NULL ) {
			ERROR("Invalid indexer list '%s'\n", indm_str);
			return 1;
		}
		free(indm_str);
	}

	if ( int_str != NULL ) {

		int err;

		iargs.int_meth = integration_method(int_str, &err);
		if ( err ) {
			ERROR("Invalid integration method '%s'\n", int_str);
			return 1;
		}
		free(int_str);
	}
	if ( integrate_saturated ) {
		/* Option provided for backwards compatibility */
		iargs.int_meth |= INTEGRATION_SATURATED;
	}

	if ( have_push_res && !(iargs.int_meth & INTEGRATION_RESCUT) ) {
		ERROR("WARNING: You used --push-res, but not -rescut, "
		      "therefore --push-res will have no effect.\n");
	}

	if ( toler != NULL ) {
		int ttt;
		ttt = sscanf(toler, "%f,%f,%f,%f",
		             &iargs.tols[0], &iargs.tols[1],
		             &iargs.tols[2], &iargs.tols[3]);
		if ( ttt != 4 ) {
			ERROR("Invalid parameters for '--tolerance'\n");
			return 1;
		}
		free(toler);
	}

	if ( intrad != NULL ) {
		int r;
		r = sscanf(intrad, "%f,%f,%f",
		           &iargs.ir_inn, &iargs.ir_mid, &iargs.ir_out);
		if ( r != 3 ) {
			ERROR("Invalid parameters for '--int-radius'\n");
			return 1;
		}
		free(intrad);
	} else {
		STATUS("WARNING: You did not specify --int-radius.\n");
		STATUS("WARNING: I will use the default values, which are"
		       " probably not appropriate for your patterns.\n");
	}

	if ( pkrad != NULL ) {
		int r;
		r = sscanf(pkrad, "%f,%f,%f",
		           &iargs.pk_inn, &iargs.pk_mid, &iargs.pk_out);
		if ( r != 3 ) {
			ERROR("Invalid parameters for '--peak-radius'\n");
			return 1;
		}
		free(pkrad);
	}

	if ( iargs.pk_inn < 0.0 ) {
		iargs.pk_inn = iargs.ir_inn;
		iargs.pk_mid = iargs.ir_mid;
		iargs.pk_out = iargs.ir_out;
	}

	if ( iargs.det == NULL ) {
		ERROR("You need to provide a geometry file (please read the"
		      " manual for more details).\n");
		return 1;
	}

	add_geom_beam_stuff_to_copy_hdf5(iargs.copyme, iargs.det, iargs.beam);

	if ( cellfile != NULL ) {
		iargs.cell = load_cell_from_file(cellfile);
		if ( iargs.cell == NULL ) {
			ERROR("Couldn't read unit cell (from %s)\n", cellfile);
			return 1;
		}
		free(cellfile);
		STATUS("This is what I understood your unit cell to be:\n");
		cell_print(iargs.cell);
	} else {
		STATUS("No unit cell given.\n");
		iargs.cell = NULL;
	}

	if ( int_diag != NULL ) {

		int r;
		signed int h, k, l;

		if ( strcmp(int_diag, "random") == 0 ) {
			iargs.int_diag = INTDIAG_RANDOM;
		}

		if ( strcmp(int_diag, "all") == 0 ) {
			iargs.int_diag = INTDIAG_ALL;
		}

		if ( strcmp(int_diag, "negative") == 0 ) {
			iargs.int_diag = INTDIAG_NEGATIVE;
		}

		if ( strcmp(int_diag, "implausible") == 0 ) {
			iargs.int_diag = INTDIAG_IMPLAUSIBLE;
		}

		if ( strcmp(int_diag, "strong") == 0 ) {
			iargs.int_diag = INTDIAG_STRONG;
		}

		r = sscanf(int_diag, "%i,%i,%i", &h, &k, &l);
		if ( r == 3 ) {
			iargs.int_diag = INTDIAG_INDICES;
			iargs.int_diag_h = h;
			iargs.int_diag_k = k;
			iargs.int_diag_l = l;
		}

		if ( (iargs.int_diag == INTDIAG_NONE)
		  && (strcmp(int_diag, "none") != 0) ) {
			ERROR("Invalid value for --int-diag.\n");
			return 1;
		}

		free(int_diag);

		if ( (n_proc > 1) && (iargs.int_diag != INTDIAG_NONE) ) {
			n_proc = 1;
			STATUS("Ignored \"-j\" because you used --int-diag.\n");
		}

	}

	st = open_stream_for_write_2(outfile, geom_filename, argc, argv);
	if ( st == NULL ) {
		ERROR("Failed to open stream '%s'\n", outfile);
		return 1;
	}
	free(outfile);

	/* Prepare the indexer */
	if ( indm != NULL ) {
		ipriv = prepare_indexing(indm, iargs.cell, iargs.det,
		                         iargs.tols);
		if ( ipriv == NULL ) {
			ERROR("Failed to prepare indexing.\n");
			return 1;
		}
	} else {
		ipriv = NULL;
	}

	gsl_set_error_handler_off();

	iargs.indm = indm;
	iargs.ipriv = ipriv;

	create_sandbox(&iargs, n_proc, prefix, config_basename, fh,
	               st, tempdir);

	free(prefix);
	free(tempdir);
	free_detector_geometry(iargs.det);
	close_stream(st);
	cleanup_indexing(indm, ipriv);

	return 0;
}
