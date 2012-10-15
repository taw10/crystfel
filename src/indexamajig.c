/*
 * indexamajig.c
 *
 * Index patterns, output hkl+intensity etc.
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
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

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

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
#include "cell-utils.h"

#include "im-sandbox.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Process and index FEL diffraction images.\n"
"\n"
" -h, --help               Display this help message.\n"
"\n"
" -i, --input=<filename>   Specify file containing list of images to process.\n"
"                           '-' means stdin, which is the default.\n"
" -o, --output=<filename>  Write output stream to this file. '-' for stdout.\n"
"                           Default: indexamajig.stream\n"
"\n"
"     --indexing=<methods> Use 'methods' for indexing.  Provide one or more\n"
"                           methods separated by commas.  Choose from:\n"
"                            none     : no indexing (default)\n"
"                            dirax    : invoke DirAx\n"
"                            mosflm   : invoke MOSFLM (DPS)\n"
"                            reax     : DPS algorithm with known unit cell\n"
" -g. --geometry=<file>    Get detector geometry from file.\n"
" -b, --beam=<file>        Get beam parameters from file (provides nominal\n"
"                           wavelength value if no per-shot value is found in\n"
"                           the HDF5 files.\n"
" -p, --pdb=<file>         PDB file from which to get the unit cell to match.\n"
"                           Default: 'molecule.pdb'.\n"
"     --basename           Remove the directory parts of the filenames.\n"
" -x, --prefix=<p>         Prefix filenames from input file with <p>.\n"
"     --peaks=<method>     Use 'method' for finding peaks.  Choose from:\n"
"                           zaef  : Use Zaefferer (2000) gradient detection.\n"
"                                    This is the default method.\n"
"                           hdf5  : Get from a table in HDF5 file.\n"
"     --hdf5-peaks=<p>     Find peaks table in HDF5 file here.\n"
"                           Default: /processing/hitfinder/peakinfo\n"
"\n\n"
"You can control what information is included in the output stream using\n"
"' --record=<flag1>,<flag2>,<flag3>' and so on.  Possible flags are:\n\n"
" integrated        Include a list of reflection intensities, produced by\n"
"                    integrating around predicted peak locations.\n"
"\n"
" peaks             Include peak locations and intensities from the peak\n"
"                    search.\n"
"\n"
" peaksifindexed    As 'peaks', but only if the pattern could be indexed.\n"
"\n"
" peaksifnotindexed As 'peaks', but only if the pattern could NOT be indexed.\n"
"\n\n"
"The default is '--record=integrated'.\n"
"\n\n"
"For more control over the process, you might need:\n\n"
"  --cell-reduction=<m>  Use <m> as the cell reduction method. Choose from:\n"
"                         none    : no matching, just use the raw cell.\n"
"                         reduce  : full cell reduction.\n"
"                         compare : match by at most changing the order of\n"
"                                   the indices.\n"
"                         compare_ab : compare 'a' and 'b' lengths only.\n"
"    --tolerance=<tol>   Set the tolerances for cell reduction.\n"
"                          Default: 5,5,5,1.5.\n"
"    --filter-cm         Perform common-mode noise subtraction on images\n"
"                         before proceeding.  Intensities will be extracted\n"
"                         from the image as it is after this processing.\n"
"    --filter-noise      Apply an aggressive noise filter which sets all\n"
"                         pixels in each 3x3 region to zero if any of them\n"
"                         have negative values.  Intensity measurement will\n"
"                         be performed on the image as it was before this.\n"
"    --no-sat-corr       Don't correct values of saturated peaks using a\n"
"                         table included in the HDF5 file.\n"
"    --threshold=<n>     Only accept peaks above <n> ADU.  Default: 800.\n"
"    --min-gradient=<n>  Minimum gradient for Zaefferer peak search.\n"
"                         Default: 100,000.\n"
"    --min-snr=<n>       Minimum signal-to-noise ratio for peaks.\n"
"                         Default: 5.\n"
"    --min-integration-snr=<n> Minimum signal-to-noise ratio for peaks\n"
"                         during integration. Default: -infinity.\n"
"    --int-radius=<r>    Set the integration radii.  Default: 4,5,7.\n"
"-e, --image=<element>   Use this image from the HDF5 file.\n"
"                          Example: /data/data0.\n"
"                          Default: The first one found.\n"
"\n"
"\nFor time-resolved stuff, you might want to use:\n\n"
"     --copy-hdf5-field <f>  Copy the value of field <f> into the stream. You\n"
"                             can use this option as many times as you need.\n"
"\n"
"\nOptions for greater performance or verbosity:\n\n"
"     --verbose            Be verbose about indexing.\n"
" -j <n>                   Run <n> analyses in parallel.  Default 1.\n"
"\n"
"\nOptions you probably won't need:\n\n"
"     --no-check-prefix    Don't attempt to correct the --prefix.\n"
"     --closer-peak        Don't integrate from the location of a nearby peak\n"
"                           instead of the predicted spot.  Don't use.\n"
"     --insane             Don't check that the reduced cell accounts for at\n"
"                           least 10%% of the located peaks.\n"
"     --no-bg-sub          Don't subtract local background estimates from\n"
"                           integrated intensities.\n"
);
}


static int parse_cell_reduction(const char *scellr, int *err,
                                int *reduction_needs_cell)
{
	*err = 0;
	if ( strcmp(scellr, "none") == 0 ) {
		*reduction_needs_cell = 0;
		return CELLR_NONE;
	} else if ( strcmp(scellr, "reduce") == 0) {
		*reduction_needs_cell = 1;
		return CELLR_REDUCE;
	} else if ( strcmp(scellr, "compare") == 0) {
		*reduction_needs_cell = 1;
		return CELLR_COMPARE;
	} else if ( strcmp(scellr, "compare_ab") == 0) {
		*reduction_needs_cell = 1;
		return CELLR_COMPARE_AB;
	} else {
		*err = 1;
		*reduction_needs_cell = 0;
		return CELLR_NONE;
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *outfile = NULL;
	FILE *fh;
	FILE *ofh;
	char *rval = NULL;
	int config_noindex = 0;
	int config_cmfilter = 0;
	int config_noisefilter = 0;
	int config_verbose = 0;
	int config_satcorr = 1;
	int config_checkprefix = 1;
	int config_closer = 0;
	int config_insane = 0;
	int config_bgsub = 1;
	int config_basename = 0;
	float threshold = 800.0;
	float min_gradient = 100000.0;
	float min_snr = 5.0;
	double min_int_snr = -INFINITY;
	struct detector *det;
	char *geometry = NULL;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	int indexer_needs_cell;
	int reduction_needs_cell;
	char *indm_str = NULL;
	UnitCell *cell;
	char *pdb = NULL;
	char *prefix = NULL;
	char *speaks = NULL;
	char *scellr = NULL;
	char *toler = NULL;
	float tols[4] = {5.0, 5.0, 5.0, 1.5}; /* a,b,c,angles (%,%,%,deg) */
	int cellr;
	int peaks;
	int n_proc = 1;
	char *prepare_line;
	char prepare_filename[1024];
	char *use_this_one_instead;
	struct index_args iargs;
	struct beam_params *beam = NULL;
	char *element = NULL;
	double nominal_photon_energy;
	int stream_flags = STREAM_INTEGRATED;
	char *hdf5_peak_path = NULL;
	struct copy_hdf5_field *copyme;
	char *intrad = NULL;
	float ir_inn = 4.0;
	float ir_mid = 5.0;
	float ir_out = 7.0;

	copyme = new_copy_hdf5_field_list();
	if ( copyme == NULL ) {
		ERROR("Couldn't allocate HDF5 field list.\n");
		return 1;
	}

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"no-index",           0, &config_noindex,     1},
		{"indexing",           1, NULL,               'z'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"filter-cm",          0, &config_cmfilter,    1},
		{"filter-noise",       0, &config_noisefilter, 1},
		{"verbose",            0, &config_verbose,     1},
		{"pdb",                1, NULL,               'p'},
		{"prefix",             1, NULL,               'x'},
		{"no-sat-corr",        0, &config_satcorr,     0},
		{"sat-corr",           0, &config_satcorr,     1}, /* Compat */
		{"threshold",          1, NULL,               't'},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{"no-closer-peak",     0, &config_closer,      0},
		{"closer-peak",        0, &config_closer,      1},
		{"insane",             0, &config_insane,      1},
		{"image",              1, NULL,               'e'},
		{"basename",           0, &config_basename,    1},
		{"bg-sub",             0, &config_bgsub,       1}, /* Compat */
		{"no-bg-sub",          0, &config_bgsub,       0},

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
		{"min-integration-snr",1, NULL,               12},
		{"tolerance",          1, NULL,               13},
		{"int-radius",         1, NULL,               14},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:o:z:p:x:j:g:t:b:e:",
	                        longopts, NULL)) != -1)
	{
		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

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
			pdb = strdup(optarg);
			break;

			case 'x' :
			prefix = strdup(optarg);
			break;

			case 'j' :
			n_proc = atoi(optarg);
			break;

			case 'g' :
			geometry = strdup(optarg);
			break;

			case 't' :
			threshold = strtof(optarg, NULL);
			break;

			case 'b' :
			beam = get_beam_parameters(optarg);
			if ( beam == NULL ) {
				ERROR("Failed to load beam parameters"
				      " from '%s'\n", optarg);
				return 1;
			}
			break;

			case 'e' :
			element = strdup(optarg);
			break;

			case 2 :
			speaks = strdup(optarg);
			break;

			case 3 :
			scellr = strdup(optarg);
			break;

			case 4 :
			min_gradient = strtof(optarg, NULL);
			break;

			case 5 :
			stream_flags = parse_stream_flags(optarg);
			if ( stream_flags < 0 ) return 1;
			break;

			case 6 :
			case 7 :
			case 8 :
			ERROR("The options --cpus, --cpugroup and --cpuoffset"
			      " are no longer used by indexamajig.\n");
			break;

			case 9 :
			hdf5_peak_path = strdup(optarg);
			break;

			case 10 :
			add_copy_hdf5_field(copyme, optarg);
			break;

			case 11 :
			min_snr = strtof(optarg, NULL);
			break;

			case 12 :
			min_int_snr = strtof(optarg, NULL);
			break;

			case 13 :
			toler = strdup(optarg);
			break;

			case 14 :
			intrad = strdup(optarg);
			break;

			case 0 :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

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

	if ( outfile == NULL ) {
		ofh = stdout;
	} else {
		ofh = fopen(outfile, "w");
		if ( ofh == NULL ) {
			ERROR("Failed to open output file '%s'\n", outfile);
			return 1;
		}
		free(outfile);
	}

	if ( hdf5_peak_path == NULL ) {
		hdf5_peak_path = strdup("/processing/hitfinder/peakinfo");
	}

	if ( speaks == NULL ) {
		speaks = strdup("zaef");
		STATUS("You didn't specify a peak detection method.\n");
		STATUS("I'm using 'zaef' for you.\n");
	}
	if ( strcmp(speaks, "zaef") == 0 ) {
		peaks = PEAK_ZAEF;
	} else if ( strcmp(speaks, "hdf5") == 0 ) {
		peaks = PEAK_HDF5;
	} else {
		ERROR("Unrecognised peak detection method '%s'\n", speaks);
		return 1;
	}
	free(speaks);

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
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

	if ( (indm_str == NULL) ||
	     ((indm_str != NULL) && (strcmp(indm_str, "none") == 0)) ) {
		STATUS("Not indexing anything.\n");
		indexer_needs_cell = 0;
		reduction_needs_cell = 0;
		indm = NULL;
		cellr = CELLR_NONE;
	} else {
		if ( indm_str == NULL ) {
			STATUS("You didn't specify an indexing method, so I "
			       " won't try to index anything.\n"
			       "If that isn't what you wanted, re-run with"
			       " --indexing=<method>.\n");
			indm = NULL;
			indexer_needs_cell = 0;
		} else {
			indm = build_indexer_list(indm_str,
			                          &indexer_needs_cell);
			if ( indm == NULL ) {
				ERROR("Invalid indexer list '%s'\n", indm_str);
				return 1;
			}
			free(indm_str);
		}

		reduction_needs_cell = 0;
		if ( scellr == NULL ) {
			STATUS("You didn't specify a cell reduction method, so"
			       " I'm going to use 'reduce'.\n");
			cellr = CELLR_REDUCE;
			reduction_needs_cell = 1;
		} else {
			int err;
			cellr = parse_cell_reduction(scellr, &err,
			                             &reduction_needs_cell);
			if ( err ) {
				ERROR("Unrecognised cell reduction '%s'\n",
			              scellr);
				return 1;
			}
			free(scellr);
		}
	}

	/* No indexing -> no reduction */
	if ( indm == NULL ) reduction_needs_cell = 0;

	if ( toler != NULL ) {
		int ttt;
		ttt = sscanf(toler, "%f,%f,%f,%f",
		             &tols[0], &tols[1], &tols[2], &tols[3] );
		if ( ttt != 4 ) {
			ERROR("Invalid parameters for '--tolerance'\n");
			return 1;
		}
		free(toler);
	}

	if ( intrad != NULL ) {
		int r;
		r = sscanf(intrad, "%f,%f,%f", &ir_inn, &ir_mid, &ir_out);
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

	if ( geometry == NULL ) {
		ERROR("You need to specify a geometry file with --geometry\n");
		return 1;
	}

	det = get_detector_geometry(geometry);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);

	cell = load_cell_from_pdb(pdb);
	if ( cell == NULL ) {
		ERROR("Couldn't read unit cell (from %s)\n", pdb);
		return 1;
	}
	free(pdb);
	cell_print(cell);

	write_stream_header(ofh, argc, argv);

	if ( beam != NULL ) {
		nominal_photon_energy = beam->photon_energy;
	} else {
		STATUS("No beam parameters file was given, so I'm taking the"
		       " nominal photon energy to be 2 keV.\n");
		ERROR("I'm also going to assume 1 ADU per photon, which is");
		ERROR(" almost certainly wrong.  Peak sigmas will be"
		      " incorrect.\n");
		nominal_photon_energy = 2000.0;
	}

	/* Get first filename and use it to set up the indexing */
	prepare_line = malloc(1024);
	rval = fgets(prepare_line, 1023, fh);
	if ( rval == NULL ) {
		ERROR("Failed to get filename to prepare indexing.\n");
		return 1;
	}
	use_this_one_instead = strdup(prepare_line);
	chomp(prepare_line);
	if ( config_basename ) {
		char *tmp;
		tmp = safe_basename(prepare_line);
		free(prepare_line);
		prepare_line = tmp;
	}
	snprintf(prepare_filename, 1023, "%s%s", prefix, prepare_line);
	free(prepare_line);

	/* Prepare the indexer */
	if ( indm != NULL ) {
		ipriv = prepare_indexing(indm, cell, prepare_filename, det,
		                         nominal_photon_energy);
		if ( ipriv == NULL ) {
			ERROR("Failed to prepare indexing.\n");
			return 1;
		}
	} else {
		ipriv = NULL;
	}

	gsl_set_error_handler_off();

	/* Static worker args */
	iargs.cell = cell;
	iargs.config_cmfilter = config_cmfilter;
	iargs.config_noisefilter = config_noisefilter;
	iargs.config_verbose = config_verbose;
	iargs.config_satcorr = config_satcorr;
	iargs.config_closer = config_closer;
	iargs.config_insane = config_insane;
	iargs.config_bgsub = config_bgsub;
	iargs.cellr = cellr;
	iargs.tols[0] = tols[0];
	iargs.tols[1] = tols[1];
	iargs.tols[2] = tols[2];
	iargs.tols[3] = tols[3];
	iargs.threshold = threshold;
	iargs.min_gradient = min_gradient;
	iargs.min_snr = min_snr;
	iargs.min_int_snr = min_int_snr;
	iargs.det = det;
	iargs.indm = indm;
	iargs.ipriv = ipriv;
	iargs.peaks = peaks;
	iargs.beam = beam;
	iargs.element = element;
	iargs.stream_flags = stream_flags;
	iargs.hdf5_peak_path = hdf5_peak_path;
	iargs.copyme = copyme;
	iargs.ir_inn = ir_inn;
	iargs.ir_mid = ir_mid;
	iargs.ir_out = ir_out;

	create_sandbox(&iargs, n_proc, prefix, config_basename, fh,
	                    use_this_one_instead, ofh);

	free(prefix);

	return 0;
}
