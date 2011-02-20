/*
 * indexamajig.c
 *
 * Index patterns, output hkl+intensity etc.
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
#include <pthread.h>
#include <sys/time.h>

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "sfac.h"
#include "filters.h"
#include "reflections.h"
#include "thread-pool.h"
#include "beam-parameters.h"
#include "geometry.h"


enum {
	PEAK_ZAEF,
	PEAK_HDF5,
};


/* Information about the indexing process which is common to all patterns */
struct static_index_args
{
	UnitCell *cell;
	int config_cmfilter;
	int config_noisefilter;
	int config_dumpfound;
	int config_verbose;
	int config_nearbragg;
	int config_polar;
	int config_satcorr;
	int config_closer;
	int config_insane;
	float threshold;
	float min_gradient;
	struct detector *det;
	IndexingMethod *indm;
	IndexingPrivate **ipriv;
	int peaks;
	int cellr;
	struct beam_params *beam;

	/* Output stream */
	pthread_mutex_t *output_mutex;  /* Protects the output stream */
	FILE *ofh;
};


/* Information about the indexing process for one pattern */
struct index_args
{
	/* "Input" */
	char *filename;
	struct static_index_args static_args;

	/* "Output" */
	int indexable;
};


/* Information needed to choose the next task and dispatch it */
struct queue_args
{
	FILE *fh;
	char *prefix;
	struct static_index_args static_args;

	int n_indexable;

	char *use_this_one_instead;
};


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
" -o, --output=<filename>  Write indexed stream to this file. '-' for stdout.\n"
"\n"
"     --indexing=<methods> Use 'methods' for indexing.  Provide one or more\n"
"                           methods separated by commas.  Choose from:\n"
"                            none     : no indexing (default)\n"
"                            dirax    : invoke DirAx\n"
"                            mosflm   : invoke MOSFLM (DPS)\n"
"                            template : index by template matching\n"
" -g. --geometry=<file>    Get detector geometry from file.\n"
" -b, --beam=<file>        Get beam parameters from file (provides nominal\n"
"                           wavelength value if no per-shot value is found in\n"
"                           the HDF5 files.\n"
" -p, --pdb=<file>         PDB file from which to get the unit cell to match.\n"
"                           Default: 'molecule.pdb'.\n"
" -x, --prefix=<p>         Prefix filenames from input file with <p>.\n"
"     --peaks=<method>     Use 'method' for finding peaks.  Choose from:\n"
"                           zaef  : Use Zaefferer (2000) gradient detection.\n"
"                                    This is the default method.\n"
"                           hdf5  : Get from /processing/hitfinder/peakinfo\n"
"                                    in the HDF5 file.\n"
"\n"
"\nWith just the above options, this program does not do much of practical use."
"\nYou should also enable some of the following:\n\n"
"     --near-bragg         Output a list of reflection intensities to stdout.\n"
"                           When pixels with fractional indices within 0.1 of\n"
"                           integer values (the Bragg condition) are found,\n"
"                           the integral of pixels within a ten pixel radius\n"
"                           of the nearest-to-Bragg pixel will be reported as\n"
"                           the intensity.  The centroid of the pixels will\n"
"                           be given as the coordinates, as well as the h,k,l\n"
"                           (integer) indices of the reflection.  If a peak\n"
"                           was located by the initial peak search close to\n"
"                           the \"near Bragg\" location, its coordinates will\n"
"                           be taken as the centre instead.\n"
"     --dump-peaks         Write the results of the peak search to stdout.\n"
"                           The intensities in this list are from the\n"
"                           centroid/integration procedure.\n"
"\n"
"\nFor more control over the process, you might need:\n\n"
"     --cell-reduction=<m> Use <m> as the cell reduction method. Choose from:\n"
"                           none    : no matching, just use the raw cell.\n"
"                           reduce  : full cell reduction.\n"
"                           compare : match by at most changing the order of\n"
"                                     the indices.\n"
"     --filter-cm          Perform common-mode noise subtraction on images\n"
"                           before proceeding.  Intensities will be extracted\n"
"                           from the image as it is after this processing.\n"
"     --filter-noise       Apply an aggressive noise filter which sets all\n"
"                           pixels in each 3x3 region to zero if any of them\n"
"                           have negative values.  Intensity measurement will\n"
"                           be performed on the image as it was before this.\n"
"     --unpolarized        Don't correct for the polarisation of the X-rays.\n"
"     --no-sat-corr        Don't correct values of saturated peaks using a\n"
"                           table included in the HDF5 file.\n"
"     --threshold=<n>      Only accept peaks above <n> ADU.  Default: 800.\n"
"     --min-gradient=<n>   Minimum gradient for Zaefferer peak search.\n"
"                           Default: 100,000.\n"
"\n"
"\nOptions for greater performance or verbosity:\n\n"
"     --verbose            Be verbose about indexing.\n"
" -j <n>                   Run <n> analyses in parallel.  Default 1.\n"
"\n"
"\nOptions you probably won't need:\n\n"
"     --no-check-prefix    Don't attempt to correct the --prefix.\n"
"     --no-closer-peak     Don't integrate from the location of a nearby peak\n"
"                           instead of the position closest to the reciprocal\n"
"                           lattice point.\n"
"     --insane             Don't check that the reduced cell accounts for at\n"
"                           least 10%% of the located peaks.\n"
);
}


static void process_image(void *pp, int cookie)
{
	struct index_args *pargs = pp;
	struct hdfile *hdfile;
	struct image image;
	float *data_for_measurement;
	size_t data_size;
	char *filename = pargs->filename;
	UnitCell *cell = pargs->static_args.cell;
	int config_cmfilter = pargs->static_args.config_cmfilter;
	int config_noisefilter = pargs->static_args.config_noisefilter;
	int config_dumpfound = pargs->static_args.config_dumpfound;
	int config_verbose = pargs->static_args.config_verbose;
	int config_nearbragg = pargs->static_args.config_nearbragg;
	int config_polar = pargs->static_args.config_polar;
	IndexingMethod *indm = pargs->static_args.indm;
	const struct beam_params *beam = pargs->static_args.beam;

	image.features = NULL;
	image.data = NULL;
	image.indexed_cell = NULL;
	image.id = cookie;
	image.filename = filename;
	image.det = pargs->static_args.det;

	STATUS("Processing '%s'\n", image.filename);

	pargs->indexable = 0;

	hdfile = hdfile_open(filename);
	if ( hdfile == NULL ) {
		return;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		return;
	}

	hdf5_read(hdfile, &image, pargs->static_args.config_satcorr,
	          beam->photon_energy);

	if ( config_cmfilter ) {
		filter_cm(&image);
	}

	/* Take snapshot of image after CM subtraction but before
	 * the aggressive noise filter. */
	data_size = image.width*image.height*sizeof(float);
	data_for_measurement = malloc(data_size);

	if ( config_noisefilter ) {
		filter_noise(&image, data_for_measurement);
	} else {
		memcpy(data_for_measurement, image.data, data_size);
	}

	switch ( pargs->static_args.peaks )
	{
	case PEAK_HDF5 :
		/* Get peaks from HDF5 */
		if ( get_peaks(&image, hdfile) ) {
			ERROR("Failed to get peaks from HDF5 file.\n");
			return;
		}
		break;
	case PEAK_ZAEF :
		search_peaks(&image, pargs->static_args.threshold,
		             pargs->static_args.min_gradient);
		break;
	}

	/* Get rid of noise-filtered version at this point
	 * - it was strictly for the purposes of peak detection. */
	free(image.data);
	image.data = data_for_measurement;

	if ( config_dumpfound ) {
		dump_peaks(&image, pargs->static_args.ofh,
		           pargs->static_args.output_mutex);
	}

	/* Not indexing? Then there's nothing left to do. */
	if ( indm == NULL ) goto done;

	/* Calculate orientation matrix (by magic) */
	index_pattern(&image, cell, indm, pargs->static_args.cellr,
		      config_verbose, pargs->static_args.ipriv,
		      pargs->static_args.config_insane);

	/* No cell at this point?  Then we're done. */
	if ( image.indexed_cell == NULL ) goto done;
	pargs->indexable = 1;

	/* Measure intensities if requested */
	if ( config_nearbragg ) {

		RefList *reflections;

		image.div = beam->divergence;
		image.bw = beam->bandwidth;
		image.profile_radius = 0.0001e9;

		//reflections = find_intersections(&image, image.indexed_cell,
		//                                 0);
		reflections = find_projected_peaks(&image, image.indexed_cell,
		                                   0, 0.1);

		output_intensities(&image, image.indexed_cell, reflections,
		                   pargs->static_args.output_mutex,
		                   config_polar,
		                   pargs->static_args.config_closer,
		                   pargs->static_args.ofh);

		reflist_free(reflections);

	}

	/* Only free cell if found */
	cell_free(image.indexed_cell);

done:
	free(image.data);
	free(image.flags);
	image_feature_list_free(image.features);
	hdfile_close(hdfile);
}


static void *get_image(void *qp)
{
	char line[1024];
	struct index_args *pargs;
	char *rval;
	struct queue_args *qargs = qp;

	/* Initialise new task arguments */
	pargs = malloc(sizeof(struct index_args));
	memcpy(&pargs->static_args, &qargs->static_args,
	       sizeof(struct static_index_args));

	/* Get the next filename */
	if ( qargs->use_this_one_instead != NULL ) {

		pargs->filename = malloc(strlen(qargs->prefix) +
		                       strlen(qargs->use_this_one_instead) + 1);

		snprintf(pargs->filename, 1023, "%s%s", qargs->prefix,
		         qargs->use_this_one_instead);

		qargs->use_this_one_instead = NULL;

	} else {

		rval = fgets(line, 1023, qargs->fh);
		if ( rval == NULL ) {
			free(pargs);
			return NULL;
		}
		chomp(line);
		pargs->filename = malloc(strlen(qargs->prefix)+strlen(line)+1);
		snprintf(pargs->filename, 1023, "%s%s", qargs->prefix, line);

	}

	return pargs;
}


static void finalise_image(void *qp, void *pp)
{
	struct queue_args *qargs = qp;
	struct index_args *pargs = pp;

	qargs->n_indexable += pargs->indexable;

	free(pargs->filename);
	free(pargs);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *outfile = NULL;
	FILE *fh;
	FILE *ofh;
	char *rval = NULL;
	int n_images;
	int config_noindex = 0;
	int config_dumpfound = 0;
	int config_nearbragg = 0;
	int config_cmfilter = 0;
	int config_noisefilter = 0;
	int config_verbose = 0;
	int config_polar = 1;
	int config_satcorr = 1;
	int config_checkprefix = 1;
	int config_closer = 1;
	int config_insane = 0;
	float threshold = 800.0;
	float min_gradient = 100000.0;
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
	int cellr;
	int peaks;
	int nthreads = 1;
	int i;
	pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;
	char prepare_line[1024];
	char prepare_filename[1024];
	struct queue_args qargs;
	struct beam_params *beam = NULL;
	double nominal_photon_energy;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"no-index",           0, &config_noindex,     1},
		{"dump-peaks",         0, &config_dumpfound,   1},
		{"peaks",              1, NULL,                2},
		{"cell-reduction",     1, NULL,                3},
		{"near-bragg",         0, &config_nearbragg,   1},
		{"indexing",           1, NULL,               'z'},
		{"geometry",           1, NULL,               'g'},
		{"beam",               1, NULL,               'b'},
		{"filter-cm",          0, &config_cmfilter,    1},
		{"filter-noise",       0, &config_noisefilter, 1},
		{"verbose",            0, &config_verbose,     1},
		{"pdb",                1, NULL,               'p'},
		{"prefix",             1, NULL,               'x'},
		{"unpolarized",        0, &config_polar,       0},
		{"no-sat-corr",        0, &config_satcorr,     0},
		{"sat-corr",           0, &config_satcorr,     1}, /* Compat */
		{"threshold",          1, NULL,               't'},
		{"min-gradient",       1, NULL,                4},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{"no-closer-peak",     0, &config_closer,      0},
		{"insane",             0, &config_insane,      1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:wp:j:x:g:t:o:b:",
	                        longopts, NULL)) != -1) {

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
			nthreads = atoi(optarg);
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

		case 2 :
			speaks = strdup(optarg);
			break;

		case 3 :
			scellr = strdup(optarg);
			break;

		case 4 :
			min_gradient = strtof(optarg, NULL);
			break;

		case 0 :
			break;

		default :
			return 1;
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
		outfile = strdup("-");
	}
	if ( strcmp(outfile, "-") == 0 ) {
		ofh = stdout;
	} else {
		ofh = fopen(outfile, "w");
	}
	if ( ofh == NULL ) {
		ERROR("Failed to open output file '%s'\n", outfile);
		return 1;
	}
	free(outfile);

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

	if ( nthreads == 0 ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	if ( indm_str == NULL ) {
		STATUS("You didn't specify an indexing method, so I won't"
		       " try to index anything.\n"
		       "If that isn't what you wanted, re-run with"
		       " --indexing=<method>.\n");
		indm = NULL;
		indexer_needs_cell = 0;
	} else {
		indm = build_indexer_list(indm_str, &indexer_needs_cell);
		if ( indm == NULL ) {
			ERROR("Invalid indexer list '%s'\n", indm_str);
			return 1;
		}
		free(indm_str);
	}

	reduction_needs_cell = 0;
	if ( scellr == NULL ) {
		STATUS("You didn't specify a cell reduction method, so I'm"
		       " going to use 'reduce'.\n");
		cellr = CELLR_REDUCE;
		reduction_needs_cell = 1;
	} else if ( strcmp(scellr, "none") == 0 ) {
		cellr = CELLR_NONE;
	} else if ( strcmp(scellr, "reduce") == 0) {
		cellr = CELLR_REDUCE;
		reduction_needs_cell = 1;
	} else if ( strcmp(scellr, "compare") == 0) {
		cellr = CELLR_COMPARE;
		reduction_needs_cell = 1;
	} else {
		ERROR("Unrecognised cell reduction method '%s'\n", scellr);
		return 1;
	}
	free(scellr);  /* free(NULL) is OK. */

	/* No indexing -> no reduction */
	if ( indm == NULL ) reduction_needs_cell = 0;

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

	if ( reduction_needs_cell || indexer_needs_cell ) {
		cell = load_cell_from_pdb(pdb);
		if ( cell == NULL ) {
			ERROR("Couldn't read unit cell (from %s)\n", pdb);
			return 1;
		}
	} else {
		STATUS("No cell needed for these choices of indexing"
		       " and reduction.\n");
		cell = NULL;
	}
	free(pdb);

	/* Start by writing the entire command line to stdout */
	fprintf(ofh, "Command line:");
	for ( i=0; i<argc; i++ ) {
		fprintf(ofh, " %s", argv[i]);
	}
	fprintf(ofh, "\n");
	fflush(ofh);

	if ( beam != NULL ) {
		nominal_photon_energy = beam->photon_energy;
	} else {
		STATUS("No beam parameters file was given, so I'm taking the"
		       " nominal photon energy to be 2 keV.\n");
		nominal_photon_energy = 2000.0;
	}

	/* Get first filename and use it to set up the indexing */
	rval = fgets(prepare_line, 1023, fh);
	if ( rval == NULL ) {
		ERROR("Failed to get filename to prepare indexing.\n");
		return 1;
	}
	chomp(prepare_line);
	snprintf(prepare_filename, 1023, "%s%s", prefix, prepare_line);
	qargs.use_this_one_instead = prepare_line;

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

	qargs.static_args.cell = cell;
	qargs.static_args.config_cmfilter = config_cmfilter;
	qargs.static_args.config_noisefilter = config_noisefilter;
	qargs.static_args.config_dumpfound = config_dumpfound;
	qargs.static_args.config_verbose = config_verbose;
	qargs.static_args.config_nearbragg = config_nearbragg;
	qargs.static_args.config_polar = config_polar;
	qargs.static_args.config_satcorr = config_satcorr;
	qargs.static_args.config_closer = config_closer;
	qargs.static_args.config_insane = config_insane;
	qargs.static_args.cellr = cellr;
	qargs.static_args.threshold = threshold;
	qargs.static_args.min_gradient = min_gradient;
	qargs.static_args.det = det;
	qargs.static_args.indm = indm;
	qargs.static_args.ipriv = ipriv;
	qargs.static_args.peaks = peaks;
	qargs.static_args.output_mutex = &output_mutex;
	qargs.static_args.ofh = ofh;
	qargs.static_args.beam = beam;

	qargs.fh = fh;
	qargs.prefix = prefix;
	qargs.n_indexable = 0;

	n_images = run_threads(nthreads, process_image, get_image,
	                       finalise_image, &qargs, 0);

	cleanup_indexing(ipriv);

	free(indm);
	free(ipriv);
	free(prefix);
	free(det->panels);
	free(det);
	cell_free(cell);
	if ( fh != stdin ) fclose(fh);
	if ( ofh != stdout ) fclose(ofh);

	STATUS("There were %i images, of which %i could be indexed.\n",
	        n_images, qargs.n_indexable);

	return 0;
}
