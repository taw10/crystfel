/*
 * reintegrate.c
 *
 * Like "indexamajig", but skip the indexing step
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
#include <pthread.h>
#include <sys/time.h>
#include <assert.h>

#include "utils.h"
#include "hdf5-file.h"
#include "symmetry.h"
#include "peaks.h"
#include "stream.h"
#include "index.h"
#include "thread-pool.h"


struct static_integration_args
{
	struct detector *det;
	pthread_mutex_t *output_mutex;  /* Protects 'stdout' */
	int config_cmfilter;
	int config_polar;
	int config_satcorr;
	int config_sa;
	int config_closer;
	int config_sanity;
};


struct integration_args
{
	char *filename;
	UnitCell *cell;

	struct static_integration_args static_args;
};


struct queue_args
{
	FILE *fh;
	const char *prefix;
	int config_basename;

	struct static_integration_args static_args;
};


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Like 'indexamajig', but skip the indexing step.\n"
"\n"
"  -h, --help               Display this help message.\n"
"\n"
"  -i, --input=<filename>   Specify the name of the input 'stream'.\n"
"                            (must be a file, not e.g. stdin)\n"
"  -g. --geometry=<file>    Get detector geometry from file.\n"
"  -x, --prefix=<p>         Prefix filenames from input file with <p>.\n"
"      --basename           Remove the directory parts of the filenames.\n"
"      --no-check-prefix    Don't attempt to correct the --prefix.\n"
"      --check-sanity       Check that indexed locations approximately correspond\n"
"                            with detected peaks.\n"
"      --filter-cm          Perform common-mode noise subtraction on images\n"
"                            before proceeding.  Intensities will be extracted\n"
"                            from the image as it is after this processing.\n"
"      --unpolarized        Don't correct for the polarisation of the X-rays.\n"
"      --no-sat-corr        Don't correct values of saturated peaks using a\n"
"                            table included in the HDF5 file.\n"
"      --no-sa              Don't correct for the differing solid angles of\n"
"                            the pixels.\n"
"      --no-closer-peak     Don't integrate from the location of a nearby peak\n"
"                            instead of the position closest to the reciprocal\n"
"                            lattice point.\n"
"  -j <n>                   Run <n> analyses in parallel.\n");
}


static void process_image(void *pg, int cookie)
{
	struct integration_args *apargs = pg;
	struct static_integration_args *pargs = &apargs->static_args;
	struct hdfile *hdfile;
	struct image image;

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.indexed_cell = NULL;
	image.filename = apargs->filename;
	image.cpeaks = NULL;
	image.n_cpeaks = 0;
	image.det = pargs->det;

	STATUS("Processing '%s'\n", apargs->filename);

	hdfile = hdfile_open(apargs->filename);
	if ( hdfile == NULL ) {
		return;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	/* FIXME: Nominal photon energy */
	hdf5_read(hdfile, &image, pargs->config_satcorr, 2000.0);

	map_all_peaks(&image);

	/* Sanity check */
	if ( pargs->config_sanity
	  && !peak_sanity_check(&image, image.indexed_cell, 1, 0.006e9) ) {

		STATUS("Failed peak sanity check.\n");

	} else {

		output_intensities(&image, apargs->cell,
		                   pargs->output_mutex, pargs->config_polar,
		                   pargs->config_sa, pargs->config_closer,
		                   stdout, 0, 0.1);
	}

	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	hdfile_close(hdfile);

	free(apargs->filename);
	cell_free(apargs->cell);
	free(apargs);
}


static void *get_image(void *qp)
{
	struct integration_args *pargs;
	struct queue_args *qargs = qp;
	UnitCell *cell;
	char *filename;

	/* Get the next filename */
	if ( find_chunk(qargs->fh, &cell, &filename) ) {
		return NULL;
	}

	pargs = malloc(sizeof(struct integration_args));

	if ( qargs->config_basename ) {
		char *tmp;
		tmp = strdup(basename(filename));
		free(filename);
		filename = tmp;
	}

	memcpy(&pargs->static_args, &qargs->static_args,
	       sizeof(struct static_integration_args));

	pargs->cell = cell;
	pargs->filename = malloc(1024);
	snprintf(pargs->filename, 1023, "%s%s", qargs->prefix, filename);

	return pargs;
}


static void integrate_all(int nthreads, struct detector *det, FILE *fh,
                          int config_basename, const char *prefix,
                          int config_cmfilter, int config_polar,
                          int config_satcorr, int config_sa, int config_closer,
                          int config_sanity)
{
	struct queue_args qargs;
	pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

	/* Information required to choose the next image */
	qargs.fh = fh;
	qargs.prefix = prefix;
	qargs.config_basename = config_basename;

	/* Information for the task which does not vary */
	qargs.static_args.det = det;
	qargs.static_args.config_cmfilter = config_cmfilter;
	qargs.static_args.config_polar = config_polar;
	qargs.static_args.config_satcorr = config_satcorr;
	qargs.static_args.config_sa = config_sa;
	qargs.static_args.config_closer = config_closer;
	qargs.static_args.config_sanity = qargs.static_args.config_sanity;
	qargs.static_args.output_mutex = &output_mutex;

	run_threads(nthreads, process_image, get_image, NULL, &qargs, 0);
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *geomfile = NULL;
	FILE *fh;
	char *prefix = NULL;
	int nthreads = 1;
	int config_basename = 0;
	int config_checkprefix = 1;
	int config_closer = 1;
	int config_polar = 1;
	int config_sanity = 0;
	int config_satcorr = 1;
	int config_sa = 1;
	int config_cmfilter = 0;
	struct detector *det;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"geometry",           1, NULL,               'g'},
		{"prefix",             1, NULL,               'x'},
		{"basename",           0, &config_basename,    1},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{"no-closer-peak",     0, &config_closer,      0},
		{"unpolarized",        0, &config_polar,       0},
		{"check-sanity",       0, &config_sanity,      1},
		{"sat-corr",           0, &config_satcorr,     1},
		{"no-sat-corr",        0, &config_satcorr,     0},
		{"no-sa",              0, &config_sa,          0},
		{"filter-cm",          0, &config_cmfilter,    1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:x:j:",
	                        longopts, NULL)) != -1)
	{

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			infile = strdup(optarg);
			break;

		case 'g' :
			geomfile = strdup(optarg);
			break;

		case 'x' :
			prefix = strdup(optarg);
			break;

		case 'j' :
			nthreads = atoi(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( infile == NULL ) {
		infile = strdup("-");
	}
	if ( strcmp(infile, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(infile, "r");
	}
	if ( fh == NULL ) {
		ERROR("Failed to open input file '%s'\n", infile);
		return 1;
	}
	free(infile);

	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
	}

	det = get_detector_geometry(geomfile);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geomfile);
		return 1;
	}
	free(geomfile);

	rewind(fh);
	integrate_all(nthreads, det, fh, config_basename, prefix,
	              config_cmfilter, config_polar, config_satcorr, config_sa,
	              config_closer, config_sanity);

	fclose(fh);
	free(prefix);

	return 0;
}
