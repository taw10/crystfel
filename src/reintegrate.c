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


#define MAX_THREADS (256)

struct process_args
{
	char *filename;
	int id;

	/* Thread control */
	pthread_mutex_t control_mutex;  /* Protects the scary stuff below */
	int start;
	int finish;
	int done;

	UnitCell *cell;
	struct detector *det;
	pthread_mutex_t *output_mutex;  /* Protects 'stdout' */
	int config_cmfilter;
	int config_polar;
	int config_satcorr;
	int config_sa;
	int config_closer;
	int config_sanity;
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
"      --sat-corr           Correct values of saturated peaks using a table\n"
"                            included in the HDF5 file.\n"
"      --no-sa              Don't correct for the differing solid angles of\n"
"                            the pixels.\n"
"      --no-closer-peak     Don't integrate from the location of a nearby peak\n"
"                            instead of the position closest to the reciprocal\n"
"                            lattice point.\n"
"  -j <n>                   Run <n> analyses in parallel.\n");
}


static void process_image(struct process_args *pargs)
{
	struct hdfile *hdfile;
	struct image image;

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.indexed_cell = NULL;
	image.id = pargs->id;
	image.filename = pargs->filename;
	image.hits = NULL;
	image.n_hits = 0;
	image.det = pargs->det;

	/* View head-on (unit cell is tilted) */
	image.orientation.w = 1.0;
	image.orientation.x = 0.0;
	image.orientation.y = 0.0;
	image.orientation.z = 0.0;

	STATUS("Processing '%s'\n", pargs->filename);

	hdfile = hdfile_open(pargs->filename);
	if ( hdfile == NULL ) {
		return;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	hdf5_read(hdfile, &image, pargs->config_satcorr);

	map_all_peaks(&image);

	/* Sanity check */
	if ( pargs->config_sanity
	  && !peak_sanity_check(&image, image.indexed_cell, 1, 0.006e9) ) {

		STATUS("Failed peak sanity check.\n");

	} else {

		output_intensities(&image, pargs->cell,
		                   pargs->output_mutex, pargs->config_polar,
		                   pargs->config_sa, pargs->config_closer,
		                   0, 0.1);
	}

	free(image.data);
	cell_free(pargs->cell);
	if ( image.flags != NULL ) free(image.flags);
	hdfile_close(hdfile);
}


static void *worker_thread(void *pargsv)
{
	struct process_args *pargs = pargsv;
	int finish;

	do {

		int wakeup;

		process_image(pargs);

		pthread_mutex_lock(&pargs->control_mutex);
		pargs->done = 1;
		pargs->start = 0;
		pthread_mutex_unlock(&pargs->control_mutex);

		/* Go to sleep until told to exit or process next image */
		do {

			pthread_mutex_lock(&pargs->control_mutex);
			/* Either of these can result in the thread waking up */
			wakeup = pargs->start || pargs->finish;
			finish = pargs->finish;
			pthread_mutex_unlock(&pargs->control_mutex);
			usleep(20000);

		} while ( !wakeup );

	} while ( !pargs->finish );

	return NULL;
}


static void integrate_all(int nthreads, struct detector *det, FILE *fh,
                          int config_basename, const char *prefix,
                          int config_cmfilter, int config_polar,
                          int config_satcorr, int config_sa, int config_closer,
                          int config_sanity)
{
	pthread_t workers[MAX_THREADS];
	struct process_args *worker_args[MAX_THREADS];
	int worker_active[MAX_THREADS];
	int i;
	int rval;
	pthread_mutex_t output_mutex = PTHREAD_MUTEX_INITIALIZER;

	/* Initialise worker arguments */
	for ( i=0; i<nthreads; i++ ) {

		worker_args[i] = malloc(sizeof(struct process_args));
		worker_args[i]->filename = malloc(1024);
		worker_active[i] = 0;
		worker_args[i]->det = det;
		worker_args[i]->config_cmfilter = config_cmfilter;
		worker_args[i]->config_polar = config_polar;
		worker_args[i]->config_sanity = config_sanity;
		worker_args[i]->config_satcorr = config_satcorr;
		worker_args[i]->config_sa = config_sa;
		worker_args[i]->config_closer = config_closer;
		pthread_mutex_init(&worker_args[i]->control_mutex, NULL);
		worker_args[i]->output_mutex = &output_mutex;

	}

	/* Start threads off */
	for ( i=0; i<nthreads; i++ ) {

		struct process_args *pargs;
		int r;
		int rval;
		char *filename;
		UnitCell *cell;

		pargs = worker_args[i];

		/* Get the next filename */
		rval = find_chunk(fh, &cell, &filename);
		if ( rval == 1 ) break;
		if ( config_basename ) {
			char *tmp;
			tmp = strdup(basename(filename));
			free(filename);
			filename = tmp;
		}
		snprintf(pargs->filename, 1023, "%s%s",
		         prefix, filename);
		pargs->cell = cell;
		free(filename);

		pthread_mutex_lock(&pargs->control_mutex);
		pargs->done = 0;
		pargs->start = 1;
		pargs->finish = 0;
		pthread_mutex_unlock(&pargs->control_mutex);

		worker_active[i] = 1;
		r = pthread_create(&workers[i], NULL, worker_thread, pargs);
		if ( r != 0 ) {
			worker_active[i] = 0;
			ERROR("Couldn't start thread %i\n", i);
		}

	}

	/* Keep threads busy until the end of the data */
	do {

		int i;
		rval = 0;

		for ( i=0; i<nthreads; i++ ) {

			struct process_args *pargs;
			int done;
			char *filename;
			UnitCell *cell;

			/* Spend time working, not managing threads */
			usleep(100000);

			/* Are we using this thread record at all? */
			if ( !worker_active[i] ) continue;

			/* Has the thread finished yet? */
			pargs = worker_args[i];
			pthread_mutex_lock(&pargs->control_mutex);
			done = pargs->done;
			pthread_mutex_unlock(&pargs->control_mutex);
			if ( !done ) continue;

			/* Get the next filename */
			rval = find_chunk(fh, &cell, &filename);
			if ( rval == 1 ) break;
			if ( config_basename ) {
				char *tmp;
				tmp = strdup(basename(filename));
				free(filename);
				filename = tmp;
			}
			snprintf(pargs->filename, 1023, "%s%s",
			         prefix, filename);
			pargs->cell = cell;
			free(filename);

			/* Wake the thread up ... */
			pthread_mutex_lock(&pargs->control_mutex);
			pargs->done = 0;
			pargs->start = 1;
			pthread_mutex_unlock(&pargs->control_mutex);

		}

	} while ( rval == 0 );

	/* Join threads */
	for ( i=0; i<nthreads; i++ ) {

		if ( !worker_active[i] ) goto free;

		/* Tell the thread to exit */
		struct process_args *pargs = worker_args[i];
		pthread_mutex_lock(&pargs->control_mutex);
		pargs->finish = 1;
		pthread_mutex_unlock(&pargs->control_mutex);

		/* Wait for it to join */
		pthread_join(workers[i], NULL);

	free:
		if ( worker_args[i]->filename != NULL ) {
			free(worker_args[i]->filename);
		}

	}
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
	int config_satcorr = 0;
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
