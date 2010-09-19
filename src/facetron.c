/*
 * facetron.c
 *
 * Profile fitting for coherent nanocrystallography
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
	char *sym;
};


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Post-refinement and profile fitting for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input 'stream'.\n"
"                              (must be a file, not e.g. stdin)\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -x, --prefix=<p>           Prefix filenames from input file with <p>.\n"
"      --basename             Remove the directory parts of the filenames.\n"
"      --no-check-prefix      Don't attempt to correct the --prefix.\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
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

	hdf5_read(hdfile, &image, 1);


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


static UnitCell *read_orientation_matrix(FILE *fh)
{
	float u, v, w;
	struct rvec as, bs, cs;
	UnitCell *cell;
	char line[1024];

	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read a-star\n");
		return NULL;
	}
	as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "bstar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read b-star\n");
		return NULL;
	}
	bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "cstar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read c-star\n");
		return NULL;
	}
	cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
	cell = cell_new_from_axes(as, bs, cs);

	return cell;
}


static int find_chunk(FILE *fh, UnitCell **cell, char **filename)
{
	char line[1024];
	char *rval = NULL;

	do {

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;

		chomp(line);

		if ( strncmp(line, "Reflections from indexing", 25) != 0 ) {
			continue;
		}

		*filename = strdup(line+29);

		/* Skip two lines (while checking for errors) */
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;

		*cell = read_orientation_matrix(fh);
		if ( *cell == NULL ) {
			STATUS("Got filename but no cell for %s\n", *filename);
			continue;
		}

		return 0;

	} while ( rval != NULL );

	return 1;
}


static void add_to_mean(UnitCell *cell, double *ast, double *bst, double *cst,
                        double *alst, double *best, double *gast)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(cell, &asx, &asy, &asz, &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	*ast += modulus(asx, asy, asz);
	*bst += modulus(bsx, bsy, bsz);
	*cst += modulus(csx, csy, csz);
	*alst += angle_between(bsx, bsy, bsz, csx, csy, csz);
	*best += angle_between(asx, asy, asz, csx, csy, csz);
	*gast += angle_between(asx, asy, asz, bsx, bsy, bsz);
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *geomfile = NULL;
	FILE *fh;
	int rval;
	int n_images;
	char *prefix = NULL;
	int nthreads = 1;
	pthread_t workers[MAX_THREADS];
	struct process_args *worker_args[MAX_THREADS];
	int worker_active[MAX_THREADS];
	int config_basename = 0;
	int config_checkprefix = 1;
	struct detector *det;
	int i;
	char *sym = NULL;
	double as, bs, cs, als, bes, gas;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"geometry",           1, NULL,               'g'},
		{"prefix",             1, NULL,               'x'},
		{"basename",           0, &config_basename,    1},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:x:j:",
	                        longopts, NULL)) != -1) {

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

	sym = strdup("6/mmm");  /* FIXME: Should be on command line */

	as = 0.0; bs = 0.0; cs = 0.0; als = 0.0; bes = 0.0; gas = 0.0;

	/* Initialise worker arguments */
	for ( i=0; i<nthreads; i++ ) {

		worker_args[i] = malloc(sizeof(struct process_args));
		worker_args[i]->filename = malloc(1024);
		worker_active[i] = 0;
		worker_args[i]->det = det;
		pthread_mutex_init(&worker_args[i]->control_mutex, NULL);
		worker_args[i]->sym = sym;

	}

	n_images = 0;

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
		add_to_mean(cell, &as, &bs, &cs, &als, &bes, &gas);
		if ( config_basename ) {
			char *tmp;
			tmp = basename(filename);
			free(filename);
			filename = tmp;
		}
		snprintf(pargs->filename, 1023, "%s%s",
		         prefix, filename);
		pargs->cell = cell;
		free(filename);

		n_images++;

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
			add_to_mean(cell, &as, &bs, &cs, &als, &bes, &gas);
			if ( config_basename ) {
				char *tmp;
				tmp = basename(filename);
				free(filename);
				filename = tmp;
			}
			snprintf(pargs->filename, 1023, "%s%s",
			         prefix, filename);
			pargs->cell = cell;
			free(filename);

			n_images++;

			STATUS("Done %i images\n", n_images);

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

	fclose(fh);

	STATUS("There were %i images.\n", n_images);

	return 0;
}
