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
#include "reflections.h"
#include "stream.h"
#include "geometry.h"


#define MAX_THREADS (256)

struct process_args
{
	struct image *image;

	/* Thread control */
	pthread_mutex_t control_mutex;  /* Protects the scary stuff below */
	int start;
	int finish;
	int done;

	/* Analysis parameters */
	const char *sym;
	ReflItemList *obs;
	double *i_full;
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
"  -o, --output=<filename>    Output filename.  Default: facetron.hkl.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -x, --prefix=<p>           Prefix filenames from input file with <p>.\n"
"      --basename             Remove the directory parts of the filenames.\n"
"      --no-check-prefix      Don't attempt to correct the --prefix.\n"
"  -y, --symmetry=<sym>       Merge according to symmetry <sym>.\n"
"  -n, --iterations=<n>       Run <n> cycles of post-refinement.\n"
"\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
}


static void refine_image(struct image *image, ReflItemList *obs, double *i_full)
{
	//struct hdfile *hdfile;

	image->features = NULL;
	image->data = NULL;
	image->flags = NULL;
	image->hits = NULL;
	image->n_hits = 0;

	/* View head-on (unit cell is tilted) */
	image->orientation.w = 1.0;
	image->orientation.x = 0.0;
	image->orientation.y = 0.0;
	image->orientation.z = 0.0;

	//hdfile = hdfile_open(pargs->filename);
	//if ( hdfile == NULL ) {
	//	return;
	//} else if ( hdfile_set_first_image(hdfile, "/") ) {
	//	ERROR("Couldn't select path\n");
	//	hdfile_close(hdfile);
	//	return;
	//}

	//hdf5_read(hdfile, &image, 1);

	free(image->data);
	if ( image->flags != NULL ) free(image->flags);
	//hdfile_close(hdfile);
}


static void *worker_thread(void *pargsv)
{
	struct process_args *pargs = pargsv;
	int finish;

	do {

		int wakeup;

		refine_image(pargs->image, pargs->obs, pargs->i_full);

		pthread_mutex_lock(&pargs->control_mutex);
		pargs->done = 1;
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


static void refine_all(struct image *images, int n_total_patterns,
                       struct detector *det, const char *sym,
                       ReflItemList *obs, double *i_full, int nthreads)
{
	pthread_t workers[MAX_THREADS];
	struct process_args *worker_args[MAX_THREADS];
	int worker_active[MAX_THREADS];
	int i;
	int n_done = 0;
	int n_started = 0;

	/* Initialise worker arguments with the unchanging data */
	for ( i=0; i<nthreads; i++ ) {

		worker_args[i] = malloc(sizeof(struct process_args));
		worker_active[i] = 0;
		pthread_mutex_init(&worker_args[i]->control_mutex, NULL);
		worker_args[i]->sym = sym;

	}

	/* Start threads off */
	for ( i=0; i<nthreads; i++ ) {

		struct process_args *pargs;
		int r;

		pargs = worker_args[i];
		pargs->image = &images[n_started++];

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

		for ( i=0; i<nthreads; i++ ) {

			struct process_args *pargs;
			int done;

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

			n_done++;
			progress_bar(n_done, n_total_patterns, "Refining");

			/* Start work on the next pattern */
			if ( n_started == n_total_patterns ) break;
			pargs->image = &images[n_started++];

			/* Wake the thread up ... */
			pthread_mutex_lock(&pargs->control_mutex);
			pargs->done = 0;
			pargs->start = 1;
			pthread_mutex_unlock(&pargs->control_mutex);

		}

	} while ( n_started < n_total_patterns );

	/* Join threads */
	for ( i=0; i<nthreads; i++ ) {

		if ( !worker_active[i] ) continue;

		/* Tell the thread to exit */
		struct process_args *pargs = worker_args[i];
		pthread_mutex_lock(&pargs->control_mutex);
		pargs->finish = 1;
		pthread_mutex_unlock(&pargs->control_mutex);

		/* Wait for it to join */
		pthread_join(workers[i], NULL);

		n_done++;
		progress_bar(n_done, n_total_patterns, "Refining");

	}
}


static double partiality(struct image *image, double x, double y)
{
	return 1.0;
}


static void estimate_full(struct image *images, int n_total_patterns,
                          struct detector *det, const char *sym,
                          ReflItemList *obs, double *i_full)
{
	int i;

	for ( i=0; i<n_total_patterns; i++ ) {

		struct reflhit *spots;
		struct image *image = &images[i];
		int j, n;

		/* Get the "spots" appearing in this pattern */
		spots = find_intersections(image, image->indexed_cell,
		                           image->div, image->bw, &n, 0);

		/* For each spot, estimate the partiality */
		for ( j=0; j<n; j++ ) {
			double p = partiality(image, spots[j].x, spots[j].y);
		}

		progress_bar(i, n_total_patterns-1, "Integrating");

	}
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *outfile = NULL;
	char *geomfile = NULL;
	char *prefix = NULL;
	char *sym = NULL;
	FILE *fh;
	int nthreads = 1;
	int config_basename = 0;
	int config_checkprefix = 1;
	struct detector *det;
	double *i_full;
	ReflItemList *obs;
	int i;
	int n_total_patterns;
	struct image *images;
	int n_iter = 10;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"geometry",           1, NULL,               'g'},
		{"prefix",             1, NULL,               'x'},
		{"basename",           0, &config_basename,    1},
		{"no-check-prefix",    0, &config_checkprefix, 0},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:x:j:y:o:",
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

		case 'y' :
			sym = strdup(optarg);
			break;

		case 'o' :
			outfile = strdup(optarg);
			break;

		case 'n' :
			n_iter = atoi(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	/* Sanitise input filename and open */
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

	/* Sanitise output filename */
	if ( outfile == NULL ) {
		outfile = strdup("facetron.hkl");
	}

	/* Sanitise prefix */
	if ( prefix == NULL ) {
		prefix = strdup("");
	} else {
		if ( config_checkprefix ) {
			prefix = check_prefix(prefix);
		}
	}

	/* Get detector geometry */
	det = get_detector_geometry(geomfile);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geomfile);
		return 1;
	}
	free(geomfile);

	/* Prepare for iteration */
	i_full = new_list_intensity();
	obs = new_items();

	n_total_patterns = count_patterns(fh);
	STATUS("There are %i patterns to process\n", n_total_patterns);

	images = malloc(n_total_patterns * sizeof(struct image));
	if ( images == NULL ) {
		ERROR("Couldn't allocate memory for images.\n");
		return 1;
	}

	/* Fill in what we know about the images so far */
	rewind(fh);
	for ( i=0; i<n_total_patterns; i++ ) {

		UnitCell *cell;
		char *filename;
		char fnamereal[1024];

		if ( find_chunk(fh, &cell, &filename) == 1 ) {
			ERROR("Couldn't get all of the filenames and cells"
			      " from the input stream.\n");
			return 1;
		}

		images[i].indexed_cell = cell;

		/* Mangle the filename now */
		if ( config_basename ) {
			char *tmp;
			tmp = strdup(basename(filename));
			free(filename);
			filename = tmp;
		}
		snprintf(fnamereal, 1023, "%s%s", prefix, filename);

		images[i].filename = fnamereal;
		images[i].div = 0.5e-3;
		images[i].bw = 0.001;

		free(filename);

		progress_bar(i, n_total_patterns-1, "Loading pattern data");

	}
	fclose(fh);
	free(prefix);

	/* Make initial estimates */
	estimate_full(images, n_total_patterns, det, sym, obs, i_full);

	/* Iterate */
	for ( i=0; i<n_iter; i++ ) {

		STATUS("Post refinement iteration %i of %i\n", i+1, n_iter);

		/* Refine the geometry of all patterns to get the best fit */
		refine_all(images, n_total_patterns, det, sym, obs, i_full,
		           nthreads);

		/* Re-estimate all the full intensities */
		estimate_full(images, n_total_patterns, det, sym, obs, i_full);

	}

	/* Output results */
	write_reflections(outfile, obs, i_full, NULL, NULL, NULL);

	/* Clean up */
	free(i_full);
	delete_items(obs);
	free(sym);
	free(outfile);

	return 0;
}
