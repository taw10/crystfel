/*
 * calibrate-detector.c
 *
 * Detector calibration
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _GNU_SOURCE 1
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>
#include <sys/time.h>

#include "utils.h"
#include "hdf5-file.h"
#include "filters.h"
#include "peaks.h"


#define INTEGRATION_RADIUS (10)
#define MAX_THREADS (96)


typedef enum
{
	SUM_THRESHOLD,
	SUM_PEAKS
} SumMethod;

struct process_args
{
	char *filename;
	int id;
	int config_cmfilter;
	int config_noisefilter;
	double *sum;
	int w;
	int h;
	SumMethod sum_method;
};


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Calibrate detector geometry from FEL diffraction images.\n"
"\n"
"  -h, --help              Display this help message.\n"
"\n"
"  -i, --input=<filename>  Specify file containing list of images to process.\n"
"                           '-' means stdin, which is the default.\n"
"  -o, --output=<filename> Output filename for summed image in HDF5 format.\n"
"                           Default: summed.h5.\n"
"\n"
"  -p, --intermediate=<P>  Stem of filename for intermediate images.\n"
"                            The filename stem <p> will be postfixed with a\n"
"                            hyphen, the current number of patterns processed\n"
"                            and '.h5'.  Such a pattern will be saved after\n"
"                            every 1000 input patterns.\n"
"                            If this option is not specified, no intermediate\n"
"                            patterns will be saved.\n"
"\n"
"  -s, --sum=<method>      Use this method for summation.  Choose from:\n"
"                           peaks : sum 10px radius circles around peaks.\n"
"                           threshold : sum thresholded images.\n"
"\n"
"      --filter-cm         Perform common-mode noise subtraction on images\n"
"                           before proceeding.\n"
"      --filter-noise      Apply an aggressive noise filter which sets all\n"
"                           pixels in each 3x3 region to zero if any of them\n"
"                           have negative values.\n"
"\n"
"  -j <n>                  Run <n> analyses in parallel.  Default 1.\n");
}


static void sum_peaks(struct image *image, double *sum)
{
	int x, y, i;
	int w = image->width;
	const int lim = INTEGRATION_RADIUS * INTEGRATION_RADIUS;

	search_peaks(image);

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f = image_get_feature(image->features, i);
		int xp, yp;

		/* This is not an error. */
		if ( f == NULL ) continue;

		xp = f->x;
		yp = f->y;

		for ( x=-INTEGRATION_RADIUS; x<+INTEGRATION_RADIUS; x++ ) {
		for ( y=-INTEGRATION_RADIUS; y<+INTEGRATION_RADIUS; y++ ) {

			/* Circular mask */
			if ( x*x + y*y > lim ) continue;

			if ( ((x+xp)>=image->width) || ((x+xp)<0) ) continue;
			if ( ((y+yp)>=image->height) || ((y+yp)<0) ) continue;

			float val = image->data[(x+xp)+w*(y+yp)];
			sum[(x+xp)+w*(y+yp)] += val;

		}
		}

	}
}


static void sum_threshold(struct image *image, double *sum)
{
	int x, y;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
		float val = image->data[x+image->width*y];
		if ( val > 400.0 ) {
			sum[x+image->width*y] += val;
		}
	}
	}
}


static void *process_image(void *pargsv)
{
	struct process_args *pargs = pargsv;
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
	image.det = NULL;

	/* View head-on (unit cell is tilted) */
	image.orientation.w = 1.0;
	image.orientation.x = 0.0;
	image.orientation.y = 0.0;
	image.orientation.z = 0.0;

	STATUS("Processing '%s'\n", pargs->filename);

	hdfile = hdfile_open(pargs->filename);
	if ( hdfile == NULL ) {
		return NULL;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return NULL;
	}

	hdf5_read(hdfile, &image);

	if ( pargs->config_cmfilter ) {
		filter_cm(&image);
	}

	if ( pargs->config_noisefilter ) {
		filter_noise(&image, NULL);
	}

	if ( (pargs->w != image.width) || (pargs->h != image.height) ) {
		ERROR("Wrong image size.\n");
		goto out;
	}

	switch ( pargs->sum_method ) {

	case SUM_THRESHOLD :
		sum_threshold(&image, pargs->sum);
		break;

	case SUM_PEAKS :
		sum_peaks(&image, pargs->sum);
		break;

	}

out:
	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	hdfile_close(hdfile);

	return NULL;
}


static void dump_to_file(struct process_args *worker_args[], int nthreads,
                         int w, int h, int n, const char *stem)
{
	int i;
	double *total;
	char outfile[256];

	total = calloc(w*h, sizeof(double));

	/* Add the individual sums to the 0th sum */
	for ( i=0; i<nthreads; i++ ) {

		int x, y;

		for ( x=0; x<w; x++ ) {
		for ( y=0; y<h; y++ ) {
			double val = worker_args[i]->sum[x+w*y];
			total[x+w*y] += val;
		}
		}

	}

	snprintf(outfile, 255, "%s-%i.h5", stem, n);

	hdf5_write(outfile, total, w, h, H5T_NATIVE_DOUBLE);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *outfile = NULL;
	FILE *fh;
	char *rval = NULL;
	int n_images;
	int config_cmfilter = 0;
	int config_noisefilter = 0;
	char *prefix = NULL;
	char *sum_str = NULL;
	char *intermediate = NULL;
	SumMethod sum;
	int nthreads = 1;
	pthread_t workers[MAX_THREADS];
	struct process_args *worker_args[MAX_THREADS];
	int worker_active[MAX_THREADS];
	int i;
	const int w = 1024;  /* FIXME! */
	const int h = 1024;  /* FIXME! */

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"filter-cm",          0, &config_cmfilter,    1},
		{"filter-noise",       0, &config_noisefilter, 1},
		{"prefix",             1, NULL,               'x'},
		{"sum",                1, NULL,               's'},
		{"intermediate",       1, NULL,               'p'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:x:j:o:s:p:",
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

		case 'x' :
			prefix = strdup(optarg);
			break;

		case 'j' :
			nthreads = atoi(optarg);
			break;

		case 's' :
			sum_str = strdup(optarg);
			break;

		case 'p' :
			intermediate = strdup(optarg);
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

	if ( sum_str == NULL ) {
		STATUS("You didn't specify a summation method, so I'm using"
		       " the 'peaks' method, which gives the best results.\n");
		sum = SUM_PEAKS;
	} else if ( strcmp(sum_str, "peaks") == 0 ) {
		sum = SUM_PEAKS;
	} else if ( strcmp(sum_str, "threshold") == 0) {
		sum = SUM_THRESHOLD;
	} else {
		ERROR("Unrecognised summation method '%s'\n", sum_str);
		return 1;
	}
	free(sum_str);

	if ( prefix == NULL ) {
		prefix = strdup("");
	}

	if ( outfile == NULL ) {
		outfile = strdup("summed.h5");
	}

	if ( (nthreads == 0) || (nthreads > MAX_THREADS) ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	/* Initialise worker arguments */
	for ( i=0; i<nthreads; i++ ) {

		worker_args[i] = malloc(sizeof(struct process_args));
		worker_args[i]->filename = malloc(1024);
		worker_args[i]->sum = calloc(w*h, sizeof(double));
		worker_active[i] = 0;

		worker_args[i]->w = w;
		worker_args[i]->h = h;
		worker_args[i]->sum_method = sum;

	}

	n_images = 0;

	/* Initially, fire off the full number of threads */
	for ( i=0; i<nthreads; i++ ) {

		char line[1024];
		struct process_args *pargs;
		int r;

		pargs = worker_args[i];

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);
		snprintf(pargs->filename, 1023, "%s%s", prefix, line);

		n_images++;

		pargs->config_cmfilter = config_cmfilter;
		pargs->config_noisefilter = config_noisefilter;

		worker_active[i] = 1;
		r = pthread_create(&workers[i], NULL, process_image, pargs);
		if ( r != 0 ) {
			worker_active[i] = 0;
			ERROR("Couldn't start thread %i\n", i);
		}

	}

	/* Start new threads as old ones finish */
	do {

		int i;

		for ( i=0; i<nthreads; i++ ) {

			char line[1024];
			int r;
			struct timespec t;
			struct timeval tv;
			struct process_args *pargs;

			if ( !worker_active[i] ) continue;

			pargs = worker_args[i];

			gettimeofday(&tv, NULL);
			t.tv_sec = tv.tv_sec;
			t.tv_nsec = tv.tv_usec * 1000 + 20000;

			r = pthread_timedjoin_np(workers[i], NULL, &t);
			if ( r != 0 ) continue; /* Not ready yet */

			worker_active[i] = 0;

			rval = fgets(line, 1023, fh);
			if ( rval == NULL ) break;
			chomp(line);
			snprintf(pargs->filename, 1023, "%s%s", prefix, line);

			worker_active[i] = 1;
			r = pthread_create(&workers[i], NULL, process_image, pargs);
			if ( r != 0 ) {
				worker_active[i] = 0;
				ERROR("Couldn't start thread %i\n", i);
			}

			n_images++;
			STATUS("Done %i images\n", n_images);

			if ( n_images % 1000 == 0 ) {
				if ( intermediate != NULL ) {
					dump_to_file(worker_args, nthreads,
					              w, h, n_images,
					              intermediate);
				}
			}
		}

	} while ( rval != NULL );

	/* Catch all remaining threads */
	for ( i=0; i<nthreads; i++ ) {

		if ( !worker_active[i] ) goto free;

		pthread_join(workers[i], NULL);

		worker_active[i] = 0;

	free:
		if ( worker_args[i]->filename != NULL ) {
			free(worker_args[i]->filename);
		}

	}

	/* Add the individual sums to the 0th sum */
	for ( i=1; i<nthreads; i++ ) {

		int x, y;

		for ( x=0; x<w; x++ ) {
		for ( y=0; y<h; y++ ) {
			double val = worker_args[i]->sum[x+w*y];
			worker_args[0]->sum[x+w*y] += val;
		}
		}
		free(worker_args[i]->sum);
		free(worker_args[i]);

	}

	hdf5_write(outfile, worker_args[0]->sum, w, h, H5T_NATIVE_DOUBLE);

	free(worker_args[0]->sum);
	free(worker_args[0]);
	free(prefix);
	free(outfile);
	if ( intermediate != NULL ) free(intermediate);
	fclose(fh);

	STATUS("There were %i images.\n", n_images);

	return 0;
}
