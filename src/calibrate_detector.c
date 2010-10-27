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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>

#include "utils.h"
#include "hdf5-file.h"
#include "filters.h"
#include "peaks.h"
#include "thread-pool.h"


#define INTEGRATION_RADIUS (10)


typedef enum
{
	SUM_THRESHOLD,
	SUM_PEAKS
} SumMethod;

struct sum_args
{
	char *filename;
	int config_cmfilter;
	int config_noisefilter;
	double *sum;
	int w;
	int h;
	SumMethod sum_method;
	double threshold;
};


struct queue_args
{
	FILE *fh;
	char *prefix;
	int config_cmfilter;
	int config_noisefilter;
	double *sum;
	int w;
	int h;
	SumMethod sum_method;
	double threshold;
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
"  -t, --threshold=<n>     Set the threshold if summing using the 'threshold'\n"
"                           method.\n"
"\n"
"      --filter-cm         Perform common-mode noise subtraction on images\n"
"                           before proceeding.\n"
"      --filter-noise      Apply an aggressive noise filter which sets all\n"
"                           pixels in each 3x3 region to zero if any of them\n"
"                           have negative values.\n"
"\n"
" -j <n>                   Run <n> analyses in parallel.  Default 1.\n"
" -x, --prefix=<p>         Prefix filenames from input file with 'p'.\n");
}


static void sum_peaks(struct image *image, double *sum)
{
	int x, y, i;
	int w = image->width;
	const int lim = INTEGRATION_RADIUS * INTEGRATION_RADIUS;

	/* FIXME: Get threshold value from command line */
	search_peaks(image, 800.0);

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


static void sum_threshold(struct image *image, double *sum, double threshold)
{
	int x, y;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
		float val = image->data[x+image->width*y];
		if ( val > threshold ) {
			sum[x+image->width*y] += val;
		}
	}
	}
}


static void add_image(void *args, int cookie)
{
	struct sum_args *pargs = args;
	struct hdfile *hdfile;
	struct image image;

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.indexed_cell = NULL;
	image.filename = pargs->filename;
	image.cpeaks = NULL;
	image.n_cpeaks = 0;
	image.det = NULL;

	/* View head-on (unit cell is tilted) */
	image.orientation.w = 1.0;
	image.orientation.x = 0.0;
	image.orientation.y = 0.0;
	image.orientation.z = 0.0;

	STATUS("%3i: Processing '%s'\n", cookie, pargs->filename);

	hdfile = hdfile_open(pargs->filename);
	if ( hdfile == NULL ) {
		return;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	/* FIXME: Nominal photon energy */
	hdf5_read(hdfile, &image, 1, 2000.0);

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
		sum_threshold(&image, pargs->sum, pargs->threshold);
		break;

	case SUM_PEAKS :
		sum_peaks(&image, pargs->sum);
		break;

	}

out:
	free(image.data);
	image_feature_list_free(image.features);
	if ( image.flags != NULL ) free(image.flags);
	hdfile_close(hdfile);

	free(pargs->filename);
	free(pargs);
}


static void *get_image(void *qp)
{
	char line[1024];
	struct sum_args *pargs;
	char *rval;
	struct queue_args *qargs = qp;

	/* Get the next filename */
	rval = fgets(line, 1023, qargs->fh);
	if ( rval == NULL ) return NULL;

	pargs = malloc(sizeof(struct sum_args));

	pargs->w = qargs->w;
	pargs->h = qargs->h;
	pargs->sum_method = qargs->sum_method;
	pargs->threshold = qargs->threshold;
	pargs->config_cmfilter = qargs->config_cmfilter;
	pargs->config_noisefilter = qargs->config_noisefilter;
	pargs->sum = qargs->sum;

	chomp(line);
	pargs->filename = malloc(strlen(qargs->prefix) + strlen(line) + 1);
	snprintf(pargs->filename, 1023, "%s%s", qargs->prefix, line);

	return pargs;
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *outfile = NULL;
	FILE *fh;
	int n_images = 0;
	int config_cmfilter = 0;
	int config_noisefilter = 0;
	char *prefix = NULL;
	char *sum_str = NULL;
	char *intermediate = NULL;
	double threshold = 400.0;
	SumMethod sum;
	int nthreads = 1;
	struct queue_args qargs;
	int n_done;
	const int chunk_size = 1000;

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
		{"threshold",          1, NULL,               't'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:x:j:o:s:p:t:",
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

		case 't' :
			threshold = atof(optarg);
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

	if ( nthreads == 0  ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	qargs.w = 1024;  /* FIXME! */
	qargs.h = 1024;  /* FIXME! */
	qargs.sum_method = sum;
	qargs.threshold = threshold;
	qargs.config_cmfilter = config_cmfilter;
	qargs.config_noisefilter = config_noisefilter;
	qargs.sum = calloc(qargs.w*qargs.h, sizeof(double));
	qargs.prefix = prefix;
	qargs.fh = fh;

	do {

		n_done = run_threads(nthreads, add_image, get_image,
		                     (void *)&qargs, NULL, chunk_size);

		n_images += n_done;

		/* Write intermediate sum if requested */
		if ( (intermediate != NULL) && (n_done == chunk_size) ) {
			char outfile[1024];
			snprintf(outfile, 1023, "%s-%i.h5",
			         intermediate, n_images);
			hdf5_write(outfile, qargs.sum, qargs.w, qargs.h,
			           H5T_NATIVE_DOUBLE);
		}

	} while ( n_done == chunk_size );

	/* Write the final output */
	hdf5_write(outfile, qargs.sum, qargs.w, qargs.h, H5T_NATIVE_DOUBLE);

	free(qargs.sum);
	free(prefix);
	free(outfile);
	if ( intermediate != NULL ) free(intermediate);
	fclose(fh);

	STATUS("There were %i images.\n", n_images);

	return 0;
}
