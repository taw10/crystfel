/*
 * sum_stack.c
 *
 * Sum a stack of images (e.g. for detector calibration)
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
	double min_gradient;
};


struct queue_args
{
	FILE *fh;
	int config_cmfilter;
	int config_noisefilter;
	double *sum;
	int w;
	int h;
	SumMethod sum_method;
	double threshold;
	double min_gradient;

	char *use_this_one_instead;
	char *prefix;
	int config_basename;
};


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Sum FEL diffraction images.\n"
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
"                           method.  Default: 400 adu\n"
"      --min-gradient=<n>   Minimum gradient for Zaefferer peak search.\n"
"                            Default: 100,000.\n"
"\n"
"      --filter-cm         Perform common-mode noise subtraction on images\n"
"                           before proceeding.\n"
"      --filter-noise      Apply an aggressive noise filter which sets all\n"
"                           pixels in each 3x3 region to zero if any of them\n"
"                           have negative values.\n"
"\n"
"  -j <n>                  Run <n> analyses in parallel.  Default 1.\n"
"  -x, --prefix=<p>        Prefix filenames from input file with 'p'.\n");
}


static void sum_peaks(struct image *image, double *sum, double threshold,
                      double min_gradient)
{
	int x, y, i;
	int w = image->width;
	const int lim = INTEGRATION_RADIUS * INTEGRATION_RADIUS;

	search_peaks(image, threshold, min_gradient);

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
	image.det = NULL;

	STATUS("%3i: Processing '%s'\n", cookie, pargs->filename);

	hdfile = hdfile_open(pargs->filename);
	if ( hdfile == NULL ) {
		return;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	hdf5_read(hdfile, &image, 1);

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
		sum_peaks(&image, pargs->sum, pargs->threshold,
		          pargs->min_gradient);
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
	char *line;
	struct sum_args *pargs;
	char *rval;
	struct queue_args *qargs = qp;

	pargs = malloc(sizeof(struct sum_args));

	pargs->w = qargs->w;
	pargs->h = qargs->h;
	pargs->sum_method = qargs->sum_method;
	pargs->threshold = qargs->threshold;
	pargs->min_gradient = qargs->min_gradient;
	pargs->config_cmfilter = qargs->config_cmfilter;
	pargs->config_noisefilter = qargs->config_noisefilter;
	pargs->sum = qargs->sum;

	/* Get the next filename */
	if ( qargs->use_this_one_instead != NULL ) {

		line = qargs->use_this_one_instead;
		qargs->use_this_one_instead = NULL;

	} else {

		line = malloc(1024*sizeof(char));
		rval = fgets(line, 1023, qargs->fh);
		if ( rval == NULL ) {
			free(pargs);
			return NULL;
		}
		chomp(line);

	}

	if ( qargs->config_basename ) {
		char *tmp;
		tmp = safe_basename(line);
		free(line);
		line = tmp;
	}

	pargs->filename = malloc(strlen(qargs->prefix)+strlen(line)+1);
	snprintf(pargs->filename, 1023, "%s%s", qargs->prefix, line);
	free(line);

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
	float min_gradient = 100000.0;
	SumMethod sum;
	int nthreads = 1;
	struct queue_args qargs;
	int n_done;
	const int chunk_size = 1000;
	struct hdfile *hdfile;
	struct image image;
	char *prepare_line;
	char prepare_filename[1024];
	char *rval;
	int config_basename = 0;

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
		{"min-gradient",       1, NULL,                4},
		{"basename",           0, &config_basename,    1},
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

	/* Get first filename and use it to set up the summed array */
	prepare_line = malloc(1024*sizeof(char));
	rval = fgets(prepare_line, 1023, fh);
	if ( rval == NULL ) {
		ERROR("Failed to get filename to prepare indexing.\n");
		return 1;
	}
	chomp(prepare_line);
	if ( config_basename ) {
		char *tmp;
		tmp = safe_basename(prepare_line);
		free(prepare_line);
		prepare_line = tmp;
	}
	snprintf(prepare_filename, 1023, "%s%s", prefix, prepare_line);
	qargs.use_this_one_instead = prepare_line;

	hdfile = hdfile_open(prepare_filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open '%s'\n", prepare_filename);
		return 1;
	}

	if ( hdf5_read(hdfile, &image, 0) ) {
		ERROR("Couldn't read '%s'\n", prepare_filename);
		return 1;
	}

	hdfile_close(hdfile);
	qargs.w = image.width;
	qargs.h = image.height;

	qargs.sum_method = sum;
	qargs.threshold = threshold;
	qargs.config_basename = config_basename;
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
