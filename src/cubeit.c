/*
 * cubeit.c
 *
 * "Full integration" of diffraction data
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
#include "diffraction.h"


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
	pthread_mutex_t vals_mutex;  /* Protects "vals" */
	double *vals;
	int xs;
	int ys;
	int zs;
	int config_angles;
	pthread_mutex_t angles_mutex;  /* Protects "angles" */
	unsigned int *angles;

};


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"'Full integration' of diffraction data.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input stream.\n"
"                              Can be '-' for stdin.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -x, --prefix=<p>           Prefix filenames from input file with <p>.\n"
"      --basename             Remove the directory parts of the filenames.\n"
"      --no-check-prefix      Don't attempt to correct the --prefix.\n"
"  -j <n>                     Run <n> analyses in parallel.\n");
}


static void interpolate_linear(double *vals, double v,
                               int xs, int ys, int zs,
                               int xv, int yv, double zv)
{
	int k;
	double val1, val2;
	float f, c;

	c = (zv+0.5)*(float)zs;
	k = (int)c;
	f = c - (float)k;
	assert(f >= 0.0);
	assert(k < zs);

	val1 = v*(1.0-f);
	val2 = v*f;

	vals[xs*ys*k + xs*yv + xv] += val1;
	vals[xs*ys*(k+1) + xs*yv + xv] += val2;
}


static void interpolate_bilinear(double *vals, double v,
                                 int xs, int ys, int zs,
                                 int xv, double yv, double zv)
{
	int k;
	double val1, val2;
	float f, c;

	c = (yv+0.5)*(float)ys;
	k = (int)c;
	f = c - (float)k;
	assert(f >= 0.0);
	assert(k < ys);

	val1 = v*(1.0-f);
	val2 = v*f;

	interpolate_linear(vals, val1, xs, ys, zs, xv, k, zv);
	interpolate_linear(vals, val2, xs, ys, zs, xv, k+1, zv);
}


static void interpolate_onto_grid(double *vals, double v,
                                  int xs, int ys, int zs,
                                  double xv, double yv, double zv)
{
	int k;
	double val1, val2;
	float f, c;

	c = (xv+0.5)*(float)xs;
	k = (int)c;
	f = c - (float)k;
	assert(f >= 0.0);
	assert(k < xs);

	val1 = v*(1.0-f);
	val2 = v*f;

	interpolate_bilinear(vals, val1, xs, ys, zs, k, yv, zv);
	interpolate_bilinear(vals, val2, xs, ys, zs, k+1, yv, zv);
}


static void process_image(struct process_args *pargs)
{
	struct hdfile *hdfile;
	struct image image;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	int x, y;
	const signed int ht = 3;
	const signed int kt = 4;
	const signed int lt = 5;

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
		return;
	} else if ( hdfile_set_first_image(hdfile, "/") ) {
		ERROR("Couldn't select path\n");
		hdfile_close(hdfile);
		return;
	}

	hdf5_read(hdfile, &image, 1);

	cell_get_cartesian(pargs->cell, &ax, &ay, &az, &bx, &by,
	                                &bz, &cx, &cy, &cz);

	for ( x=0; x<image.width; x++ ) {
	for ( y=0; y<image.height; y++ ) {

		double hd, kd, ld;
		signed int h, k, l;
		double dh, dk, dl;
		struct rvec q;

		q = get_q(&image, x, y, 1, NULL, 1.0/image.lambda);

		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;

		h = (signed int)rint(hd);
		k = (signed int)rint(kd);
		l = (signed int)rint(ld);

		if ( (h!=ht) || (k!=kt) || (l!=lt) ) continue;

		dh = hd - h;
		dk = kd - k;
		dl = ld - l;

		double v = image.data[x+image.width*y];

		pthread_mutex_lock(&pargs->vals_mutex);
		interpolate_onto_grid(pargs->vals, v, pargs->xs, pargs->ys,
		                      pargs->zs, dh, dk, dl);
		pthread_mutex_unlock(&pargs->vals_mutex);

	}
	}

	if ( pargs->config_angles ) {

		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		double ang;
		int bin;

		cell_get_reciprocal(pargs->cell, &asx, &asy, &asz,
		                                 &bsx, &bsy, &bsz,
		                                 &csx, &csy, &csz);
		ang = angle_between(csx, csy, csz, 0.0, 0.0, 1.0);
		ang = rad2deg(ang);  /* 0->180 deg */
		bin = rint(ang);
		pthread_mutex_lock(&pargs->vals_mutex);
		pargs->angles[bin]++;
		pthread_mutex_unlock(&pargs->vals_mutex);

	}

	free(image.data);
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


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	char *geomfile = NULL;
	FILE *fh;
	char *rval = NULL;
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
	double *vals;
	const int gs = 16;
	unsigned int angles[180];
	int config_angles = 0;

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

	/* Initialise histogram */
	for ( i=0; i<180; i++ ) angles[i] = 0;

	/* Initialise shape transform array */
	vals = calloc(gs*gs*gs, sizeof(double));

	if ( (nthreads == 0) || (nthreads > MAX_THREADS) ) {
		ERROR("Invalid number of threads.\n");
		return 1;
	}

	/* Initialise worker arguments */
	for ( i=0; i<nthreads; i++ ) {

		worker_args[i] = malloc(sizeof(struct process_args));
		worker_args[i]->filename = malloc(1024);
		worker_active[i] = 0;
		worker_args[i]->xs = gs;
		worker_args[i]->ys = gs;
		worker_args[i]->zs = gs;
		worker_args[i]->config_angles = config_angles;
		worker_args[i]->vals = vals;
		worker_args[i]->angles = angles;
		pthread_mutex_init(&worker_args[i]->control_mutex, NULL);
		pthread_mutex_init(&worker_args[i]->vals_mutex, NULL);
		pthread_mutex_init(&worker_args[i]->angles_mutex, NULL);

	}

	n_images = 0;

	/* Start threads off */
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

			char line[1024];
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

			/* Get the next filename */
			rval = fgets(line, 1023, fh);
			if ( rval == NULL ) break;
			chomp(line);
			snprintf(pargs->filename, 1023, "%s%s", prefix, line);

			n_images++;

			STATUS("Done %i images\n", n_images);

			/* Wake the thread up ... */
			pthread_mutex_lock(&pargs->control_mutex);
			pargs->done = 0;
			pargs->start = 1;
			pthread_mutex_unlock(&pargs->control_mutex);

		}

	} while ( rval != NULL );

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
