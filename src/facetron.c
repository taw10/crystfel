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
#include <errno.h>

#include "image.h"
#include "cell.h"
#include "hdf5-file.h"
#include "detector.h"
#include "geometry.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Profile fitting for coherent nanocrystallography.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<filename>     Specify the name of the input stream.\n"
"                              Can be '-' for stdin.\n"
"  -g. --geometry=<file>   Get detector geometry from file.\n"
);
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


static void get_trial_cell(UnitCell *out, UnitCell *in, int i, double step)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(in, &asx, &asy, &asz, &bsx, &bsy, &bsz,
	                    &csx, &csy, &csz);

	switch ( i ) {
	case  0 : asx += step; break;
	case  1 : asx -= step; break;
	case  2 : asy += step; break;
	case  3 : asy -= step; break;
	case  4 : asz += step; break;
	case  5 : asz -= step; break;
	case  6 : bsx += step; break;
	case  7 : bsx -= step; break;
	case  8 : bsy += step; break;
	case  9 : bsy -= step; break;
	case 10 : bsz += step; break;
	case 11 : bsz -= step; break;
	case 12 : csx += step; break;
	case 13 : csx -= step; break;
	case 14 : csy += step; break;
	case 15 : csy -= step; break;
	case 16 : csz += step; break;
	case 17 : csz -= step; break;
	default : break;
	}

	cell_set_reciprocal(out, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);
}


static void try_refine(struct image *image, UnitCell *cell,
                       double da, double dw, double step)
{
	struct reflhit *hits;
	int np;
	double itot;
	UnitCell *trial_cell;
	int did_something;

	trial_cell = cell_new();

	hits = find_intersections(image, cell, da, dw, &np, 0);
	itot = integrate_all(image, hits, np);
	STATUS("%f\n", itot);

	do {

		int i;

		did_something = 0;

		for ( i=0; i<18; i++ ) {

			struct reflhit *trial_hits;
			double trial_itot;

			get_trial_cell(trial_cell, cell, i, step);

			trial_hits = find_intersections(image, trial_cell,
				                        da, dw, &np, 0);
			trial_itot = integrate_all(image, hits, np);

			if ( trial_itot > itot ) {

				double asx, asy, asz;
				double bsx, bsy, bsz;
				double csx, csy, csz;

				cell_get_reciprocal(trial_cell, &asx, &asy,
				                    &asz, &bsx, &bsy, &bsz,
				                    &csx, &csy, &csz);
				cell_set_reciprocal(cell, asx, asy, asz, bsx,
				                    bsy, bsz, csx, csy, csz);

				itot = trial_itot;
				free(hits);
				hits = trial_hits;

				did_something = 1;

			} else {

				free(trial_hits);

			}
		}

	} while ( did_something );

	free(trial_cell);
}


/* Predict peaks */
static void pre_refine(struct image *image, UnitCell *cell,
                       double *da, double *dw)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double istep, step;

	/* Start by changing parameters by 1% */
	cell_get_reciprocal(cell, &asx, &asy,
	                    &asz, &bsx, &bsy, &bsz,
	                    &csx, &csy, &csz);
	istep = (asx+asy+asz+bsx+bsy+bsz+csx+csy+csz)/900.0;

	for ( step=istep; step>istep/100.0; step-=istep/100.0 ) {
		try_refine(image, cell, *da, *dw, step);
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	FILE *fh;
	UnitCell *cell;
	char *filename;
	struct detector *det;
	char *geometry = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"geometry",           1, NULL,               'g'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:g:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			infile = strdup(optarg);
			break;

		case 'g' :
			geometry = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( infile == NULL ) infile = strdup("-");
	if ( strcmp(infile, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(infile, "r");
	}
	free(infile);
	if ( fh == NULL ) {
		ERROR("Couldn't open input stream '%s'\n", infile);
		return ENOENT;
	}

	det = get_detector_geometry(geometry);
	if ( det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);

	/* Loop over all successfully indexed patterns */
	while ( find_chunk(fh, &cell, &filename) == 0 ) {

		struct image image;
		struct hdfile *hdfile;
		double da, dw;
		int np;

		STATUS("Integrating intensities from '%s'\n", filename);

		image.det = det;

		hdfile = hdfile_open(filename);
		if ( hdfile == NULL ) {
			return ENOENT;
		} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
			ERROR("Couldn't select path\n");
			return ENOENT;
		}

		hdf5_read(hdfile, &image, 0);

		da = 5.0e-3;     /* Initial divergence */
		dw = 3.0/100.0;  /* Initial bandwidth */

		pre_refine(&image, cell, &da, &dw);

		find_intersections(&image, cell, da, dw, &np, 1);

		hdfile_close(hdfile);
		cell_free(cell);
		free(filename);
		free(image.data);
		free(image.flags);

	}

	free(det->panels);
	free(det);
	fclose(fh);

	return 0;
}
