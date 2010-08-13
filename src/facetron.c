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


#define MAX_HITS (1024)


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


/* Predict peaks */
static int find_intersections(struct image *image, UnitCell *cell,
                              double divergence, double bandwidth)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct reflhit *hits;
	int np = 0;
	int hmax, kmax, lmax;
	double mres;
	signed int h, k, l;

	hits = malloc(sizeof(struct reflhit)*MAX_HITS);
	if ( hits == NULL ) return 0;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	mres = 1.0 / 8.0e-10;  /* 8 Angstroms */
	hmax = mres / modulus(asx, asy, asz);
	kmax = mres / modulus(bsx, bsy, bsz);
	lmax = mres / modulus(csx, csy, csz);

	for ( h=-hmax; h<hmax; h++ ) {
	for ( k=-kmax; k<kmax; k++ ) {
	for ( l=-lmax; l<lmax; l++ ) {

		double xl, yl, zl;
		double ds_sq, dps_sq;
		double delta, divfact;
		double llow, lhigh;
		double xd, yd, cl;
		double xda, yda;
		int p;
		int found = 0;

		if ( (h==0) && (k==0) && (l==0) ) continue;

		llow = image->lambda - image->lambda*bandwidth/2.0;
		lhigh = image->lambda + image->lambda*bandwidth/2.0;

		/* Get the coordinates of the reciprocal lattice point */
		zl = h*asz + k*bsz + l*csz;
		if ( zl < 0.0 ) continue;  /* Do this check very early */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;

		ds_sq = modulus_squared(xl, yl, zl);  /* d*^2 */
		delta = divergence/image->lambda;
		dps_sq = ds_sq + pow(delta, 2.0);  /* d'*^2 */

		/* In range? */
		divfact = 2.0 * delta * sqrt(xl*xl + yl*yl);
		if ( ds_sq - 2.0*zl/llow > 0.0 ) continue;
		if ( ds_sq - 2.0*zl/lhigh < 0.0 ) continue;

		/* Work out which panel this peak would fall on */
		for ( p=0; p<image->det->n_panels; p++ ) {

			/* Camera length for this panel */
			cl = image->det->panels[p].clen;

			/* Coordinates of peak relative to central beam, in m */
			xd = cl*xl / (ds_sq/(2.0*zl) - zl);
			yd = cl*yl / (ds_sq/(2.0*zl) - zl);

			/* Convert to pixels */
			xd *= image->det->panels[p].res;
			yd *= image->det->panels[p].res;

			/* Add the coordinates of the central beam */
			xda = xd + image->det->panels[p].cx;
			yda = yd + image->det->panels[p].cy;

			/* Now, is this on this panel? */
			if ( xda < image->det->panels[p].min_x ) continue;
			if ( xda > image->det->panels[p].max_x ) continue;
			if ( yda < image->det->panels[p].min_y ) continue;
			if ( yda > image->det->panels[p].max_y ) continue;

			/* Woohoo! */
			found = 1;
			break;

		}

		if ( !found ) continue;

		hits[np].h = h;
		hits[np].k = k;
		hits[np].l = l;
		np++;
		printf("%i %i %i 0.0 (at %f,%f) %i\n", h, k, l, xda, yda, p);

	}
	}
	}

	free(hits);  /* FIXME! */

	return np;
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
		find_intersections(&image, cell, 5.0e-3, 3.0/100.0);

		hdfile_close(hdfile);
		free(cell);
		free(filename);
		free(image.data);
		free(image.flags);

	}

	free(det->panels);
	free(det);
	fclose(fh);

	return 0;
}
