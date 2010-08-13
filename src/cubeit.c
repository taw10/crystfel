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
#include <errno.h>

#include "image.h"
#include "cell.h"
#include "hdf5-file.h"


#define MAX_HITS (1024)


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


int main(int argc, char *argv[])
{
	int c;
	char *infile = NULL;
	FILE *fh;
	UnitCell *cell;
	char *filename;
	unsigned int angles[180];
	int i;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			infile = strdup(optarg);
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
	if ( fh == NULL ) {
		ERROR("Couldn't open input stream '%s'\n", infile);
		return ENOENT;
	}

	/* Initialise histogram */
	for ( i=0; i<180; i++ ) angles[i] = 0;

	/* Loop over all successfully indexed patterns */
	while ( find_chunk(fh, &cell, &filename) == 0 ) {

		struct image image;
		struct hdfile *hdfile;
		double ang;
		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		unsigned int bin;

		//STATUS("Processing '%s'\n", filename);

		#if 0
		hdfile = hdfile_open(filename);
		if ( hdfile == NULL ) {
			return ENOENT;
		} else if ( hdfile_set_image(hdfile, "/data/data0") ) {
			ERROR("Couldn't select path\n");
			return ENOENT;
		}

		hdf5_read(hdfile, &image, 0);
		#endif

		cell_get_reciprocal(cell, &asx, &asy, &asz,
		                          &bsx, &bsy, &bsz,
		                          &csx, &csy, &csz);
		ang = angle_between(csx, csy, csz, 0.0, 0.0, 1.0);  /* 0->pi */
		ang = rad2deg(ang);  /* 0->180 deg */
		bin = rint(ang);
		angles[bin]++;

	}

	for ( i=0; i<180; i++ ) {
		STATUS("%i %i\n", i, angles[i]);
	}

	return 0;
}
