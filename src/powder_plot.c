/*
 * powder_plot.c
 *
 * Plot powder patterns
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

#include "utils.h"
#include "image.h"
#include "detector.h"
#include "index.h"
#include "hdf5-file.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.h5>\n\n", s);
	printf(
"Plot a powder pattern as a 1D graph using the detector geometry.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -g. --geometry=<file>      Get detector geometry from file.\n"
"  -i, --input=<file>         Input filename.\n"
"\n");
}


int main(int argc, char *argv[])
{
	int c;
	struct image image;
	int x, y;
	struct hdfile *hdfile;
	char *filename = NULL;
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

		case 0 :
			break;

		case 'i' :
			filename = strdup(optarg);
			break;

		case 'g' :
			geometry = strdup(optarg);
			break;

		default :
			return 1;
		}

	}

	if ( filename == NULL ) {
		ERROR("You must specify the input filename with -i\n");
		return 1;
	}

	if ( geometry == NULL ) {
		ERROR("You need to specify a geometry file with --geometry\n");
		return 1;
	}

	image.det = get_detector_geometry(geometry);
	if ( image.det == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n", geometry);
		return 1;
	}
	free(geometry);

	hdfile = hdfile_open(filename);
	hdfile_set_image(hdfile, "/data/data");
	hdf5_read(hdfile, &image, 1, 2000.0); /* FIXME: Hardcoded */

	for ( x=0; x<image.width; x++ ) {
	for ( y=0; y<image.height; y++ ) {

		double q;
		int intensity;
		struct rvec r;

		r = get_q(&image, x, y, NULL, 1.0/image.lambda);
		q = modulus(r.u, r.v, r.w);

		intensity = image.data[x + image.width*y];

		printf("%5i\t%5i\t%7.3f\t%7i\n", x, y, q/1.0e9, intensity);

	}
	}

	return 0;
}
