/*
 * sim-main.c
 *
 * Simulate test data
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
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

#include "image.h"
#include "relrod.h"
#include "cell.h"
#include "utils.h"
#include "hdf5-file.h"


/* Crystal size in metres */
#define CRYSTAL_SIZE (500.0e-9)


static void main_show_help(const char *s)
{
	printf("Syntax: %s <file1.h5> <file2.h5> [...]\n\n", s);
	printf("Index diffraction patterns\n\n");
	printf("  -h              Display this help message\n");
}


int main(int argc, char *argv[])
{
	int c, i;
	ImageList *list;
	UnitCell *cell;
	struct image image;
	int nrefl;

	while ((c = getopt(argc, argv, "h")) != -1) {

		switch ( c ) {

			case 'h' : {
				main_show_help(argv[0]);
				return 0;
			}

		}

	}

	/* Define unit cell */
	cell = cell_new_from_parameters(28.10e-9,
	                                28.10e-9,
	                                16.52e-9,
	                                deg2rad(90.0),
	                                deg2rad(90.0),
	                                deg2rad(120.0));

	/* Define image parameters */
	list = image_list_new();
	image.width = 512;
	image.height = 512;
	image.tilt = deg2rad(32.0);
	image.omega = deg2rad(40.0);
	image.fmode = FORMULATION_CLEN;
	image.x_centre = 255.5;
	image.y_centre = 255.5;
	image.camera_len = 0.2;  /* 20 cm */
	image.resolution = 5120; /* 512 pixels in 10 cm */
	image.lambda = 0.2e-9;   /* LCLS wavelength */
	image.data = malloc(512*512*2);
	image_add(list, &image);

	/* Calculate reflections */
	get_reflections(&image, cell, 1.0/CRYSTAL_SIZE);

	/* Construct the image */
	nrefl = image_feature_count(image.rflist);
	for ( i=0; i<nrefl; i++ ) {

		struct imagefeature *f;
		int x, y;

		f = image_get_feature(image.rflist, i);

		x = f->x;
		y = f->y;	/* Discards digits after the decimal point */

		image.data[y*image.width+x] = 1;

	}

	/* Write the output file */
	hdf5_write("simulated.h5", image.data, image.width, image.height);

	return 0;
}
