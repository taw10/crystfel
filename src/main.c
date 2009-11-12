/*
 * main.c
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
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
#include "diffraction.h"
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
	int c;
	UnitCell *cell;
	struct image image;

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
	image.width = 512;
	image.height = 512;
	image.omega = deg2rad(40.0);
	image.tilt = deg2rad(0.0);
	image.fmode = FORMULATION_CLEN;
	image.x_centre = 255.5;
	image.y_centre = 255.5;
	image.camera_len = 0.2;  /* 20 cm */
	image.resolution = 5120; /* 512 pixels in 10 cm */
	image.lambda = 0.2e-9;   /* LCLS wavelength */
	image.qvecs = NULL;
	image.sfacs = NULL;
	image.data = NULL;
	
	get_diffraction(&image, cell);

	/* Write the output file */
	hdf5_write("results/sim.h5", image.sfacs, image.width, image.height);

	return 0;
}
