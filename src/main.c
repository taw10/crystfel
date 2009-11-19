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
#include "detector.h"


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
	int c, done;
	UnitCell *cell;
	struct image image;
	char filename[1024];
	int number = 1;

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

again:

	/* Read quaternion from stdin */
	done = 0;
	do {

		int r;
		float w, x, y, z;
		char line[1024];
		char *rval;

		printf("Please input quaternion: w x y z\n");
		rval = fgets(line, 1023, stdin);
		if ( rval == NULL ) return 0;
		chomp(line);

		r = sscanf(line, "%f %f %f %f", &w, &x, &y, &z);
		if ( r == 4 ) {

			printf("Rotation is: %f %f %f %f (modulus=%f)\n",
			        w, x, y, z, sqrtf(w*w + x*x + y*y + z*z));

			image.orientation.w = w;
			image.orientation.x = x;
			image.orientation.y = y;
			image.orientation.z = z;

			done = 1;

		} else {
			fprintf(stderr, "Invalid rotation '%s'\n", line);
		}

	} while ( !done );

	/* Define image parameters */
	image.width = 512;
	image.height = 512;
	image.fmode = FORMULATION_CLEN;
	image.x_centre = 255.5;
	image.y_centre = 255.5;
	image.camera_len = 0.2;  /* 20 cm */
	image.resolution = 5120; /* 512 pixels in 10 cm */
	image.xray_energy = eV_to_J(2.0e3); /* 2 keV energy */
	image.lambda = ph_en_to_lambda(image.xray_energy);  /* Wavelength */
	image.qvecs = NULL;
	image.sfacs = NULL;
	image.data = NULL;
	image.twotheta = NULL;

	/* Splurge a few useful numbers */
	printf("Wavelength is %f nm\n", image.lambda/1.0e-9);

	get_diffraction(&image, cell);
	record_image(&image);

	snprintf(filename, 1023, "results/sim-%i.h5", number);
	number++;

	/* Write the output file */
	hdf5_write(filename, image.data, image.width, image.height);
	printf("Mock write: %i\n", number-1);
	goto again;
}
