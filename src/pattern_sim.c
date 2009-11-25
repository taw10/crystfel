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
#include <getopt.h>

#include "image.h"
#include "diffraction.h"
#include "cell.h"
#include "utils.h"
#include "hdf5-file.h"
#include "detector.h"


/* Crystal size in metres */
#define CRYSTAL_SIZE (500.0e-9)


static void show_help(const char *s)
{
	printf("Syntax: %s\n\n", s);
	printf("Simulate diffraction patterns from small crystals\n");
	printf(" probed with femosecond pulses from a free electron laser\n\n");
	printf("  -h, --help            Display this help message\n");
	printf("  --simulation-details  Show details of the simulation\n");
}


static void show_details(const char *s)
{
	printf("%s: Simulation Details\n\n", s);
	printf("Simulates diffraction patterns from small crystals\n");
	printf("probed with femtosecond pulses from a free electron laser\n\n");
}


int main(int argc, char *argv[])
{
	int c, done;
	struct image image;
	char filename[1024];
	int number = 1;
	int config_simdetails = 0;

	const struct option longopts[] = {
	      {"help", 0, NULL, 'h'},
	      {"simulation-details", 0, &config_simdetails, 1},
	      {0, 0, NULL, 0}
	};

	while ((c = getopt_long(argc, argv, "h", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 0 : {
			break;
		}

		default : {
			return 1;
		}
		}

	}

	if ( config_simdetails ) {
		show_details(argv[0]);
		return 0;
	}

	/* Define image parameters */
	image.width = 1024;
	image.height = 1024;
	image.fmode = FORMULATION_CLEN;
	image.x_centre = 512.5;
	image.y_centre = 512.5;
	image.camera_len = 0.05;  /* 5 cm (front CCD can move from 5cm-20cm) */
	image.resolution = 13333.3; /* 75 micron pixel size */
	image.xray_energy = eV_to_J(2.0e3); /* 2 keV energy */
	image.lambda = ph_en_to_lambda(image.xray_energy);  /* Wavelength */
	image.molecule = NULL;

	/* Splurge a few useful numbers */
	printf("Wavelength is %f nm\n", image.lambda/1.0e-9);

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

	/* Ensure no residual information */
	image.qvecs = NULL;
	image.sfacs = NULL;
	image.data = NULL;
	image.twotheta = NULL;
	image.hdr = NULL;

	get_diffraction(&image);
	record_image(&image);

	snprintf(filename, 1023, "results/sim-%i.h5", number);
	number++;

	/* Write the output file */
	hdf5_write(filename, image.data, image.width, image.height);

	/* Clean up */
	free(image.data);
	free(image.qvecs);
	free(image.hdr);
	free(image.sfacs);
	free(image.twotheta);

	/* Do it all again */
	goto again;
}
