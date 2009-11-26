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


static void show_help(const char *s)
{
	printf("Syntax: %s\n\n", s);
	printf("Simulate diffraction patterns from small crystals\n");
	printf(" probed with femosecond pulses from a free electron laser.\n\n");
	printf(" -h, --help                Display this help message\n");
	printf(" --simulation-details      Show details of the simulation\n");
	printf(" --near-bragg              Output h,k,l,I near Bragg conditions\n");
	printf(" -r, --random-orientation  Use a randomly generated orientation\n");
	printf("                            (a new orientation will be used for each image)\n");
	printf(" -n, --number=<N>          Generate N images.  Default 1\n");
}


static void show_details()
{
	printf("This program simulates diffraction patterns from small crystals illuminated\n");
	printf("with femtosecond X-ray pulses from a free electron laser.\n\n");
	printf("Scattering Factors\n");
	printf("------------------\n");
	printf("Scattering factors\n");
}


static struct quaternion read_quaternion()
{
	do {

		int r;
		float w, x, y, z;
		char line[1024];
		char *rval;

		printf("Please input quaternion: w x y z\n");
		rval = fgets(line, 1023, stdin);
		if ( rval == NULL ) return invalid_quaternion();
		chomp(line);

		r = sscanf(line, "%f %f %f %f", &w, &x, &y, &z);
		if ( r == 4 ) {

			struct quaternion quat;

			quat.w = w;
			quat.x = x;
			quat.y = y;
			quat.z = z;

			return quat;

		} else {
			fprintf(stderr, "Invalid rotation '%s'\n", line);
		}

	} while ( 1 );
}


int main(int argc, char *argv[])
{
	int c;
	struct image image;
	char filename[1024];
	int config_simdetails = 0;
	int config_nearbragg = 0;
	int config_randomquat = 0;
	int number = 1;  /* Index for the current image */
	int n_images = 1;  /* Generate one image by default */
	int done = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"simulation-details", 0, &config_simdetails,  1},
		{"near-bragg",         0, &config_nearbragg,   1},
		{"random-orientation", 0, NULL,               'r'},
		{"number",             1, NULL,               'n'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hrn:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'r' : {
			config_randomquat = 1;
			break;
		}

		case 'n' : {
			n_images = atoi(optarg);
			break;
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
		show_details();
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

	do {

		/* Read quaternion from stdin */
		if ( config_randomquat ) {
			image.orientation = random_quaternion();
		} else {
			image.orientation = read_quaternion();
		}

		if ( quaternion_valid(image.orientation) ) {
			fprintf(stderr, "Orientation modulus is not zero!\n");
			return 1;
		}

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

		if ( n_images && (number >= n_images) ) done = 1;

	} while ( !done );
}
