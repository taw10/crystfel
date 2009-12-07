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
#include "intensities.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Simulate diffraction patterns from small crystals probed with femtosecond\n"
"pulses of X-rays from a free electron laser.\n"
"\n"
" -h, --help                Display this help message.\n"
"     --simulation-details  Show technical details of the simulation.\n"
"\n"
"     --near-bragg          Output h,k,l,I near Bragg conditions.\n"
" -n, --number=<N>          Generate N images.  Default 1.\n"
"     --no-images           Do not output any HDF5 files.\n"
" -r, --random-orientation  Use a randomly generated orientation\n"
"                            (a new orientation will be used for each image).\n"
"\n"
"By default, the simulation aims to be as accurate as possible.  For greater\n"
"speed, or for testing, you can choose to disable certain things using the\n"
"following options.\n"
"\n"
"     --no-water            Do not simulate water background.\n"
"     --no-noise            Do not calculate Poisson noise.\n"
"     --no-bloom            Do not calculate CCD bloom (intensities which are\n"
"                            above the recordable range will be clamped to\n"
"                            the maximum allowable value).\n"
);
}


static void show_details()
{
	printf(
"This program simulates diffraction patterns from small crystals illuminated\n"
"with femtosecond X-ray pulses from a free electron laser.\n"
"\n"
"The lattice transform from the specified number of unit cells is calculated\n"
"using the closed-form solution for a truncated lattice:\n"
"\n"
"F(q) =  sin(pi*na*q.a)/sin(pi*q.a)\n"
"      * sin(pi*nb*q.b)/sin(pi*q.b)\n"
"      * sin(pi*nc*q.c)/sin(pi*q.c)\n"
"\n"
"na = number of unit cells in 'a' direction (likewise nb, nc)\n"
" q = reciprocal vector (1/d convention, not 2pi/d)\n"
"\n"
"This value is multiplied by the complex structure factor at the nearest\n"
"Bragg position, i.e. the gradient of the shape transform across each\n"
"appearance of the shape transform is not included, for speed of calculation.\n"
"\n"
"Complex structure factors are calculated using a combination of the Henke\n"
"and Waasmeier-Kirfel scattering factors. The Henke factors are complex\n"
"and energy dependence, whereas the Waas-Kirf values are real-valued and\n"
"|q|-dependent.  The difference between the Waas-Kirf value at the\n"
"appropriate |q| and the same value at |q|=0 is subtracted from the Henke\n"
"value.  The Henke values are linearly interpolated from the provided tables\n"
"(note that the interpolation should really be exponential).\n"
"\n"
"The modulus of the structure factor is taken and squared.  Intensity from\n"
"water is then added according to the first term of equation 5 from\n"
"Phys Chem Chem Phys 2003 (5) 1981--1991.\n"
"\n"
"Expected intensities at the CCD are then calculated using:\n"
"\n"
"I(q) = I0 * r^2 * |F(q)|^2 * S\n"
"\n"
"I0 = number of photons per unit area in the incident beam\n"
" r = Thomson radius\n"
" S = solid angle of corresponding pixel\n"
"\n"
"Poisson counts are generated from the expected intensities using Knuth's\n"
"algorithm.\n"
"\n"
"Bloom of the CCD is included.  Any excess intensity in a particular pixel\n"
"is divided between the neighbouring pixels.  Diagonal neighbours receive\n"
"half the contribution of adjacent pixels.  This process is repeated for\n"
"every pixel until all pixels are below the saturation value.  Note that this\n"
"process is slow for very saturated images.\n");
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
			ERROR("Invalid rotation '%s'\n", line);
		}

	} while ( 1 );
}


int main(int argc, char *argv[])
{
	int c;
	struct image image;
	int config_simdetails = 0;
	int config_nearbragg = 0;
	int config_randomquat = 0;
	int config_noimages = 0;
	int config_nowater = 0;
	int config_nonoise = 0;
	int config_nobloom = 0;
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
		{"no-images",          0, &config_noimages,    1},
		{"no-water",           0, &config_nowater,     1},
		{"no-noise",           0, &config_nonoise,     1},
		{"no-bloom",           0, &config_nobloom,     1},
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
	STATUS("Wavelength is %f nm\n", image.lambda/1.0e-9);

	do {

		int na, nb, nc;

		na = 8*random()/RAND_MAX + 4;
		nb = 8*random()/RAND_MAX + 4;
		nc = 16*random()/RAND_MAX + 30;
		STATUS("Particle size = %i x %i x %i\n", na, nb, nc);

		/* Read quaternion from stdin */
		if ( config_randomquat ) {
			image.orientation = random_quaternion();
		} else {
			image.orientation = read_quaternion();
		}

		STATUS("Orientation is %5.3f %5.3f %5.3f %5.3f\n",
		       image.orientation.w, image.orientation.x,
		       image.orientation.y, image.orientation.z);

		if ( !quaternion_valid(image.orientation) ) {
			ERROR("Orientation modulus is not zero!\n");
			return 1;
		}

		/* Ensure no residual information */
		image.qvecs = NULL;
		image.sfacs = NULL;
		image.data = NULL;
		image.twotheta = NULL;
		image.hdr = NULL;

		get_diffraction(&image, na, nb, nc);
		record_image(&image, !config_nowater, !config_nonoise,
		             !config_nobloom);

		if ( config_nearbragg ) {
			output_intensities(&image);
		}

		if ( !config_noimages ) {

			char filename[1024];

			snprintf(filename, 1023, "results/sim-%i.h5", number);
			number++;

			/* Write the output file */
			hdf5_write(filename, image.data,
			           image.width, image.height);

		}

		/* Clean up */
		free(image.data);
		free(image.qvecs);
		free(image.hdr);
		free(image.sfacs);
		free(image.twotheta);

		if ( n_images && (number >= n_images) ) done = 1;

	} while ( !done );

	return 0;
}
