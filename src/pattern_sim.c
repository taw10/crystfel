/*
 * pattern_sim.c
 *
 * Simulate diffraction patterns from small crystals
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

#include "image.h"
#include "diffraction.h"
#include "diffraction-gpu.h"
#include "cell.h"
#include "utils.h"
#include "hdf5-file.h"
#include "detector.h"
#include "intensities.h"
#include "sfac.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Simulate diffraction patterns from small crystals probed with femtosecond\n"
"pulses of X-rays from a free electron laser.\n"
"\n"
" -h, --help                Display this help message.\n"
"     --simulation-details  Show technical details of the simulation.\n"
"     --gpu                 Use the GPU to speed up the calculation.\n"
"\n"
"     --near-bragg          Output h,k,l,I near Bragg conditions.\n"
" -n, --number=<N>          Generate N images.  Default 1.\n"
"     --no-images           Do not output any HDF5 files.\n"
" -r, --random-orientation  Use a randomly generated orientation\n"
"                            (a new orientation will be used for each image).\n"
"     --powder              Write a summed pattern of all images simulated by\n"
"                            this invocation to results/integr.h5.\n"
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
"     --no-sfac             Pretend that all structure factors are 1.\n"
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
"algorithm.  When the intensity is sufficiently high that Knuth's algorithm\n"
"would result in machine precision problems, a normal distribution with\n"
"standard deviation sqrt(I) is used instead.\n"
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
	struct gpu_context *gctx = NULL;
	long long int *powder;
	int config_simdetails = 0;
	int config_nearbragg = 0;
	int config_randomquat = 0;
	int config_noimages = 0;
	int config_nowater = 0;
	int config_nonoise = 0;
	int config_nobloom = 0;
	int config_nosfac = 0;
	int config_gpu = 0;
	int config_powder = 0;
	int ndone = 0;    /* Number of simulations done (images or not) */
	int number = 1;   /* Number used for filename of image */
	int n_images = 1; /* Generate one image by default */
	int done = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"simulation-details", 0, &config_simdetails,  1},
		{"gpu",                0, &config_gpu,         1},
		{"near-bragg",         0, &config_nearbragg,   1},
		{"random-orientation", 0, NULL,               'r'},
		{"number",             1, NULL,               'n'},
		{"no-images",          0, &config_noimages,    1},
		{"no-water",           0, &config_nowater,     1},
		{"no-noise",           0, &config_nonoise,     1},
		{"no-bloom",           0, &config_nobloom,     1},
		{"no-sfac",            0, &config_nosfac,      1},
		{"powder",             0, &config_powder,      1},
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
	image.lambda = ph_en_to_lambda(eV_to_J(2.0e3));  /* Wavelength */
	image.molecule = load_molecule();

	#include "geometry-lcls.tmp"

	powder = calloc(image.width*image.height, sizeof(*powder));

	/* Splurge a few useful numbers */
	STATUS("Wavelength is %f nm\n", image.lambda/1.0e-9);

	do {

		int na, nb, nc;
		double a, b, c, d;

		//na = 8*random()/RAND_MAX + 4;
		//nb = 8*random()/RAND_MAX + 4;
		//nc = 16*random()/RAND_MAX + 30;
		na = 8;
		nb = 8;
		nc = 16;

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
		image.sfacs = NULL;
		image.data = NULL;
		image.twotheta = NULL;
		image.hdr = NULL;

		cell_get_parameters(image.molecule->cell, &a, &b, &c, &d, &d, &d);
		STATUS("Particle size = %i x %i x %i (=%5.2f x %5.2f x %5.2f nm)\n",
	               na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

		if ( config_gpu ) {
			if ( gctx == NULL ) {
				gctx = setup_gpu(config_nosfac, &image,
				                 image.molecule);
			}
			get_diffraction_gpu(gctx, &image, na, nb, nc);
		} else {
			get_diffraction(&image, na, nb, nc, config_nosfac);
		}
		if ( image.molecule == NULL ) {
			ERROR("Couldn't open molecule.pdb\n");
			return 1;
		}
		if ( image.sfacs == NULL ) {
			ERROR("Diffraction calculation failed.\n");
			goto skip;
		}

		record_image(&image, !config_nowater, !config_nonoise,
		             !config_nobloom);

		if ( config_nearbragg ) {
			output_intensities(&image, image.molecule->cell);
		}

		if ( config_powder ) {

			int x, y, w;

			w = image.width;

			for ( x=0; x<image.width; x++ ) {
			for ( y=0; y<image.height; y++ ) {
				powder[x+w*y] += image.data[x+w*y];
			}
			}

			if ( !(ndone % 10) ) {
				hdf5_write("results/integr.h5", powder,
				           image.width, image.height,
				           H5T_NATIVE_LLONG);
			}
		}

		if ( !config_noimages ) {

			char filename[1024];

			snprintf(filename, 1023, "results/sim-%i.h5", number);
			number++;

			/* Write the output file */
			hdf5_write(filename, image.data,
			           image.width, image.height, H5T_NATIVE_INT16);

		}

		/* Clean up */
		free(image.data);
		free(image.hdr);
		free(image.sfacs);
		free(image.twotheta);

skip:
		ndone++;

		if ( n_images && (ndone >= n_images) ) done = 1;

	} while ( !done );

	if ( gctx != NULL ) {
		cleanup_gpu(gctx);
	}

	free(image.det.panels);
	free(powder);
	free(image.molecule->reflections);
	free(image.molecule);

	return 0;
}
