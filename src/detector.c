/*
 * detector.c
 *
 * Detector properties
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "image.h"
#include "utils.h"


/* Pulse energy density in J/m^2 */
#define PULSE_ENERGY_DENSITY (30.0e7)

/* Detector's quantum efficiency */
#define DQE (0.9)

/* Detector's saturation value */
#define SATURATION (60000)


/* Bleed excess intensity into neighbouring pixels */
static void bloom_values(double *tmp, int x, int y,
                          int width, int height, double val)
{
	double overflow;

	overflow = val - SATURATION;

	/* Intensity which bleeds off the edge of the detector is lost */
	if ( x > 0 ) {
		tmp[x-1 + width*y] += overflow / 6.0;
		if ( y > 0 ) {
			tmp[x-1 + width*(y-1)] += overflow / 12.0;
		}
		if ( y < height-1 ) {
			tmp[x-1 + width*(y+1)] += overflow / 12.0;
		}
	}

	if ( x < width-1 ) {
		tmp[x+1 + width*y] += overflow / 6.0;
		if ( y > 0 ) {
			tmp[x+1 + width*(y-1)] += overflow / 12.0;
		}
		if ( y < height-1 ) {
			tmp[x+1 + width*(y+1)] += overflow / 12.0;
		}
	}

	if ( y > 0 ) {
		tmp[x + width*(y-1)] += overflow / 6.0;
	}
	if ( y < height-1 ) {
		tmp[x + width*(y+1)] += overflow / 6.0;
	}
}


static uint16_t *bloom(double *hdr_in, int width, int height)
{
	int x, y;
	uint16_t *data;
	double *tmp;
	double *hdr;
	int did_something;

	data = malloc(width * height * sizeof(uint16_t));
	tmp = malloc(width * height * sizeof(double));
	hdr = malloc(width * height * sizeof(double));

	memcpy(hdr, hdr_in, width*height*sizeof(double));

	/* Apply DQE (once only) */
	for ( x=0; x<width; x++ ) {
	for ( y=0; y<height; y++ ) {
		 hdr[x + width*y] *= DQE;
	}
	}

	do {

		memset(tmp, 0, width*height*sizeof(double));
		did_something = 0;

		for ( x=0; x<width; x++ ) {
		for ( y=0; y<height; y++ ) {

			double hdval;

			hdval = hdr[x + width*y];

			/* Pixel not saturated? */
			if ( hdval <= SATURATION ) {
				tmp[x + width*y] += hdval;
				continue;
			}

			bloom_values(tmp, x, y, width, height, hdval);
			tmp[x + width*y] += SATURATION;
			did_something = 1;

		}
		}

		/* Prepare new input if we're going round for another pass */
		if ( did_something ) {
			memcpy(hdr, tmp, width*height*sizeof(double));
		}

	} while ( did_something );

	/* Turn into integer array of counts */
	for ( x=0; x<width; x++ ) {
	for ( y=0; y<height; y++ ) {
		data[x + width*y] = (uint16_t)tmp[x + width*y];
	}
	}

	free(tmp);
	free(hdr);

	return data;
}


void record_image(struct image *image)
{
	int x, y;
	double ph_per_e;
	double sa_per_pixel;
	const int do_bloom = 0;

	/* How many photons are scattered per electron? */
	ph_per_e = PULSE_ENERGY_DENSITY * pow(THOMSON_LENGTH, 2.0)
	          / image->xray_energy;
	printf("%e photons are scattered per electron\n", ph_per_e);

	/* Solid angle subtended by a single pixel at the centre of the CCD */
	sa_per_pixel = pow(1.0/image->resolution, 2.0) / image->camera_len;
	printf("Solid angle of one pixel (at centre of CCD) = %e sr\n",
	       sa_per_pixel);

	image->hdr = malloc(image->width * image->height * sizeof(double));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double counts;
		double intensity;
		double sa;
		double complex val;

		val = image->sfacs[x + image->width*y];
		intensity = (double)(val * conj(val));

		/* What solid angle is subtended by this pixel? */
		sa = sa_per_pixel * cos(image->twotheta[x + image->width*y]);

		counts = intensity * ph_per_e * sa;

		image->hdr[x + image->width*y] = counts;

	}
	}

	if ( do_bloom ) {
		image->data = bloom(image->hdr, image->width, image->height);
	} else {
		image->data = malloc(image->width * image->height
		                                            * sizeof(uint16_t));
		for ( x=0; x<image->width; x++ ) {
		for ( y=0; y<image->height; y++ ) {
			double val;
			val = image->hdr[x + image->width*y];
	                image->data[x + image->width*y] = (uint16_t)val;
		}
		}
	}
}
