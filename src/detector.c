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
#include "diffraction.h"
#include "detector.h"


/* Number of photons in pulse */
#define FLUENCE (1.0e13)

/* Detector's quantum efficiency */
#define DQE (0.9)

/* Detector's saturation value */
#define SATURATION (5000)

/* Radius of the water column */
#define WATER_RADIUS (3.0e-6 / 2.0)

/* Radius of X-ray beam */
#define BEAM_RADIUS (3.0e-6 / 2.0)


/* Bleed excess intensity into neighbouring pixels */
static void bloom_values(int *tmp, int x, int y,
                          int width, int height, int val)
{
	int overflow;

	overflow = val - SATURATION;

	/* Intensity which bleeds off the edge of the detector is lost */
	if ( x > 0 ) {
		tmp[x-1 + width*y] += overflow / 6;
		if ( y > 0 ) {
			tmp[x-1 + width*(y-1)] += overflow / 12;
		}
		if ( y < height-1 ) {
			tmp[x-1 + width*(y+1)] += overflow / 12;
		}
	}

	if ( x < width-1 ) {
		tmp[x+1 + width*y] += overflow / 6;
		if ( y > 0 ) {
			tmp[x+1 + width*(y-1)] += overflow / 12;
		}
		if ( y < height-1 ) {
			tmp[x+1 + width*(y+1)] += overflow / 12;
		}
	}

	if ( y > 0 ) {
		tmp[x + width*(y-1)] += overflow / 6;
	}
	if ( y < height-1 ) {
		tmp[x + width*(y+1)] += overflow / 6;
	}
}


static uint16_t *bloom(int *hdr_in, int width, int height)
{
	int x, y;
	uint16_t *data;
	int *tmp;
	int *hdr;
	int did_something;

	data = malloc(width * height * sizeof(uint16_t));
	tmp = malloc(width * height * sizeof(int));
	hdr = malloc(width * height * sizeof(int));

	memcpy(hdr, hdr_in, width*height*sizeof(int));

	/* Apply DQE (once only) */
	for ( x=0; x<width; x++ ) {
	for ( y=0; y<height; y++ ) {
		 hdr[x + width*y] *= DQE;
	}
	}

	do {

		memset(tmp, 0, width*height*sizeof(int));
		did_something = 0;

		for ( x=0; x<width; x++ ) {
		for ( y=0; y<height; y++ ) {

			int hdval;

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
			memcpy(hdr, tmp, width*height*sizeof(int));
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


void record_image(struct image *image, int do_water, int do_poisson,
                  int do_bloom)
{
	int x, y;
	double total_energy, energy_density;
	double ph_per_e;
	double pix_area, Lsq;
	double area;

	/* How many photons are scattered per electron? */
	area = M_PI*pow(BEAM_RADIUS, 2.0);
	total_energy = FLUENCE * image->xray_energy;
	energy_density = total_energy / area;
	ph_per_e = (FLUENCE/area) * pow(THOMSON_LENGTH, 2.0);
	STATUS("Fluence = %8.2e photons, "
	       "Energy density = %5.3f kJ/cm^2, "
	       "Total energy = %5.3f microJ\n",
	       FLUENCE, energy_density/1e7, total_energy*1e6);

	/* Area of one pixel */
	pix_area = pow(1.0/image->resolution, 2.0);
	Lsq = pow(image->camera_len, 2.0);

	image->hdr = malloc(image->width * image->height * sizeof(double));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		int counts;
		double cf;
		double intensity, sa, water;
		double complex val;
		double dsq, proj_area;

		val = image->sfacs[x + image->width*y];
		intensity = pow(cabs(val), 2.0);

		if ( do_water ) {

			/* Add intensity contribution from water */
			water = water_intensity(image->qvecs[x + image->width*y],
			                        image->xray_energy,
			                        BEAM_RADIUS, WATER_RADIUS);
			intensity += water;

		}

		/* Area of pixel as seen from crystal (approximate) */
		proj_area = pix_area * cos(image->twotheta[x + image->width*y]);

		/* Calculate distance from crystal to pixel */
		dsq = pow(((double)x - image->x_centre)/image->resolution, 2.0);
		dsq += pow(((double)y - image->y_centre)/image->resolution, 2.0);

		/* Projected area of pixel divided by distance squared */
		sa = proj_area / (dsq + Lsq);

		if ( do_poisson ) {
			counts = poisson_noise(intensity * ph_per_e * sa * DQE);
		} else {
			double rounded;
			cf = intensity * ph_per_e * sa * DQE;
			rounded = rint(cf);
			counts = (int)rounded;
		}

		image->hdr[x + image->width*y] = counts;

	}
	progress_bar(x, image->width-1, "Post-processing");
	}

	if ( do_bloom ) {
		image->data = bloom(image->hdr, image->width, image->height);
	} else {
		image->data = malloc(image->width * image->height
		                                            * sizeof(uint16_t));
		for ( x=0; x<image->width; x++ ) {
		for ( y=0; y<image->height; y++ ) {
			int val;
			val = image->hdr[x + image->width*y];
			if ( val > SATURATION ) val = SATURATION;
	                image->data[x + image->width*y] = (uint16_t)val;
		}
		}
	}
}
