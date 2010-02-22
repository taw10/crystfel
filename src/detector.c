/*
 * detector.c
 *
 * Detector properties
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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


static int16_t *bloom(int *hdr_in, int width, int height)
{
	int x, y;
	int16_t *data;
	int *tmp;
	int *hdr;
	int did_something;

	data = malloc(width * height * sizeof(int16_t));
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
		data[x + width*y] = (int16_t)tmp[x + width*y];
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
	double area;

	/* How many photons are scattered per electron? */
	area = M_PI*pow(BEAM_RADIUS, 2.0);
	total_energy = FLUENCE * ph_lambda_to_en(image->lambda);
	energy_density = total_energy / area;
	ph_per_e = (FLUENCE/area) * pow(THOMSON_LENGTH, 2.0);
	STATUS("Fluence = %8.2e photons, "
	       "Energy density = %5.3f kJ/cm^2, "
	       "Total energy = %5.3f microJ\n",
	       FLUENCE, energy_density/1e7, total_energy*1e6);

	image->hdr = malloc(image->width * image->height * sizeof(double));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		int counts;
		double cf;
		double intensity, sa, water;
		double complex val;
		double pix_area, Lsq;
		double dsq, proj_area;
		int p;
		int found = 0;

		val = image->sfacs[x + image->width*y];
		intensity = pow(cabs(val), 2.0);

		for ( p=0; p<image->det.n_panels; p++ ) {
			if ( (x >= image->det.panels[p].min_x)
			  && (x <= image->det.panels[p].max_x)
			  && (y >= image->det.panels[p].min_y)
			  && (y <= image->det.panels[p].max_y) ) {
				found = 1;
				break;
			}
		}
		if ( !found ) {
			ERROR("No mapping found for %i,%i\n", x, y);
			return;
		}

		/* FIXME: Move to diffraction.c somehow */
		if ( do_water ) {

			struct rvec q;

			q = get_q(image, x, y, 1, NULL, 1.0/image->lambda);

			/* Add intensity contribution from water */
			water = water_intensity(q,
			                        ph_lambda_to_en(image->lambda),
			                        BEAM_RADIUS, WATER_RADIUS);
			intensity += water;

		}

		/* Area of one pixel */
		pix_area = pow(1.0/image->det.panels[p].res, 2.0);
		Lsq = pow(image->det.panels[p].clen, 2.0);

		/* Area of pixel as seen from crystal (approximate) */
		proj_area = pix_area * cos(image->twotheta[x + image->width*y]);

		/* Calculate distance from crystal to pixel */
		dsq = pow(((double)x - image->det.panels[p].cx)
		                               / image->det.panels[p].res, 2.0);
		dsq += pow(((double)y - image->det.panels[p].cy)
		                               / image->det.panels[p].res, 2.0);

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
		                                            * sizeof(int16_t));
		for ( x=0; x<image->width; x++ ) {
		for ( y=0; y<image->height; y++ ) {
			int val;
			val = image->hdr[x + image->width*y];
			if ( val > SATURATION ) val = SATURATION;
	                image->data[x + image->width*y] = (int16_t)val;
		}
		}
	}
}
