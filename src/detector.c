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

/* Detector's quantum efficiency (ADU per photon, front detector) */
#define DQE (167.0)

/* Radius of the water column */
#define WATER_RADIUS (3.0e-6 / 2.0)

/* Radius of X-ray beam */
#define BEAM_RADIUS (3.0e-6 / 2.0)


void record_image(struct image *image, int do_water, int do_poisson)
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
		struct panel *p;

		val = image->sfacs[x + image->width*y];
		intensity = pow(cabs(val), 2.0);

		p = find_panel(&image->det, x, y);

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
		pix_area = pow(1.0/p->res, 2.0);
		Lsq = pow(p->clen, 2.0);

		/* Area of pixel as seen from crystal (approximate) */
		proj_area = pix_area * cos(image->twotheta[x + image->width*y]);

		/* Calculate distance from crystal to pixel */
		dsq = pow(((double)x - p->cx) / p->res, 2.0);
		dsq += pow(((double)y - p->cy) / p->res, 2.0);

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

	image->data = malloc(image->width * image->height * sizeof(float));
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
		int val;
		val = image->hdr[x + image->width*y];
		image->data[x + image->width*y] = val;
	}
	}
}


struct panel *find_panel(struct detector *det, int x, int y)
{
	int p;

	for ( p=0; p<det->n_panels; p++ ) {
		if ( (x >= det->panels[p].min_x)
		  && (x <= det->panels[p].max_x)
		  && (y >= det->panels[p].min_y)
		  && (y <= det->panels[p].max_y) ) {
			return &det->panels[p];
		}
	}
	ERROR("No mapping found for %i,%i\n", x, y);

	return NULL;
}
