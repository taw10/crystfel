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
#include "parameters-lcls.tmp"


/* x,y in pixels relative to central beam */
int map_position(struct image *image, double dx, double dy,
                 double *rx, double *ry, double *rz)
{
	double d;
	double twotheta, psi;
	const double k = 1.0 / image->lambda;
	struct panel *p;
	double x = 0.0;
	double y = 0.0;

	p = find_panel(&image->det, dx, dy);

	x = ((double)dx - p->cx);
	y = ((double)dy - p->cy);

	/* Convert pixels to metres */
	x /= p->res;
	y /= p->res;	/* Convert pixels to metres */
	d = sqrt((x*x) + (y*y));
	twotheta = atan2(d, p->clen);

	psi = atan2(y, x);

	*rx = k*sin(twotheta)*cos(psi);
	*ry = k*sin(twotheta)*sin(psi);
	*rz = k - k*cos(twotheta);

	return 0;
}


void record_image(struct image *image, int do_poisson)
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

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double counts;
		double cf;
		double intensity, sa;
		double pix_area, Lsq;
		double dsq, proj_area;
		struct panel *p;

		intensity = (double)image->data[x + image->width*y];
		if ( isinf(intensity) ) {
			ERROR("Infinity at %i,%i\n", x, y);
		}

		p = find_panel(&image->det, x, y);

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
			cf = intensity * ph_per_e * sa * DQE;
			counts = rint(cf);
			if ( counts < 0.0 ) {
				ERROR("Negative at %i,%i %f\n", x, y, counts);
			}
		}

		image->data[x + image->width*y] = counts * DETECTOR_GAIN;
		if ( isinf(image->data[x+image->width*y]) ) {
			ERROR("Processed infinity at %i,%i\n", x, y);
		}
		if ( image->data[x+image->width*y] < 0.0 ) {
			ERROR("Processed negative at %i,%i %f\n", x, y, counts);
		}

	}
	progress_bar(x, image->width-1, "Post-processing");
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
