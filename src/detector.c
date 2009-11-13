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

#include "image.h"
#include "utils.h"


/* Pulse energy density in J/m^2 */
#define PULSE_ENERGY_DENSITY (30.0e7)


void record_image(struct image *image)
{
	int x, y;
	double ph_per_e;
	double twotheta_max, np, sa_per_pixel;
	
	/* How many photons are scattered per electron? */
	ph_per_e = PULSE_ENERGY_DENSITY * pow(THOMSON_LENGTH, 2.0)
	          / image->xray_energy;
	printf("%e photons are scattered per electron\n", ph_per_e);
	
	twotheta_max = image->twotheta[0];
	
	np = sqrt(pow(image->x_centre, 2.0) + pow(image->y_centre, 2.0));
	sa_per_pixel = pow(2.0 * twotheta_max / np, 2.0);
	printf("sa per pixel=%e\n", sa_per_pixel);
	
	image->data = malloc(image->width * image->height * sizeof(uint16_t));
	
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
	
		double counts;
		double val, intensity;
		double sa;
		
		val = image->sfacs[x + image->width*y];
		
		/* What solid angle is subtended by this pixel? */
		sa = sa_per_pixel * cos(image->twotheta[x + image->width*y]);
		
		intensity = pow(val, 2.0);
		counts = intensity * ph_per_e * sa;
		
		image->data[x + image->width*y] = counts;
	
	}
	}
}
