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
	
	/* How many photons are scattered per electron? */
	ph_per_e = PULSE_ENERGY_DENSITY * pow(THOMSON_LENGTH, 2.0)
	          / image->xray_energy;
	
	printf("%e photons are scattered per electron\n", ph_per_e);
	
	image->data = malloc(image->width * image->height * sizeof(uint16_t));
	
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
	
		uint16_t counts;
		double val, intensity;
		double sa;
		
		val = image->sfacs[x + image->width*y];
		
		/* What solid angle is subtended by this pixel? */
		sa = 1.0;
		
		intensity = pow(val, 2.0);
		counts = intensity * ph_per_e * sa;
		
		image->data[x + image->width*y] = counts;
	
	}
	}
}
