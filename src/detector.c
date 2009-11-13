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


void record_image(struct image *image)
{
	int x, y;
	
	image->data = malloc(image->width * image->height * sizeof(uint16_t));
	
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
	
		uint16_t counts;
		double val;
		double intensity;
		double phac;
		
		val = image->sfacs[x + image->width*y];
		phac = image->phactors[x + image->width*y];
		
		intensity = pow(val, 2.0);
		counts = intensity * phac;
		
		image->data[x + image->width*y] = counts;
	
	}
	}
}
