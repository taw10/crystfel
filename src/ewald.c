/*
 * ewald.c
 *
 * Calculate q-vector arrays
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
#include "cell.h"


void get_ewald(struct image *image)
{
	int x, y;
	double k;  /* Wavenumber */
	double i0fac;
	
	k = 1/image->lambda;
	
	image->qvecs = malloc(image->width * image->height
	                       * sizeof(struct threevec));
	
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {
	
		double rx, ry, r;
		double twothetax, twothetay, twotheta;
		double qx, qy, qz;
		
		/* Calculate q vectors for Ewald sphere */
		rx = ((double)x - image->x_centre) / image->resolution;
		ry = ((double)y - image->y_centre) / image->resolution;
		r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
		
		twothetax = atan2(rx, image->camera_len);
		twothetay = atan2(ry, image->camera_len);	
		twotheta = atan2(r, image->camera_len);
		
		qx = k * sin(twothetax);
		qy = k * sin(twothetay);
		qz = k - k * cos(twotheta);
		
		image->qvecs[x + image->width*y].u = qx;
		image->qvecs[x + image->width*y].v = qy;
		image->qvecs[x + image->width*y].w = qz;
	
	}
	}
}
