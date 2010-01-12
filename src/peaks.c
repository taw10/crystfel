/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "image.h"
#include "utils.h"


#define PEAK_WINDOW_SIZE (10)

void search_peaks(struct image *image, int dump_peaks)
{
	int x, y, width, height;
	int16_t *data;

	data = image->data;
	width = image->width;
	height = image->height;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	for ( x=1; x<image->width-1; x++ ) {
	for ( y=1; y<image->height-1; y++ ) {

		double dx1, dx2, dy1, dy2;
		double dxs, dys;
		double grad;
		int mask_x, mask_y;
		int sx, sy;
		double max;
		unsigned int did_something = 1;

		/* Overall threshold */
		if ( data[x+width*y] < 800 ) continue;

		/* Ignore streak */
		if ( abs(x-image->x_centre) < 15 ) continue;

		/* Get gradients */
		dx1 = data[x+width*y] - data[(x+1)+width*y];
		dx2 = data[(x-1)+width*y] - data[x+width*y];
		dy1 = data[x+width*y] - data[(x+1)+width*(y+1)];
		dy2 = data[x+width*(y-1)] - data[x+width*y];

		/* Average gradient measurements from both sides */
		dxs = ((dx1*dx1) + (dx2*dx2)) / 2;
		dys = ((dy1*dy1) + (dy2*dy2)) / 2;

		/* Calculate overall gradient */
		grad = dxs + dys;

		if ( grad < 2000000 ) continue;

		mask_x = x;
		mask_y = y;

		while ( (did_something)
		     && (distance(mask_x, mask_y, x, y)<50) ) {

			max = data[mask_x+width*mask_y];
			did_something = 0;

			for ( sy=biggest(mask_y-PEAK_WINDOW_SIZE/2, 0);
			      sy<smallest(mask_y+PEAK_WINDOW_SIZE/2, height);
			      sy++ ) {
			for ( sx=biggest(mask_x-PEAK_WINDOW_SIZE/2, 0);
			      sx<smallest(mask_x+PEAK_WINDOW_SIZE/2, width);
			      sx++ ) {

				if ( data[sx+width*sy] > max ) {
					max = data[sx+width*sy];
					mask_x = sx;
					mask_y = sy;
					did_something = 1;
				}

			}
			}

		}

		if ( !did_something ) {

			double d;
			int idx;

			assert(mask_x<image->width);
			assert(mask_y<image->height);
			assert(mask_x>=0);
			assert(mask_y>=0);

			/* Too far from foot point? */
			if ( distance(mask_x, mask_y, x, y) > 50.0 )  continue;

			/* Check for a feature at exactly the
			 * same coordinates */
			image_feature_closest(image->features, mask_x, mask_y,
			                      &d, &idx);

			if ( d > 1.0 ) {

				/* Map and record reflection */
				if ( dump_peaks ) {
					printf("%i %i\n", x, y);
				}

				image_add_feature(image->features,
				                  mask_x, mask_y, image, 1.0);

			}

		}


	}
	}
}
