/*
 * filters.c
 *
 * Image filtering
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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
#include <gsl/gsl_statistics_int.h>

#include "image.h"


static int compare_vals(const void *ap, const void *bp)
{
	const signed int a = *(signed int *)ap;
	const signed int b = *(signed int *)bp;

	if ( a > b ) return 1;
	if ( a < b ) return -1;
	return 0;
}


static void clean_panel(struct image *image, int sx, int sy)
{
	int x, y;
	const int s = sizeof(signed int);

	for ( x=0; x<512; x++ ) {

		signed int vals[128];
		double m;

		for ( y=0; y<128; y++ ) {
			vals[y] = image->data[(x+sx)+(y+sy)*image->width];
		}

		qsort(&vals[0], 128, s, compare_vals);

		m = gsl_stats_int_median_from_sorted_data(vals, 1, 128);

		for ( y=0; y<128; y++ ) {
			image->data[(x+sx)+(y+sy)*image->width] -= m;
		}

	}
}


static void noise_filter(struct image *image)
{
	int x, y;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		int dx, dy;
		int val = image->data[x+image->width*y];

		if ( (x==0) || (x==image->width-1)
		  || (y==0) || (y==image->height-1) ) {
			if ( val < 0 ) val = 0;
			continue;
		}

		for ( dx=-1; dx<=+1; dx++ ) {
		for ( dy=-1; dy<=+1; dy++ ) {

			int val2;

			val2 = image->data[(x+dx)+image->width*(y+dy)];

			if ( val2 < 0 ) val = 0;

		}
		}

		image->data[x+image->width*y] = val;

	}
	}
}


/* Pre-processing to make life easier */
void clean_image(struct image *image)
{
	int px, py;

	if ( (image->width != 1024) || (image->height != 1024) ) return;

	for ( px=0; px<2; px++ ) {
	for ( py=0; py<8; py++ ) {

		clean_panel(image, 512*px, 128*py);

	}
	}

	noise_filter(image);
}
