/*
 * filters.c
 *
 * Image filtering
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
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
#include <gsl/gsl_blas.h>

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


/* Pre-processing to make life easier */
void filter_cm(struct image *image)
{
	int px, py;

	if ( (image->width != 1024) || (image->height != 1024) ) return;

	for ( px=0; px<2; px++ ) {
	for ( py=0; py<8; py++ ) {

		clean_panel(image, 512*px, 128*py);

	}
	}

}


void filter_noise(struct image *image, float *old)
{
	int x, y;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		int dx, dy;
		int val = image->data[x+image->width*y];

		if ( old != NULL ) old[x+image->width*y] = val;

		/* FIXME: This isn't really the right thing to do
		 * at the edges. */
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


/* Force the linker to bring in CBLAS to make GSL happy */
void filters_fudge_gslcblas()
{
        STATUS("%p\n", cblas_sgemm);
}
