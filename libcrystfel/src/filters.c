/*
 * filters.c
 *
 * Image filtering
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
 *   2013      Anton Barty <anton.barty@desy.de>
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


void filter_noise(struct image *image)
{
	int x, y;

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		int dx, dy;
		float val = image->data[x+image->width*y];

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


#define SWAP(a,b) { float t=(a);(a)=(b);(b)=t; }
static float kth_smallest(float *a, int n, int k)
{
	long i, j, l, m;
	float x;

	l = 0;
	m = n-1;

	while ( l < m ) {
		x=a[k];
		i=l;
		j=m;
		do {
			while (a[i]<x) i++;
			while (x<a[j]) j--;
			if ( i<=j ) {
				SWAP(a[i],a[j]);
				i++;
				j--;
			}
		} while (i<=j);
		if ( j<k ) l = i;
		if ( k<i ) m = j;
	}
	return a[k];
}
#undef SWAP


void filter_median(struct image *image, int size)
{
	int counter;
	int nn;
	float *buffer;
	float *localBg;
	int pn;
	int i;

	if ( size <= 0 ) return;

	nn = 2*size+1;
	nn = nn*nn;

	buffer = calloc(nn, sizeof(float));
	localBg = calloc(image->width*image->height, sizeof(float));
	if ( (buffer == NULL) || (localBg == NULL) ) {
		ERROR("Failed to allocate LB buffers.\n");
		return;
	}

	/* Determine local background
	 * (median over window width either side of current pixel) */
	for ( pn=0; pn<image->det->n_panels; pn++ ) {

		int fs, ss;
		struct panel *p;
		int p_w, p_h;

		p = &image->det->panels[pn];
		p_w = (p->max_fs-p->min_fs)+1;
		p_h = (p->max_ss-p->min_ss)+1;

		for ( fs=0; fs<p_w; fs++ ) {
		for ( ss=0; ss<p_h; ss++ ) {

			int ifs, iss;
			int e;

			counter = 0;
			e = fs+p->min_fs;
			e += (ss+p->min_ss)*image->width;

			// Loop over median window
			for ( ifs=-size; ifs<=size; ifs++ ) {
			for ( iss=-size; iss<=size; iss++ ) {

				int idx;

				if ( (fs+ifs) < 0 ) continue;
				if ( (fs+ifs) >= p_w ) continue;
				if ( (ss+iss) < 0 ) continue;
				if ( (ss+iss) >= p_h ) continue;

				idx = fs+ifs+p->min_fs;
				idx += (ss+iss+p->min_ss)*image->width;

				buffer[counter++] = image->data[idx];

			}
			}

			// Find median value
			localBg[e] = kth_smallest(buffer, counter, counter/2);

		}
		}
	}


	/* Do the background subtraction */
	for ( i=0; i<image->width*image->height; i++ ) {
		image->data[i] -= localBg[i];
	}

	free(localBg);
	free(buffer);
}
