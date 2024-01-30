/*
 * filters.c
 *
 * Image filtering
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_blas.h>

#include "image.h"

/** \file filters.h */

static void filter_noise_in_panel(float *data, int width, int height)
{
	int x, y;

	for ( y=0; y<height; y++ ) {
	for ( x=0; x<width; x++ ) {

		int dx, dy;
		float val = data[x+width*y];

		/* This isn't really the right thing to do at the edges,
		 * but this filter is so horrible it's unlikely to matter. */
		if ( (x==0) || (x==width-1)
		  || (y==0) || (y==height-1) ) {
			if ( val < 0 ) {
				data[x+width*y] = 0.0;
			}
			continue;
		}

		for ( dy=-1; dy<=+1; dy++ ) {
		for ( dx=-1; dx<=+1; dx++ ) {
			if ( data[(x+dx)+width*(y+dy)] ) val = 0.0;
		}
		}

		data[x+width*y] = val;

	}
	}
}


void filter_noise(struct image *image)
{
	int i;

	for ( i=0; i<image->detgeom->n_panels; i++ ) {
		struct detgeom_panel *p = &image->detgeom->panels[i];
		filter_noise_in_panel(image->dp[i], p->w, p->h);
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
	long l, m;

	l = 0;
	m = n-1;

	while ( l < m ) {
		long i, j;
		float x;
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
	int pn;

	if ( size <= 0 ) return;

	nn = 2*size+1;
	nn = nn*nn;

	/* "localBg" is way too big, but guaranteed big enough */
	buffer = cfcalloc(nn, sizeof(float));
	if ( buffer == NULL ) {
		ERROR("Failed to allocate LB buffer.\n");
		return;
	}

	/* Determine local background
	 * (median over window width either side of current pixel) */
	for ( pn=0; pn<image->detgeom->n_panels; pn++ ) {

		int fs, ss;
		int i;
		struct detgeom_panel *p;
		float *localBg;

		p = &image->detgeom->panels[pn];

		localBg = cfcalloc(p->w*p->h, sizeof(float));
		if ( localBg == NULL ) {
			ERROR("Failed to allocate LB buffer.\n");
			return;
		}

		for ( ss=0; ss<p->h; ss++ ) {
		for ( fs=0; fs<p->w; fs++ ) {

			int ifs, iss;

			counter = 0;

			// Loop over median window
			for ( iss=-size; iss<=size; iss++ ) {
			for ( ifs=-size; ifs<=size; ifs++ ) {

				int idx;

				if ( (fs+ifs) < 0 ) continue;
				if ( (fs+ifs) >= p->w ) continue;
				if ( (ss+iss) < 0 ) continue;
				if ( (ss+iss) >= p->h ) continue;

				idx = fs+ifs + (ss+iss)*p->w;
				buffer[counter++] = image->dp[pn][idx];

			}
			}

			// Find median value
			localBg[fs+p->w*ss] = kth_smallest(buffer, counter,
			                                   counter/2);

		}
		}

		/* Do the background subtraction */
		for ( i=0; i<p->w*p->h; i++ ) {
			image->dp[pn][i] -= localBg[i];
		}

		cffree(localBg);
	}

	cffree(buffer);
}
