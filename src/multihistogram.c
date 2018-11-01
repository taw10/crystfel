/*
 * multimultihistogram.c
 *
 * MultiHistogram with categories
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2013-2014 Thomas White <taw@physics.org>
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

#include "multihistogram.h"

struct _multihistogram
{
	double min;
	double max;

	int n_bins;
	double bin_width;

	int *bins[32];
};


void multihistogram_delete_all_values(MultiHistogram *hi)
{
	int i;

	for ( i=0; i<32; i++ ) {
		if ( hi->bins[i] != NULL ) free(hi->bins[i]);
	}

	for ( i=0; i<32; i++ ) {
		hi->bins[i] = calloc(hi->n_bins, sizeof(int));
	}
}


MultiHistogram *multihistogram_new()
{
	MultiHistogram *hi;
	int i;

	hi = malloc(sizeof(struct _multihistogram));
	if ( hi == NULL ) return NULL;

	hi->max = -INFINITY;
	hi->min = INFINITY;

	hi->n_bins = 50;

	for ( i=0; i<32; i++ ) hi->bins[i] = NULL;

	return hi;
}


void multihistogram_free(MultiHistogram *hi)
{
	int i;
	for ( i=0; i<32; i++ ) {
		if ( hi->bins[i] != NULL ) free(hi->bins[i]);
	}
	free(hi);
}


void multihistogram_add_value(MultiHistogram *hi, double val, unsigned int cat)
{
	int i, j;

	j = (val - hi->min) / hi->bin_width;

	/* Tidy up rounding errors */
	if ( j < 0 ) j = 0;
	if ( j >= hi->n_bins ) j = hi->n_bins - 1;

	for ( i=0; i<32; i++ ) {
		if ( cat & (unsigned)1<<i ) hi->bins[i][j]++;
	}
}


int *multihistogram_get_data(MultiHistogram *hi, int cat)
{
	if ( cat < 0 ) return NULL;
	if ( cat > 31 ) return NULL;
	return hi->bins[cat];
}


void multihistogram_set_min(MultiHistogram *hi, double min)
{
	hi->min = min;
	hi->bin_width = (hi->max - hi->min)/hi->n_bins;
	multihistogram_delete_all_values(hi);
}


void multihistogram_set_max(MultiHistogram *hi, double max)
{
	hi->max = max;
	hi->bin_width = (hi->max - hi->min)/hi->n_bins;
	multihistogram_delete_all_values(hi);
}


void multihistogram_set_num_bins(MultiHistogram *hi, int n)
{
	hi->n_bins = n;
	hi->bin_width = (hi->max - hi->min)/hi->n_bins;
	multihistogram_delete_all_values(hi);
}
