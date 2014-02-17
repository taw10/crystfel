/*
 * histogram.c
 *
 * Quick histogram functions
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2013 Thomas White <taw@physics.org>
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

#include "histogram.h"

struct _histogram
{
	int n_vals;
	int max_vals;
	double *vals;

	double min;
	double max;

	int full;

	int n_bins;
	double bin_width;

	int *bins;

	int plot_height;
	int peak;
	double mean;
	double stddev;
};


Histogram *histogram_init()
{
	Histogram *hi;

	hi = malloc(sizeof(struct _histogram));
	if ( hi == NULL ) return NULL;

	hi->n_vals = 0;
	hi->max_vals = 0;
	hi->vals = NULL;
	hi->full = 0;

	hi->max = -INFINITY;
	hi->min = INFINITY;

	hi->bins = NULL;
	hi->n_bins = 50;

	hi->plot_height = 10;

	return hi;
}


void histogram_free(Histogram *hi)
{
	free(hi);
}


void histogram_add_value(Histogram *hi, double val)
{
	if ( hi->n_vals == hi->max_vals ) {
		double *new_vals;
		new_vals = realloc(hi->vals, (hi->max_vals+512)*sizeof(double));
		if ( new_vals == NULL ) {
			hi->full = 1;
			fprintf(stderr, "Histogram %p full.\n", hi);
			return;
		}
		hi->max_vals += 512;
		hi->vals = new_vals;
	}

	hi->vals[hi->n_vals++] = val;

	if ( val > hi->max ) hi->max = val;
	if ( val < hi->min ) hi->min = val;
}


static void calc_bins(Histogram *hi)
{
	int i;

	if ( hi->bins != NULL ) {
		free(hi->bins);
	}

	hi->bins = calloc(hi->n_bins, sizeof(int));
	if ( hi->bins == NULL ) {
		fprintf(stderr, "Failed to allocate bins.\n");
		return;
	}

	hi->bin_width = (hi->max - hi->min)/hi->n_bins;

	for ( i=0; i<hi->n_vals; i++ ) {

		int j;

		j = (hi->vals[i] - hi->min) / hi->bin_width;

		/* Tidy up rounding errors */
		if ( j < 0 ) j = 0;
		if ( j >= hi->n_bins ) j = hi->n_bins - 1;

		hi->bins[j]++;

	}

	hi->peak = 0;
	for ( i=0; i<hi->n_bins; i++ ) {
		if ( hi->bins[i] > hi->peak ) hi->peak = hi->bins[i];
	}
}


static void calc_stddev(Histogram *hi)
{
	double sum = 0.0;
	int i;

	for ( i=0; i<hi->n_vals; i++ ) {
		sum += hi->vals[i];
	}
	hi->mean = sum / hi->n_vals;

	sum = 0.0;
	for ( i=0; i<hi->n_vals; i++ ) {
		sum += pow(hi->vals[i] - hi->mean, 2.0);
	}
	hi->stddev = sqrt(sum / hi->n_vals);
}


int *histogram_get_data(Histogram *hi, double *min, double *max, int *n)
{
	calc_bins(hi);

	*n = hi->n_bins;
	*min = hi->min;
	*max = hi->max;

	return hi->bins;
}


void histogram_show(Histogram *hi)
{
	int i;
	int *bins;
	double lines_per_count;
	char tmp[32];
	size_t n_left, n_mid, n_right;

	calc_bins(hi);
	calc_stddev(hi);

	bins = malloc(hi->n_bins*sizeof(int));
	if ( bins == NULL ) {
		fprintf(stderr, "Couldn't scale bins\n");
		return;
	}

	lines_per_count = (double)hi->plot_height / hi->peak;

	for ( i=0; i<hi->n_bins; i++ ) {
		bins[i] = hi->bins[i] * lines_per_count;
	}

	for ( i=hi->plot_height-1; i>=0; i-- ) {

		int j;

		for ( j=0; j<hi->n_bins; j++ ) {
			if ( bins[j] > i ) {
				printf("*");
			} else {
				printf(" ");
			}
		}
		printf("\n");

	}

	printf("|");
	for ( i=1; i<hi->n_bins-1; i++ ) {
		double bin_low, bin_high;
		bin_low = hi->min + i*hi->bin_width;
		bin_high = hi->min + (i+1)*hi->bin_width;
		if ( (bin_low < 0.0) && (bin_high > 0.0) ) {
			printf("+");
		} else {
			printf("-");
		}
	}
	printf("|\n");

	snprintf(tmp, 31, "%.2f", hi->min);
	n_left = strlen(tmp);
	snprintf(tmp, 31, "%.2f", hi->max);
	n_right = strlen(tmp);
	n_mid = hi->n_bins - (n_left + n_right);

	printf("%.2f", hi->min);
	for ( i=0; i<n_mid; i++ ) printf(" ");
	printf("%.2f\n", hi->max);

	printf("Mean = %.2f, Std dev = %.2f\n", hi->mean, hi->stddev);

	free(bins);
}


double histogram_get_min(Histogram *hi)
{
	return hi->min;
}


double histogram_get_max(Histogram *hi)
{
	return hi->max;
}


void histogram_set_min(Histogram *hi, double min)
{
	hi->min = min;
}


void histogram_set_max(Histogram *hi, double max)
{
	hi->max = max;
}


void histogram_set_num_bins(Histogram *hi, int n)
{
	hi->n_bins = n;
}
