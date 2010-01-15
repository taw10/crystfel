/*
 * reflections.c
 *
 * Utilities for handling reflections
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "utils.h"
#include "cell.h"
#include "reflections.h"


void write_reflections(const char *filename, unsigned int *counts,
                       double *ref, int zone_axis, UnitCell *cell)
{
	FILE *fh;
	signed int h, k, l;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "w");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open output file!\n");
		return;
	}

	/* Write spacings and angle if zone axis pattern */
	if ( zone_axis ) {
		double a, b, c, alpha, beta, gamma;
		cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);
		fprintf(fh, "a %5.3f nm\n", a*1e9);
		fprintf(fh, "b %5.3f nm\n", b*1e9);
		fprintf(fh, "angle %5.3f deg\n", rad2deg(gamma));
		fprintf(fh, "scale 10\n");
	}

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {

		int N;
		double F;

		if ( counts ) {
			N = lookup_count(counts, h, k, l);
			if ( N == 0 ) continue;
		} else {
			N = 1;
		}

		F = lookup_intensity(ref, h, k, l) / N;
		if ( zone_axis && (l != 0) ) continue;

		fprintf(fh, "%3i %3i %3i %f\n", h, k, l, F);

	}
	}
	}
	fclose(fh);
}


double *ideal_intensities(double complex *sfac)
{
	double *ref;
	signed int h, k, l;

	ref = new_list_intensity();

	/* Generate ideal reflections from complex structure factors */
	for ( h=-INDMAX; h<=INDMAX; h++ ) {
	for ( k=-INDMAX; k<=INDMAX; k++ ) {
	for ( l=-INDMAX; l<=INDMAX; l++ ) {
		double complex F = lookup_sfac(sfac, h, k, l);
		double intensity = pow(cabs(F), 2.0);
		set_intensity(ref, h, k, l, intensity);
	}
	}
	}

	return ref;
}
