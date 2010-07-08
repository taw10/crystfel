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
                       double *ref, double *phases, int zone_axis,
                       UnitCell *cell, unsigned int min_counts)
{
	FILE *fh;
	signed int h, k, l;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "w");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open output file '%s'.\n", filename);
		return;
	}

	/* Write spacings and angle if zone axis pattern */
	if ( zone_axis ) {
		double a, b, c, alpha, beta, gamma;
		cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);
		fprintf(fh, "a %5.3f nm\n",
		            (0.5/resolution(cell, 1, 0, 0))*1e9);
		fprintf(fh, "b %5.3f nm\n",
		            (0.5/resolution(cell, 0, 1, 0))*1e9);
		fprintf(fh, "angle %5.3f deg\n", rad2deg(alpha));
		fprintf(fh, "scale 10\n");
	} else {
		fprintf(fh, "  h   k   l          I    phase   sigma(I) "
		            " 1/d(nm^-1)  counts\n");
	}

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {

		int N;
		double intensity, s;
		char ph[32];

		if ( counts ) {
			N = lookup_count(counts, h, k, l);
			if ( N < min_counts ) continue;
		} else {
			N = 1;
		}
		if ( zone_axis && (l != 0) ) continue;

		/* Divide measured intensity by the number of counts */
		intensity = lookup_intensity(ref, h, k, l) / N;
		if ( phases != NULL ) {
			double p;
			p = lookup_phase(phases, h, k, l);
			snprintf(ph, 31, "%8.6f", p);
		} else {
			strncpy(ph, "       -", 31);
		}

		if ( cell != NULL ) {
			s = 2.0*resolution(cell, h, k, l);
		} else {
			s = 0.0;
		}

		/* h, k, l, I, sigma(I), s */
		fprintf(fh, "%3i %3i %3i %10.2f %s %10.2f  %10.2f %7i\n",
		        h, k, l, intensity, ph, 0.0, s/1.0e9, N);

	}
	}
	}
	fclose(fh);
}


double *read_reflections(const char *filename, unsigned int *counts,
                         double *phases)
{
	double *ref;
	FILE *fh;
	char *rval;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open input file '%s'\n", filename);
		return NULL;
	}

	ref = new_list_intensity();

	do {

		char line[1024];
		signed int h, k, l;
		float intensity, ph, res, sigma;
		char phs[1024];
		int r;
		int cts;

		rval = fgets(line, 1023, fh);
		r = sscanf(line, "%i %i %i %f %s %f %f %i",
		           &h, &k, &l, &intensity, phs, &sigma, &res, &cts);
		if ( r >= 8 ) {
			/* Woohoo */
		} else if ( r >= 7 ) {
			/* No "counts", that's fine.. */
			cts = 1;
		} else if ( r >= 6 ) {
			/* No resolution.  Didn't want it anyway. */
			res = 0.0;
		} else if ( r >= 5 ) {
			/* No sigma.  It's OK today, but one
			 * day I'll get you... */
			sigma = 0.0;
		} else if ( r >= 4 ) {
			/* No phase.  Better not need it.. */
			if ( phases != NULL ) {
				ERROR("Need phases and none were specified!\n");
				abort();
			}
		} else {
			/* You lose. */
			continue;
		}

		set_intensity(ref, h, k, l, intensity);
		if ( phases != NULL ) {
			ph = atof(phs);
			set_phase(phases, h, k, l, ph);
		}
		if ( counts != NULL ) {
			set_count(counts, h, k, l, cts);
			/* In this case, the intensity must be multiplied up
			 * because other parts of the program will try to
			 * divide it down. */
			set_intensity(ref, h, k, l, intensity*(double)cts);
		}

	} while ( rval != NULL );

	fclose(fh);

	return ref;
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


void divide_down(double *intensities, unsigned int *counts)
{
	int i;

	for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
		if ( counts[i] > 0 ) {
			intensities[i] /= (double)counts[i];
			counts[i] = 1;
		}
	}
}
