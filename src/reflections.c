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


void write_reflections(const char *filename, ReflItemList *items,
                       double *intensities, double *phases,
                       unsigned int *counts, UnitCell *cell)
{
	FILE *fh;
	int i;

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
	fprintf(fh, "  h   k   l          I    phase   sigma(I) "
		    " 1/d(nm^-1)  counts\n");

	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		int N;
		double intensity, s, sigma;
		char ph[32];

		if ( counts ) {
			N = lookup_count(counts, h, k, l);
		} else {
			N = 1;
		}

		intensity = lookup_intensity(intensities, h, k, l);

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

		if ( intensity > 0.0 ) {
			sigma = sqrt(intensity);
		} else {
			sigma = 0.0;
		}

		/* h, k, l, I, sigma(I), s */
		fprintf(fh, "%3i %3i %3i %10.2f %s %10.2f  %10.2f %7i\n",
		        h, k, l, intensity, ph, sigma, s/1.0e9, N);

	}
	fclose(fh);
}


/* Read reflections from file.  Returns the list of reflections successfully
 * read in.  "intensities", "phases" and "counts" are lists which will be
 * populated with the values read from the file.  Existing values in either list
 * will be overwritten if the reflection is read from the file, but other values
 * will be left intact.
 *
 * "intensities", "phases" or "counts" can be NULL, if you don't need them.
 */
ReflItemList *read_reflections(const char *filename,
                               double *intensities, double *phases,
                               unsigned int *counts, double *esds)
{
	FILE *fh;
	char *rval;
	ReflItemList *items;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Failed to open input file '%s'\n", filename);
		return NULL;
	}

	items = new_items();

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

		add_item(items, h, k, l);

		if ( intensities != NULL ) {
			set_intensity(intensities, h, k, l, intensity);
		}
		if ( phases != NULL ) {
			ph = atof(phs);
			set_phase(phases, h, k, l, ph);
		}
		if ( counts != NULL ) {
			set_count(counts, h, k, l, cts);
		}
		if ( esds != NULL ) {
			set_sigma(esds, h, k, l, sigma);
		}

	} while ( rval != NULL );

	fclose(fh);

	return items;
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
