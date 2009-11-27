/*
 * statistics.c
 *
 * Structure-factor statistics
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * integr_sim - Test relrod integration
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdlib.h>

#include "statistics.h"


/* By what (best-fitted) factor must the list "list" be multiplied by,
 * if it were to be merged with "target"? */
static double stat_scale_intensity(double *obs, double *calc, unsigned int *c,
                                   int size)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<size; i++ ) {

		if ( c[i] > 0 ) {
			double obsi;
			obsi = obs[i] / (double)c[i];
			top += obsi * calc[i];
			bot += calc[i] * calc[i];
		} /* else reflection not measured so don't worry about it */

	}

	return top/bot;
}


double stat_r2(double *obs, double *calc, unsigned int *c, int size,
               double *scalep)
{
	double top = 0.0;
	double bot = 0.0;
	double scale;
	int i;
	scale = stat_scale_intensity(obs, calc, c, size);
	*scalep = scale;

	for ( i=0; i<size; i++ ) {

		if ( c[i] > 0 ) {
			double obsi;
			obsi = obs[i] / (double)c[i];
			top += fabs(obsi/scale - calc[i]);
			bot += obsi/scale;
		}

	} /* else reflection not measured so don't worry about it */

	return top/bot;
}
