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

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<size; i++ ) {

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

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<size; i++ ) {

		if ( c[i] > 0 ) {

			double obsi;

			obsi = obs[i] / (double)c[i];
			obsi = obsi / scale;

			top += pow(fabs(obsi - calc[i]), 2.0);
			bot += pow(obsi, 2.0);

		} /* else reflection not measured so don't worry about it */

	}

	return sqrt(top/bot);
}
