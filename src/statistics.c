/*
 * statistics.c
 *
 * Structure-factor statistics
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdlib.h>

#include "statistics.h"
#include "utils.h"


double stat_scale_intensity(const double *ref1, const unsigned int *c1,
                            const double *ref2, const unsigned int *c2)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<LIST_SIZE; i++ ) {

		if ( c1[i] && c2[i] ) {
			double i1, i2;
			i1 = ref1[i] / (double)c1[i];
			i2 = ref2[i] / (double)c2[i];
			top += i1 * i2;
			bot += i2 * i2;
		} /* else reflection not common so don't worry about it */

	}

	return top/bot;
}


double stat_r2(const double *ref1, const unsigned int *c1,
               const double *ref2, const unsigned int *c2, double *scalep)
{
	double top = 0.0;
	double bot = 0.0;
	double scale;
	int i;
	scale = stat_scale_intensity(ref1, c1, ref2, c2);
	*scalep = scale;

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<LIST_SIZE; i++ ) {

		if ( c1[i] && c2[i] ) {

			double i1, i2;
			i1 = ref1[i] / (scale*(double)c1[i]);
			i2 = ref2[i] / (double)c2[i];

			top += pow(fabs(i1 - i2), 2.0);
			bot += pow(i1, 2.0);

		} /* else reflection not measured so don't worry about it */

	}

	return sqrt(top/bot);
}


double stat_scale_sqrti(const double *ref1, const unsigned int *c1,
                            const double *ref2, const unsigned int *c2)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<LIST_SIZE; i++ ) {

		if ( c1[i] && c2[i] ) {

			double f1, f2;

			if ( (ref1[i]<0.0) || (ref2[i]<0.0) ) continue;

			f1 = sqrt(ref1[i]) / (double)c1[i];
			f2 = sqrt(ref2[i]) / (double)c2[i];

			top += f1 * f2;
			bot += f2 * f2;

		} /* else reflection not common so don't worry about it */

	}

	return top/bot;
}


double stat_rmerge(const double *ref1, const unsigned int *c1,
                   const double *ref2, const unsigned int *c2, double *scalep)
{
	double top = 0.0;
	double bot = 0.0;
	double scale;
	int i;
	scale = stat_scale_sqrti(ref1, c1, ref2, c2);
	*scalep = scale;

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<LIST_SIZE; i++ ) {

		if ( c1[i] && c2[i] ) {

			double f1, f2;

			if ( (ref1[i]<0.0) || (ref2[i]<0.0) ) continue;

			f1 = sqrt(ref1[i]) / (scale*(double)c1[i]);
			f2 = sqrt(ref2[i]) / (double)c2[i];

			top += fabs(f1 - f2);
			bot += f1 + f2;

		} /* else reflection not measured so don't worry about it */

	}

	return 2.0*top/bot;
}
