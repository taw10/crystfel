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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_statistics.h>

#include "statistics.h"
#include "utils.h"


struct r_params {
	const double *ref1;
	const unsigned int *c1;
	const double *ref2;
	const unsigned int *c2;
	int fom;
};

enum {
	R_2,
	R_MERGE,
};


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


static double internal_r2(const double *ref1, const unsigned int *c1,
                          const double *ref2, const unsigned int *c2,
                          double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	/* Start from i=1 -> skip central beam */
	for ( i=1; i<LIST_SIZE; i++ ) {

		if ( c1[i] && c2[i] ) {

			double i1, i2;
			i1 = ref1[i] / (scale*(double)c1[i]);
			i2 = ref2[i] / (double)c2[i];

			top += pow(i1 - i2, 2.0);
			bot += pow(i1, 2.0);

		} /* else reflection not measured so don't worry about it */

	}

	return sqrt(top/bot);
}


static double internal_rmerge(const double *ref1, const unsigned int *c1,
                              const double *ref2, const unsigned int *c2,
                              double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

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


static double calc_r(double scale, void *params)
{
	struct r_params *rp = params;

	switch ( rp->fom) {
	case R_MERGE :
		return internal_rmerge(rp->ref1, rp->c1,
		                       rp->ref2, rp->c2, scale);
	case R_2 :
		return internal_r2(rp->ref1, rp->c1,
		                   rp->ref2, rp->c2, scale);
	}

	ERROR("No such FoM!\n");
	abort();
}


static double r_minimised(const double *ref1, const unsigned int *c1,
                          const double *ref2, const unsigned int *c2,
                          double *scalep, int fom)
{
	gsl_function F;
	gsl_min_fminimizer *s;
	int status;
	double scale = 1.0;
	struct r_params rp;
	int iter = 0;

	rp.ref1 = ref1;
	rp.ref2 = ref2;
	rp.c1 = c1;
	rp.c2 = c2;
	rp.fom = fom;

	F.function = &calc_r;
	F.params = &rp;

	s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

	/* Initial guess */
	switch ( fom ) {
	case R_MERGE :
		scale = stat_scale_sqrti(ref1, c1, ref2, c2);
		break;
	case R_2 :
		scale = stat_scale_intensity(ref1, c1, ref2, c2);
		break;
	}
	//STATUS("Initial scale factor estimate: %5.2e\n", scale);

	/* Probably within an order of magnitude either side */
	gsl_min_fminimizer_set(s, &F, scale, scale/10.0, scale*10.0);

	do {

		double lo, up;

		/* Iterate */
		gsl_min_fminimizer_iterate(s);

		/* Get the current estimate */
		scale = gsl_min_fminimizer_x_minimum(s);
		lo = gsl_min_fminimizer_x_lower(s);
		up = gsl_min_fminimizer_x_upper(s);

		/* Check for convergence */
		status = gsl_min_test_interval(lo, up, 0.001, 0.0);

		iter++;

	} while ( status == GSL_CONTINUE );

	if (status != GSL_SUCCESS) {
		ERROR("Scale factor minimisation failed.\n");
	}

	gsl_min_fminimizer_free(s);

	//STATUS("Final scale factor: %5.2e\n", scale);
	*scalep = scale;
	return calc_r(scale, &rp);
}


double stat_rmerge(const double *ref1, const unsigned int *c1,
                   const double *ref2, const unsigned int *c2,
                   double *scalep)
{
	return r_minimised(ref1, c1, ref2, c2, scalep, R_MERGE);
}


double stat_r2(const double *ref1, const unsigned int *c1,
               const double *ref2, const unsigned int *c2,
               double *scalep)
{
	return r_minimised(ref1, c1, ref2, c2, scalep, R_2);
}


double stat_pearson(const double *ref1, const unsigned int *c1,
                    const double *ref2, const unsigned int *c2)
{
	double vec1[4096];
	double vec2[4096];
	signed int h, k, l;
	int i = 0;

	for ( l=-INDMAX; l<INDMAX; l++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( h=-INDMAX; h<INDMAX; h++ ) {

		double i1, i2;
		unsigned int c1s, c2s;

		c1s = lookup_count(c1, h, k, l);
		c2s = lookup_count(c2, h, k, l);

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = lookup_intensity(ref2, h, k, l);

		if ( c1s && c2s && (i1>0.0) && (i2>0.0) ) {
			vec1[i] = sqrt(i1 / (double)c1s);
			vec2[i] = sqrt(i2 / (double)c2s);
			i++;
		}

	}
	}
	}

	return gsl_stats_correlation(vec1, 1, vec2, 1, i);
}
