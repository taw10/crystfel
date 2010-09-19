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
	const double *ref2;
	ReflItemList *items;  /* Which reflections to use */
	int fom;              /* Which FoM to use (see the enum just below) */
};

enum {
	R_1_ZERO,
	R_1_IGNORE,
	R_2,
	R_INT,
	R_DIFF_ZERO,
	R_DIFF_IGNORE,
};


/* Return the least squares optimal scaling factor when comparing intensities.
 * ref1,ref2 are the two intensity lists to compare.  "items" is a ReflItemList
 * containing the reflections which should be taken into account.
 */
double stat_scale_intensity(const double *ref1, const double *ref2,
                            ReflItemList *items)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, i2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = lookup_intensity(ref2, h, k, l);

		top += i1 * i2;
		bot += i2 * i2;

	}

	return top/bot;
}


/* Return the least squares optimal scaling factor when comparing the square
 * roots of the intensities (i.e. one approximation to the structure factor
 * moduli).
 * ref1,ref2 are the two intensity lists to compare (they contain intensities,
 * not square rooted intensities).  "items" is a ReflItemList containing the
 * reflections which should be taken into account.
 */
static double stat_scale_sqrti(const double *ref1, const double *ref2,
                               ReflItemList *items)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, i2, f1, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);
		i2 = lookup_intensity(ref2, h, k, l);
		if ( i2 < 0.0 ) continue;
		f2 = sqrt(i2);

		top += f1 * f2;
		bot += f2 * f2;

	}

	return top/bot;
}


static double internal_r1_ignorenegs(const double *ref1, const double *ref2,
                                     ReflItemList *items, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, f1, i2, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);
		i2 = lookup_intensity(ref2, h, k, l);
		if ( i2 < 0.0 ) continue;
		f2 = sqrt(i2);
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1;

	}

	return top/bot;
}


static double internal_r1_negstozero(const double *ref1, const double *ref2,
                                     ReflItemList *items, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, f1, i2, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;
		i2 = lookup_intensity(ref2, h, k, l);
		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1;

	}

	return top/bot;
}


static double internal_r2(const double *ref1, const double *ref2,
                          ReflItemList *items, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, i2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = scale * lookup_intensity(ref2, h, k, l);

		top += pow(i1 - i2, 2.0);
		bot += pow(i1, 2.0);

	}

	return sqrt(top/bot);
}


static double internal_rint(const double *ref1, const double *ref2,
                          ReflItemList *items, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, i2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = scale * lookup_intensity(ref2, h, k, l);

		top += fabs(i1-i2);
		bot += fabs(i1);

	}

	return top/bot;
}


static double internal_rdiff_negstozero(const double *ref1, const double *ref2,
                                        ReflItemList *items, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, i2, f1, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;
		i2 = lookup_intensity(ref2, h, k, l);
		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1 + f2;

	}

	return 2.0*top/bot;
}


static double internal_rdiff_ignorenegs(const double *ref1, const double *ref2,
                                        ReflItemList *items, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	int i;

	for ( i=0; i<num_items(items); i++ ) {

		double i1, i2, f1, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);
		i2 = lookup_intensity(ref2, h, k, l);
		if ( i2 < 0.0 ) continue;
		f2 = sqrt(i2);
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1 + f2;

	}

	return 2.0*top/bot;
}


static double calc_r(double scale, void *params)
{
	struct r_params *rp = params;

	switch ( rp->fom) {
	case R_1_ZERO :
		return internal_r1_negstozero(rp->ref1, rp->ref2,
		                              rp->items, scale);
	case R_1_IGNORE :
		return internal_r1_ignorenegs(rp->ref1, rp->ref2,
		                              rp->items, scale);
	case R_2 :
		return internal_r2(rp->ref1, rp->ref2, rp->items, scale);

	case R_INT :
		return internal_rint(rp->ref1, rp->ref2, rp->items, scale);

	case R_DIFF_ZERO :
		return internal_rdiff_negstozero(rp->ref1, rp->ref2,
		                                 rp->items, scale);
	case R_DIFF_IGNORE :
		return internal_rdiff_ignorenegs(rp->ref1, rp->ref2,
		                                 rp->items, scale);
	}

	ERROR("No such FoM!\n");
	abort();
}


static double r_minimised(const double *ref1, const double *ref2,
                          ReflItemList *items, double *scalep, int fom)
{
	gsl_function F;
	gsl_min_fminimizer *s;
	int status;
	double scale = 1.0;
	struct r_params rp;
	int iter = 0;

	rp.ref1 = ref1;
	rp.ref2 = ref2;
	rp.items = items;
	rp.fom = fom;

	F.function = &calc_r;
	F.params = &rp;

	s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

	/* Initial guess */
	switch ( fom ) {
	case R_1_ZERO :
	case R_1_IGNORE :
	case R_DIFF_ZERO :
	case R_DIFF_IGNORE :
		scale = stat_scale_sqrti(ref1, ref2, items);
		break;
	case R_2 :
	case R_INT :
		scale = stat_scale_intensity(ref1, ref2, items);
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


double stat_r1_ignore(const double *ref1, const double *ref2,
                      ReflItemList *items, double *scalep)
{
	return r_minimised(ref1, ref2, items, scalep, R_1_IGNORE);
}


double stat_r1_zero(const double *ref1, const double *ref2,
                    ReflItemList *items, double *scalep)
{
	return r_minimised(ref1, ref2, items, scalep, R_1_ZERO);
}


double stat_r2(const double *ref1, const double *ref2,
               ReflItemList *items, double *scalep)
{
	return r_minimised(ref1, ref2, items, scalep, R_2);
}


double stat_rint(const double *ref1, const double *ref2,
                 ReflItemList *items, double *scalep)
{
	return r_minimised(ref1, ref2, items, scalep, R_INT);
}


double stat_rdiff_zero(const double *ref1, const double *ref2,
                       ReflItemList *items, double *scalep)
{
	return r_minimised(ref1, ref2, items, scalep, R_DIFF_ZERO);
}


double stat_rdiff_ignore(const double *ref1, const double *ref2,
                         ReflItemList *items, double *scalep)
{
	return r_minimised(ref1, ref2, items, scalep, R_DIFF_IGNORE);
}


double stat_pearson_i(const double *ref1, const double *ref2,
                      ReflItemList *items)
{
	double *vec1, *vec2;
	int i = 0;
	int ni = num_items(items);
	double val;

	vec1 = malloc(ni*sizeof(double));
	vec2 = malloc(ni*sizeof(double));

	for ( i=0; i<ni; i++ ) {

		double i1, i2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = lookup_intensity(ref2, h, k, l);

		vec1[i] = i1;
		vec2[i] = i2;

		printf("%f %f\n", i1, i2);

	}

	val = gsl_stats_correlation(vec1, 1, vec2, 1, i);
	free(vec1);
	free(vec2);

	return val;
}


double stat_pearson_f_ignore(const double *ref1, const double *ref2,
                             ReflItemList *items)
{
	double *vec1, *vec2;
	int i = 0;
	int ni = num_items(items);
	double val;

	vec1 = malloc(ni*sizeof(double));
	vec2 = malloc(ni*sizeof(double));

	for ( i=0; i<ni; i++ ) {

		double i1, i2, f1, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);
		i2 = lookup_intensity(ref2, h, k, l);
		if ( i2 < 0.0 ) continue;
		f2 = sqrt(i2);

		vec1[i] = f1;
		vec2[i] = f2;

	}

	val = gsl_stats_correlation(vec1, 1, vec2, 1, i);
	free(vec1);
	free(vec2);

	return val;
}


double stat_pearson_f_zero(const double *ref1, const double *ref2,
                      ReflItemList *items)
{
	double *vec1, *vec2;
	int i = 0;
	int ni = num_items(items);
	double val;

	vec1 = malloc(ni*sizeof(double));
	vec2 = malloc(ni*sizeof(double));

	for ( i=0; i<ni; i++ ) {

		double i1, i2, f1, f2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;
		i2 = lookup_intensity(ref2, h, k, l);
		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;

		vec1[i] = f1;
		vec2[i] = f2;

	}

	val = gsl_stats_correlation(vec1, 1, vec2, 1, i);
	free(vec1);
	free(vec2);

	return val;
}
