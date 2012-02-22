/*
 * statistics.c
 *
 * Structure-factor statistics
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
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

#define _ISOC99_SOURCE
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_statistics.h>

#include "statistics.h"
#include "utils.h"

/**
 * SECTION:statistics
 * @short_description: Intensity statistics and R-factors
 * @title: Statistics
 * @section_id:
 * @see_also:
 * @include: "statistics.h"
 * @Image:
 *
 * These functions are for calculating various figures of merit.
 */


struct r_params {
	RefList *list1;
	RefList *list2;
	int fom;              /* Which FoM to use (see the enum just below) */
};

enum {
	R_1_ZERO,
	R_1_IGNORE,
	R_2,
	R_1_I,
	R_DIFF_ZERO,
	R_DIFF_IGNORE,
	R_DIFF_INTENSITY,
};


/* Return the least squares optimal scaling factor when comparing intensities.
 * list1,list2 are the two intensity lists to compare.
 */
double stat_scale_intensity(RefList *list1, RefList *list2)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		top += i1 * i2;
		bot += i2 * i2;

	}

	return top/bot;
}


/* Return the least squares optimal scaling factor when comparing the square
 * roots of the intensities (i.e. one approximation to the structure factor
 * moduli).
 * list1,list2 are the two intensity lists to compare (they contain intensities,
 * not square rooted intensities).
 */
static double stat_scale_sqrti(RefList *list1, RefList *list2)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);

		if ( i2 < 0.0 ) continue;
		f2 = sqrt(i2);

		top += f1 * f2;
		bot += f2 * f2;

	}

	return top/bot;
}


static double internal_r1_ignorenegs(RefList *list1, RefList *list2,
                                     double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);

		if ( i2 < 0.0 ) continue;
		f2 = sqrt(i2);
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1;

	}

	return top/bot;
}


static double internal_r1_negstozero(RefList *list1, RefList *list2,
                                     double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;

		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1;

	}

	return top/bot;
}


static double internal_r2(RefList *list1, RefList *list2, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		i2 *= scale;

		top += pow(i1 - i2, 2.0);
		bot += pow(i1, 2.0);

	}

	return sqrt(top/bot);
}


static double internal_r_i(RefList *list1, RefList *list2, double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		i2 *= scale;

		top += fabs(i1-i2);
		bot += fabs(i1);

	}

	return top/bot;
}


static double internal_rdiff_intensity(RefList *list1, RefList *list2,
                                       double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		i2 *= scale;

		top += fabs(i1 - i2);
		bot += i1 + i2;

	}

	return 2.0*top/bot;
}


static double internal_rdiff_negstozero(RefList *list1, RefList *list2,
                                        double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;

		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;
		f2 *= scale;

		top += fabs(f1 - f2);
		bot += f1 + f2;

	}

	return 2.0*top/bot;
}


static double internal_rdiff_ignorenegs(RefList *list1, RefList *list2,
                                        double scale)
{
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl1;
	RefListIterator *iter;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		if ( i1 < 0.0 ) continue;
		f1 = sqrt(i1);

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

	switch ( rp->fom ) {
	case R_1_ZERO :
		return internal_r1_negstozero(rp->list1, rp->list2, scale);
	case R_1_IGNORE :
		return internal_r1_ignorenegs(rp->list1, rp->list2, scale);
	case R_2 :
		return internal_r2(rp->list1, rp->list2, scale);

	case R_1_I :
		return internal_r_i(rp->list1, rp->list2, scale);

	case R_DIFF_ZERO :
		return internal_rdiff_negstozero(rp->list1, rp->list2,scale);
	case R_DIFF_IGNORE :
		return internal_rdiff_ignorenegs(rp->list1, rp->list2, scale);
	case R_DIFF_INTENSITY :
		return internal_rdiff_intensity(rp->list1, rp->list2, scale);
	}

	ERROR("No such FoM!\n");
	abort();
}


static double r_minimised(RefList *list1, RefList *list2, double *scalep, int fom,
                          int u)
{
	gsl_function F;
	gsl_min_fminimizer *s;
	int status;
	double scale = 1.0;
	struct r_params rp;
	int iter = 0;

	rp.list1 = list1;
	rp.list2 = list2;
	rp.fom = fom;

	if ( u ) {

		scale = 1.0;

	} else {

		F.function = &calc_r;
		F.params = &rp;

		s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

		/* Initial guess */
		switch ( fom ) {
		case R_1_ZERO :
		case R_1_IGNORE :
		case R_DIFF_ZERO :
		case R_DIFF_IGNORE :
			scale = stat_scale_sqrti(list1, list2);
			break;
		case R_2 :
		case R_1_I :
		case R_DIFF_INTENSITY :
			scale = stat_scale_intensity(list1, list2);
			break;
		}
		//STATUS("Initial scale factor estimate: %5.2e\n", scale);

		/* Probably within an order of magnitude either side */
		gsl_min_fminimizer_set(s, &F, scale, scale/10.0, scale*10.0);

		do {

			double lo, up;

			/* Iterate */
			if ( gsl_min_fminimizer_iterate(s) ) {
				ERROR("Failed to find scale factor.\n");
				return NAN;
			}

			/* Get the current estimate */
			scale = gsl_min_fminimizer_x_minimum(s);
			lo = gsl_min_fminimizer_x_lower(s);
			up = gsl_min_fminimizer_x_upper(s);

			/* Check for convergence */
			status = gsl_min_test_interval(lo, up, 0.001, 0.0);

			iter++;

		} while ( status == GSL_CONTINUE );

		if ( status != GSL_SUCCESS ) {
			ERROR("Scale factor minimisation failed.\n");
		}

		gsl_min_fminimizer_free(s);

	}

	//STATUS("Final scale factor: %5.2e\n", scale);
	*scalep = scale;
	return calc_r(scale, &rp);
}


double stat_r1_ignore(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2, scalep, R_1_IGNORE, u);
}


double stat_r1_zero(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2, scalep, R_1_ZERO, u);
}


double stat_r2(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2, scalep, R_2, u);
}


double stat_r1_i(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2, scalep, R_1_I, u);
}


double stat_rdiff_zero(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2, scalep, R_DIFF_ZERO, u);
}


double stat_rdiff_ignore(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2, scalep, R_DIFF_IGNORE, u);
}


double stat_rdiff_intensity(RefList *list1, RefList *list2, double *scalep, int u)
{
	return r_minimised(list1, list2,  scalep, R_DIFF_INTENSITY, u);
}


double stat_pearson_i(RefList *list1, RefList *list2)
{
	double *vec1, *vec2;
	int ni = num_reflections(list1);
	double val;
	int nacc = 0;
	Reflection *refl1;
	RefListIterator *iter;

	vec1 = malloc(ni*sizeof(double));
	vec2 = malloc(ni*sizeof(double));

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		vec1[nacc] = i1;
		vec2[nacc] = i2;
		nacc++;
	}

	val = gsl_stats_correlation(vec1, 1, vec2, 1, nacc);
	free(vec1);
	free(vec2);

	return val;
}


double stat_pearson_f_ignore(RefList *list1, RefList *list2)
{
	double *vec1, *vec2;
	int ni = num_reflections(list1);
	double val;
	int nacc = 0;
	Reflection *refl1;
	RefListIterator *iter;

	vec1 = malloc(ni*sizeof(double));
	vec2 = malloc(ni*sizeof(double));

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		if ( i1 < 0.0 ) continue;
		if ( i2 < 0.0 ) continue;

		f1 = sqrt(i1);
		f2 = sqrt(i2);

		vec1[nacc] = f1;
		vec2[nacc] = f2;
		nacc++;

	}

	val = gsl_stats_correlation(vec1, 1, vec2, 1, nacc);
	free(vec1);
	free(vec2);

	return val;
}


double stat_pearson_f_zero(RefList *list1, RefList *list2)
{
	double *vec1, *vec2;
	int ni = num_reflections(list1);
	double val;
	int nacc = 0;
	Reflection *refl1;
	RefListIterator *iter;

	vec1 = malloc(ni*sizeof(double));
	vec2 = malloc(ni*sizeof(double));

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		double i1, i2;
		double f1, f2;
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);

		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;
		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;

		vec1[nacc] = f1;
		vec2[nacc] = f2;
		nacc++;

	}

	val = gsl_stats_correlation(vec1, 1, vec2, 1, nacc);
	free(vec1);
	free(vec2);

	return val;
}
