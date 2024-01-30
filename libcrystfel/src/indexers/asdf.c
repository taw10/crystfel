/*
 * asdf.c
 *
 * Alexandra's Superior Direction Finder, or
 * Algorithm Similar to DirAx, FFT-based
 *
 * Copyright © 2014-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014-2015 Alexandra Tolstikova <alexandra.tolstikova@desy.de>
 *   2015-2017 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <profile.h>
#include <assert.h>

#include "index.h"
#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell-utils.h"
#include "asdf.h"
#include "cell.h"

/**
 * \file asdf.h
 */

#ifdef HAVE_FFTW

#define FFTW_NO_Complex  /* Please use "double[2]", not C99 "complex",
                          * despite complex.h possibly already being
                          * included.  For more information, refer to:
                          * http://www.fftw.org/doc/Complex-numbers.html */

#include <fftw3.h>

struct fftw_vars {
	int N;
	fftw_plan p;
	double *in;
	fftw_complex *out;
};


struct asdf_private {
	IndexingMethod          indm;
	UnitCell                *template;
	struct fftw_vars        fftw;
	int                     fast_execution;
};


/* Possible direct vector */
struct tvector {
	gsl_vector *t;
	int n; // number of fitting reflections
	int *fits;
};


struct fftw_vars fftw_vars_new()
{
	struct fftw_vars fftw;
	int N = 1024;

	fftw.N = N;
	fftw.in = fftw_alloc_real(N);
	fftw.out = fftw_alloc_complex(N);
	fftw.p = fftw_plan_dft_r2c_1d(N, fftw.in, fftw.out, FFTW_MEASURE);

	return fftw;
}


static void fftw_vars_free(struct fftw_vars fftw)
{
	fftw_free(fftw.in);
	fftw_free(fftw.out);
	fftw_destroy_plan(fftw.p);
}


struct asdf_cell {
	gsl_vector *axes[3];
	gsl_vector *reciprocal[3];

	int n; // number of fitting reflections
	double volume;

	int N_refls; // total number of reflections
	int *reflections; // reflections[i] = 1 if reflections fits
	double **indices; // indices[i] = [h, k, l] for all reflections (not rounded)

	int acl; // minimum number of reflections fitting to one of the axes[]
	int n_max; // maximum number of reflections fitting to some t-vector
};


struct tvector tvector_new(int n)
{
	struct tvector t;

	t.t = gsl_vector_alloc(3);
	t.n = 0;
	t.fits = cfmalloc(sizeof(int) * n);

	return t;
}


static int tvector_free(struct tvector t)
{
	gsl_vector_free(t.t);
	cffree(t.fits);

	return 1;
}


static int asdf_cell_free(struct asdf_cell *c)
{
	int i;
	for ( i = 0; i < 3; i++ ) {
		gsl_vector_free(c->axes[i]);
		gsl_vector_free(c->reciprocal[i]);
	}

	cffree(c->reflections);
	for ( i = 0; i < c->N_refls; i++ ) {
		cffree(c->indices[i]);
	}
	cffree(c->indices);
	cffree(c);

	return 1;
}


static struct asdf_cell *asdf_cell_new(int n)
{
	struct asdf_cell *c;
	c = cfmalloc(sizeof(struct asdf_cell));

	int i;
	for ( i = 0; i < 3; i++ ) {
		c->axes[i] = gsl_vector_alloc(3);
		c->reciprocal[i] = gsl_vector_alloc(3);
	}

	c->N_refls = n;
	c->reflections = cfmalloc(sizeof(int) * n);
	if (c->reflections == NULL) return NULL;

	c->indices = cfmalloc(sizeof(double *) * n);
	if (c->indices == NULL) {
		cffree(c->reflections);
		return NULL;
	}

	for ( i = 0; i < n; i++ ) {
		c->indices[i] = cfmalloc(sizeof(double) * 3);
		if (c->indices[i] == NULL) {
			cffree(c->reflections);
			cffree(c->indices);
			return NULL;
		}
	}

	c->n = 0;

	c->acl = 0;
	c->n_max = 0;

	return c;
}


static int asdf_cell_memcpy(struct asdf_cell *dest, struct asdf_cell *src)
{
	int i;
	for ( i = 0; i < 3; i++ ) {
		gsl_vector_memcpy(dest->axes[i], src->axes[i]);
		gsl_vector_memcpy(dest->reciprocal[i], src->reciprocal[i]);
	}

	dest->volume = src->volume;

	int n = src->N_refls;

	dest->n = src->n;
	memcpy(dest->reflections, src->reflections, sizeof(int) * n);

	for (i  = 0; i < n; i++ ) {
		memcpy(dest->indices[i], src->indices[i], sizeof(double) * 3);
	}
	for ( i=n; i<dest->N_refls; i++ ) {
		cffree(dest->indices[i]);
	}
	dest->N_refls = n;

	dest->acl = src->acl;
	dest->n_max = src->n_max;
	return 1;
}


/* result = vec1 cross vec2 */
static int cross_product(gsl_vector *vec1, gsl_vector *vec2,
                         gsl_vector **result)
{
	double c1[3], c2[3], p[3];
	int i;
	for ( i = 0; i < 3; i++ ) {
		c1[i] = gsl_vector_get(vec1, i);
		c2[i] = gsl_vector_get(vec2, i);
	}

	p[0] = c1[1] * c2[2] - c1[2] * c2[1];
	p[1] = - c1[0] * c2[2] + c1[2] * c2[0];
	p[2] = c1[0] * c2[1] - c1[1] * c2[0];

	for ( i = 0; i < 3; i++ ) {
		gsl_vector_set(*result, i, p[i]);
	}

	return 1;
}


/* Returns triple product of three gsl_vectors */
static double calc_volume(gsl_vector *vec1, gsl_vector *vec2, gsl_vector *vec3)
{
	double volume;
	gsl_vector *cross = gsl_vector_alloc(3);

	cross_product(vec1, vec2, &cross);
	gsl_blas_ddot(vec3, cross, &volume);

	gsl_vector_free(cross);
	return volume;
}


static int calc_reciprocal(gsl_vector **direct, gsl_vector **reciprocal)
{
	double volume;

	cross_product(direct[1], direct[2], &reciprocal[0]);
	gsl_blas_ddot(direct[0], reciprocal[0], &volume);
	gsl_vector_scale(reciprocal[0], 1/volume);

	cross_product(direct[2], direct[0], &reciprocal[1]);
	gsl_vector_scale(reciprocal[1], 1/volume);

	cross_product(direct[0], direct[1], &reciprocal[2]);
	gsl_vector_scale(reciprocal[2], 1/volume);

	return 1;
}


static int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}


static double max(double a, double b, double c)
{
     double m = a;
     if ( m < b ) m = b;
     if ( m < c ) m = c;
     return m;
}


/* Compares tvectors by length */
static int compare_tvectors(const void *a, const void *b)
{
	struct tvector *ta = (struct tvector *) a;
	struct tvector *tb = (struct tvector *) b;

	//~ if (ta->n == tb->n) {
	return (gsl_blas_dnrm2(ta->t) > gsl_blas_dnrm2(tb->t)) -
	       (gsl_blas_dnrm2(ta->t) < gsl_blas_dnrm2(tb->t));
	//~ }
	//~
	//~ return (ta->n > tb->n) - (ta->n < tb->n);
}


/* Calculates normal to a triplet c1, c2, c3. Returns 0 if reflections are on
 * the same line */
static int calc_normal(gsl_vector *c1, gsl_vector *c2, gsl_vector *c3,
                       gsl_vector *normal)
{
	gsl_vector *c12 = gsl_vector_alloc(3);
	gsl_vector *c23 = gsl_vector_alloc(3);
	gsl_vector *c31 = gsl_vector_alloc(3);
	gsl_vector *res = gsl_vector_alloc(3);

	cross_product(c1, c2, &c12);
	cross_product(c2, c3, &c23);
	cross_product(c3, c1, &c31);

	int i;
	for ( i = 0; i < 3; i++ ) {
		gsl_vector_set(res, i, gsl_vector_get(c12, i) +
		                       gsl_vector_get(c23, i) +
		                       gsl_vector_get(c31, i));
	}

	gsl_vector_free(c12);
	gsl_vector_free(c23);
	gsl_vector_free(c31);

	double norm = gsl_blas_dnrm2(res);
	if ( norm < 0.0001 ) {
		gsl_vector_free(res);
		return 0;
	} else {
		gsl_vector_scale(res, 1/norm);
		gsl_vector_memcpy(normal, res);
		gsl_vector_free(res);
	}

	return 1;
}


static float find_ds_fft(double *projections, int N_projections, double d_max,
                         struct fftw_vars fftw)
{
	int n = N_projections;
	double projections_sorted[n];
	memcpy(projections_sorted, projections, sizeof(double) * n);
	qsort(projections_sorted, n, sizeof(double), compare_doubles);

	int i;

	int N = fftw.N; // number of points in fft calculation

	for ( i=0; i<N; i++ ) {
		fftw.in[i] = 0;
	}

	for ( i=0; i<n; i++ ) {
		int k;
		k = (int)((projections_sorted[i] - projections_sorted[0]) /
		          (projections_sorted[n - 1] - projections_sorted[0]) *
		          (N - 1));
		if ( (k>=N) || (k<0) ) {
			ERROR("Bad k value in find_ds_fft() (k=%i, N=%i)\n",
			      k, N);
			return -1.0;
		}
		fftw.in[k]++;
	}

	fftw_execute(fftw.p);

	int i_max = (int)(d_max * (projections_sorted[n - 1] -
	                           projections_sorted[0]));

	if ( i_max > N / 2 ) i_max = N / 2;

	int d = 1;
	double maxval = 0;
	for ( i=1; i<=i_max; i++ ) {
		double a;
		a = sqrt(fftw.out[i][0] * fftw.out[i][0] + fftw.out[i][1] * fftw.out[i][1]);
		if (a > maxval) {
			maxval = a;
			d = i;
		}
	}

	double ds = (projections_sorted[n - 1] - projections_sorted[0]) / d;

	return ds;
}


/* Returns number of reflections fitting ds.
 * A projected reflection fits a one-dimensional lattice with elementary
 * lattice vector d* if its absolute distance to the nearest lattice
 * point is less than LevelFit. */
static int check_refl_fitting_ds(double *projections, int N_projections,
                                 double ds, double LevelFit)
{
	if ( ds == 0 ) return 0;

	int i;
	int n = 0;
	for ( i = 0; i < N_projections; i++ ) {
		if ( fabs(projections[i] -
		     ds * round(projections[i]/ds)) < LevelFit )
		{
			n += 1;
		}
	}

	return n;
}


/* Refines d*, writes 1 to fits[i] if the i-th projection fits d* */
static float refine_ds(double *projections, int N_projections, double ds,
                       double LevelFit, int *fits)
{
	double fit_refls[N_projections];
	double indices[N_projections];

	int i;

	int N_fits = 0;
	int N_new = N_projections;

	double c1, cov11, sumsq;
	double ds_new = ds;
	while ( N_fits < N_new ) {
		N_fits = 0;
		for ( i = 0; i < N_projections; i++ ) {
			if ( fabs(projections[i] - ds_new *
			     round(projections[i] / ds_new)) < LevelFit )
			{
				fit_refls[N_fits] = projections[i];
				indices[N_fits] = round(projections[i]/ds_new);
				N_fits ++;
				fits[i] = 1;
			} else {
				fits[i] = 0;
			}
		}


		gsl_fit_mul(indices, 1, fit_refls, 1, N_fits, &c1, &cov11,
		            &sumsq);
		N_new = check_refl_fitting_ds(projections, N_projections, c1,
		                              LevelFit);
		if ( N_new >= N_fits ) {
			ds_new = c1;
		}
	}

	return ds_new;
}


static int check_refl_fitting_cell(struct asdf_cell *c,
                                   gsl_vector **reflections,
                                   int N_reflections, double IndexFit)
{
	double dist[3];

	calc_reciprocal(c->axes, c->reciprocal);
	c->n = 0;
	int i, j, k;
	for( i = 0; i < N_reflections; i += 1 ) {

		for ( j = 0; j < 3; j++ ) dist[j] = 0;

		for ( j = 0; j < 3; j++ ) {
			gsl_blas_ddot(reflections[i], c->axes[j],
			              &c->indices[i][j]);

			for ( k = 0; k < 3; k++ ) {
				dist[k] += gsl_vector_get(c->reciprocal[j], k) *
				(c->indices[i][j] - round(c->indices[i][j]));
			}
		}

		/* A reflection fits if the distance (in reciprocal space)
		 * between the observed and calculated reflection position
		 * is less than Indexfit */

		if ( dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2] <
							 IndexFit * IndexFit ) {
			c->reflections[i] = 1;
			c->n++;
		} else {
			c->reflections[i] = 0;
		}
	}

	return 1;
}


/* Returns 0 when refinement doesn't converge (i.e. all fitting reflections
 * lie in the same plane) */
static int refine_asdf_cell(struct asdf_cell *c, gsl_vector **reflections,
                            int N_reflections, double IndexFit)
{
	gsl_matrix *X = gsl_matrix_alloc(c->n, 3);

	gsl_vector *r[] = {gsl_vector_alloc(c->n),
	                   gsl_vector_alloc(c->n),
	                   gsl_vector_alloc(c->n)};

	gsl_vector *res = gsl_vector_alloc(3);
	gsl_matrix *cov = gsl_matrix_alloc (3, 3);
	double chisq;

	int i, j;
	int n = 0;
	for ( i = 0; i < N_reflections; i++ ) if ( c->reflections[i] == 1 )
	{
		for ( j = 0; j < 3; j++ ) {
			gsl_matrix_set(X, n, j, round(c->indices[i][j]));
			gsl_vector_set(r[j], n,
			               gsl_vector_get(reflections[i], j));
		}
		n++;
	}

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(c->n,
	                                                                3);

	for ( i = 0; i < 3; i++ ) {
		gsl_multifit_linear (X, r[i], res, cov, &chisq, work);

		for (j = 0; j < 3; j++ ) {
			gsl_vector_set(c->reciprocal[j], i,
			               gsl_vector_get(res, j));
		}
	}

	calc_reciprocal(c->reciprocal, c->axes);

	double a[3];
	for ( i = 0; i < 3; i++ ) {
		a[i] = gsl_blas_dnrm2(c->axes[i]);
	}

	gsl_multifit_linear_free(work);
	gsl_vector_free(res);
	gsl_matrix_free(cov);
	gsl_matrix_free(X);
	for ( i = 0; i < 3; i++ ) {
		gsl_vector_free(r[i]);
	}

	if ( fabs(a[0]) > 10000 || fabs(a[1]) > 10000 ||
	     fabs(a[2]) > 10000 || isnan(a[0]) )
	{
		return 0;
	}

	return 1;
}


static int reduce_asdf_cell(struct asdf_cell *cl)
{
	gsl_vector *va = gsl_vector_alloc(3);
	gsl_vector *vb = gsl_vector_alloc(3);
	gsl_vector *vc = gsl_vector_alloc(3);

	int changed = 1;

	int n = 0;
	while ( changed ) {

		double a, b, c, alpha, beta, gamma, ab, bc, ca, bb;

		n += 1;
		changed = 0;

		gsl_vector_memcpy(va, cl->axes[0]);
		gsl_vector_memcpy(vb, cl->axes[1]);
		gsl_vector_memcpy(vc, cl->axes[2]);

		a = gsl_blas_dnrm2(va);
		b = gsl_blas_dnrm2(vb);
		c = gsl_blas_dnrm2(vc);

		gsl_blas_ddot(va, vb, &ab);
		gsl_blas_ddot(vb, vc, &bc);
		gsl_blas_ddot(vc, va, &ca);

		alpha = acos(bc/b/c)/M_PI*180;
		beta = acos(ca/a/c)/M_PI*180;
		gamma = acos(ab/a/b)/M_PI*180;

		if ( gamma < 90 ) {
			gsl_vector_scale(vb, -1);
			gamma = 180 - gamma;
			alpha = 180 - alpha;
		}

		gsl_vector_add(vb, va);
		bb = gsl_blas_dnrm2(vb);
		if ( bb < b ) {
			b = bb;
			if ( a < b ) {
				gsl_vector_memcpy(cl->axes[1], vb);
			} else {
				gsl_vector_memcpy(cl->axes[1], va);
				gsl_vector_memcpy(cl->axes[0], vb);
			}
			changed = 1;
		}

		if ( changed == 0 ) {

			double cc;

			if ( beta < 90 ) {
				gsl_vector_scale(vc, -1);
				beta = 180 - beta;
				alpha = 180 - alpha;
			}

			gsl_vector_add(vc, va);
			cc = gsl_blas_dnrm2(vc);
			if ( cc < c ) {
				c = cc;
				if ( b < c ) {
					gsl_vector_memcpy(cl->axes[2], vc);
				} else if ( a < c ) {
					gsl_vector_memcpy(cl->axes[1], vc);
					gsl_vector_memcpy(cl->axes[2], vb);
				} else {
					gsl_vector_memcpy(cl->axes[0], vc);
					gsl_vector_memcpy(cl->axes[1], va);
					gsl_vector_memcpy(cl->axes[2], vb);
				}
				changed = 1;
			}
		}

		if ( changed == 0 ) {

			double cc;

			if ( alpha < 90 ) {
				gsl_vector_scale(vc, -1);
				beta = 180 - beta;
				alpha = 180 - alpha;
			}

			gsl_vector_add(vc, vb);
			cc = gsl_blas_dnrm2(vc);
			if ( cc < c ) {
				c = cc;
				if ( b < c ) {
					gsl_vector_memcpy(cl->axes[2], vc);
				} else if ( a < c ) {
					gsl_vector_memcpy(cl->axes[1], vc);
					gsl_vector_memcpy(cl->axes[2], vb);
				} else {
					gsl_vector_memcpy(cl->axes[0], vc);
					gsl_vector_memcpy(cl->axes[1], va);
					gsl_vector_memcpy(cl->axes[2], vb);
				}
				changed = 1;
			}
		}

		if (n > 30) changed = 0;
	}

	cross_product(cl->axes[0], cl->axes[1], &vc);
	gsl_blas_ddot(vc, cl->axes[2], &cl->volume);
	if ( cl->volume < 0 ) {
		gsl_vector_scale(cl->axes[2], -1);
		cl->volume *= -1;
	}

	gsl_vector_free(va);
	gsl_vector_free(vb);
	gsl_vector_free(vc);

	return 1;
}


static int check_cell_angles(gsl_vector *va, gsl_vector *vb, gsl_vector *vc,
                             double max_cos)
{
	double a, b, c, cosa, cosb, cosg, ab, bc, ca;

	a = gsl_blas_dnrm2(va);
	b = gsl_blas_dnrm2(vb);
	c = gsl_blas_dnrm2(vc);

	gsl_blas_ddot(va, vb, &ab);
	gsl_blas_ddot(vb, vc, &bc);
	gsl_blas_ddot(vc, va, &ca);

	cosa = bc/b/c;
	cosb = ca/a/c;
	cosg = ab/a/b;

	if ( fabs(cosa) > max_cos || fabs(cosb) > max_cos ||
	                             fabs(cosg) > max_cos ) {
		return 0;
	}

	return 1;
}


/* Returns min(t1.n, t2.n, t3.n) */
static int find_acl(struct tvector t1, struct tvector t2, struct tvector t3)
{
	int i = t1.n, j = t2.n, k = t3.n;
	if ( i <= j && i <= k ) return i;
	if ( j <= i && j <= k ) return j;
	if ( k <= i && k <= j ) return k;
	ERROR("This point never reached!\n");
	abort();
}


static int create_cell(struct tvector tvec1, struct tvector tvec2,
                       struct tvector tvec3, struct asdf_cell *c,
                       double IndexFit, double volume_min, double volume_max,
                       gsl_vector **reflections, int N_reflections)
{

	double volume = calc_volume(tvec1.t, tvec2.t, tvec3.t);
	if ( fabs(volume) < volume_min || fabs(volume) > volume_max ) return 0;

	gsl_vector_memcpy(c->axes[0], tvec1.t);
	gsl_vector_memcpy(c->axes[1], tvec2.t);
	gsl_vector_memcpy(c->axes[2], tvec3.t);

	c->volume = volume;
	check_refl_fitting_cell(c, reflections, N_reflections, IndexFit);

	if ( c->n < 6 ) return 0;

	reduce_asdf_cell(c);

	/* If one of the cell angles > 135 or < 45 return 0 */
	if ( !check_cell_angles(c->axes[0], c->axes[1],
	                        c->axes[2], 0.71) ) return 0;

	/* Index reflections with new cell axes */
	check_refl_fitting_cell(c, reflections, N_reflections, IndexFit);

	/* Refine cell until the number of fitting
	 * reflections stops increasing */
	int n = 0;
	int cell_correct = 1;
	while ( c->n - n > 0 && cell_correct ) {

		n = c->n;
		cell_correct = refine_asdf_cell(c, reflections, N_reflections,
					        IndexFit);
		check_refl_fitting_cell(c, reflections, N_reflections,
					IndexFit);
	}

	return cell_correct;
}


static int find_cell(struct tvector *tvectors, int N_tvectors, double IndexFit,
                     double volume_min, double volume_max, int n_max,
                     gsl_vector **reflections, int N_reflections,
                     struct asdf_cell *result)
{
	int i, j, k;

	/* Only tvectors with the number of fitting reflections > acl are
	 * considered */
	int acl = N_reflections < 18 ? 6 : N_reflections/3;

	struct asdf_cell *c = asdf_cell_new(N_reflections);
	if (c == NULL) {
		ERROR("Failed to allocate asdf_cell in find_cell!\n");
		return 0;
	}

	/* Traversing a 3d array in slices perpendicular to the main diagonal */
	int sl;
	for ( sl = 0; sl < 3 * N_tvectors - 1; sl++ ) {

		int i_min = sl < 2 * N_tvectors ? 0 : sl - 2 * N_tvectors;
		int i_max = sl < N_tvectors ? sl : N_tvectors;

		for ( i = i_min; i < i_max; i++) {

			if (tvectors[i].n <= acl ) continue;

			int j_min = sl - N_tvectors - 2 * i - 1 < 0 ?
			                            i + 1 : sl - N_tvectors - i;
			int j_max = sl - N_tvectors - i < 0 ?
			                                    sl - i : N_tvectors;

			for ( j = j_min; j < j_max; j++ ) {

				if ( tvectors[j].n <= acl ) continue;

				k = sl - i - j - 1;

				if ( k > j && tvectors[k].n > acl ) {

					if ( !create_cell(tvectors[i],
						          tvectors[j],
						          tvectors[k],
						          c, IndexFit,
						          volume_min,
						          volume_max,
						          reflections,
						          N_reflections) )
					{
						break;
					}

					acl = find_acl(tvectors[i],
					               tvectors[j],
					               tvectors[k]);
					c->acl = acl;
					c->n_max = n_max;

					/* If the new cell has more fitting
					 * reflections save it to result */
					if ( result->n < c->n ) {
						asdf_cell_memcpy(result, c);
					}
					acl++;

					if ( acl > n_max ) break;
					if ( tvectors[j].n <= acl ||
					     tvectors[i].n <= acl ) break;
				}
			}
			if ( acl > n_max ) break;
			if ( tvectors[i].n <= acl ) break;
		}
		if ( acl > n_max ) break;
	}

	asdf_cell_free(c);

	if ( result->n ) return 1;
	return 0;
}


void swap(int *a, int *b) {
	int c = *a;
	*a = *b;
	*b = c;
}


long nCk(int n, int k) {
	// computes nCk, the number of combinations n choose k for 0 <= k <= 3
	assert(k>=0 && k<4);
	switch ( k ) {
		case 0 : return 1;
		case 1 : return (long)n;
		case 2 : return (long)n*(n-1)/2;
		case 3 : return (long)n*(n-1)*(n-2)/6;
	}
	return 0;
}


void get_triplet_by_index(int index, int n, int *triplet) {
	// finds the i-th combination of 3 numbers chosen from 0,1,2,...,n-1
	int r = index + 1;
	int j = 0;
	int s, cs;
	for ( s = 1; s < 4; s++ )
	{
		cs = j + 1;
		while ( r - nCk(n - cs, 3 - s) > 0 )
		{
			r -= nCk(n - cs, 3 - s);
			cs += 1;
		}
		triplet[s - 1] = cs - 1;
		j = cs;
	}
}


/* Generate triplets of integers < N_reflections in random sequence */
static int **generate_triplets(int N_reflections, int N_triplets_max, int *N)
{
	int i, n, ri;
	long int N_triplets_tot = nCk(N_reflections, 3);

	int N_triplets;
	if ( N_triplets_tot > N_triplets_max || N_reflections > 1000 ) {
		N_triplets = N_triplets_max;
	} else {
		N_triplets = N_triplets_tot;
	}
	*N = N_triplets;

	int **triplets = cfmalloc(N_triplets * sizeof(int *));

	if ( triplets == NULL ) {
		ERROR("Failed to allocate triplets in generate_triplets!\n");
		return 0;
	}

	n = 0;
	if ( N_triplets_tot / N_triplets < 7 ) {
		// Reservoir sampling:
		for ( i = 0; i < N_triplets_tot; i++ ) {
			if ( n < N_triplets ) {
				triplets[n] = (int *)cfmalloc(3 * sizeof(int));
				if (triplets[n] == NULL) {
					ERROR("Failed to allocate triplet in generate_triplets!\n");
					return NULL;
				}
				get_triplet_by_index(i, N_reflections, triplets[n]);
				n += 1;
			} else {
				ri = rand() % N_triplets;
				if (ri < N_triplets)
				{
					get_triplet_by_index(i, N_reflections, triplets[ri]);
				}
			}
		}
	} else {
		// Random selection from the whole set:
		int *tidx = (int *)cfmalloc(N_triplets * sizeof(int));
		if ( tidx == NULL ) {
			ERROR("Failed to allocate tidx in generate_triplets_2!\n");
			cffree(triplets);
			return NULL;
		}
		while ( n < N_triplets ) {
			int already_in_triplets = 1;
			while ( already_in_triplets ) {
				ri = rand() % N_triplets_tot;
				already_in_triplets = 0;
				for ( i = 0; i < n; i++ ) {
					if (tidx[i] == ri) already_in_triplets = 1;
				}
			}
			tidx[n] = ri;
			triplets[n] = (int *)cfmalloc(3 * sizeof(int));
			if ( triplets[n] == NULL ) {
				ERROR("Failed to allocate triplet in generate_triplets!\n");
				return NULL;
			}
			get_triplet_by_index(ri, N_reflections, triplets[n]);
			n += 1;
		}
		cffree(tidx);
	}
	return triplets;
}


static int index_refls(gsl_vector **reflections, int N_reflections, int N_refl_max,
                       double d_max, double volume_min, double volume_max,
                       double LevelFit, double IndexFit, int N_triplets_max,
                       struct asdf_cell *c, struct fftw_vars fftw)
{

	int i, k, n;

	int N_triplets;
	int **triplets;

	if ( N_reflections < 4 ) return 0;

	gsl_vector **refl_sample;
	if ( N_reflections > N_refl_max ) {
		refl_sample = (gsl_vector **)cfmalloc(N_refl_max * sizeof(gsl_vector *));
		n = 0;
		for ( i = 0; i < N_reflections; i++ ) {
			if (i < N_refl_max) {
				refl_sample[n] = reflections[i];
				n += 1;
			} else {
				k = rand() % N_refl_max;
				if (k < N_refl_max) refl_sample[k] = reflections[i];
			}
		}
	} else {
		refl_sample = reflections;
		N_refl_max = N_reflections;
	}

	profile_start("asdf-triplets");
	triplets = generate_triplets(N_refl_max, N_triplets_max, &N_triplets);
	profile_end("asdf-triplets");

	if ( N_triplets == 0 ) return 0;

	gsl_vector *normal = gsl_vector_alloc(3);

	double projections[N_refl_max];
	double ds;

	int *fits = cfmalloc(N_refl_max * sizeof(int));
	if ( fits == NULL ) {
		ERROR("Failed to allocate fits in index_refls!\n");
		if ( N_reflections > N_refl_max ) cffree(refl_sample);
		return 0;
	}

	struct tvector *tvectors = cfmalloc(N_triplets * sizeof(struct tvector));
	if ( tvectors == NULL ) {
		ERROR("Failed to allocate tvectors in index_refls!\n");
		if ( N_reflections > N_refl_max ) cffree(refl_sample);
		cffree(fits);
		return 0;
	}

	int N_tvectors = 0;

	int n_max = 0; // maximum number of reflections fitting one of tvectors
	profile_start("asdf-search");
	for ( i = 0; i < N_triplets; i++ ) {
		if ( calc_normal(refl_sample[triplets[i][0]],
				 refl_sample[triplets[i][1]],
				 refl_sample[triplets[i][2]],
				 normal) )
		{

			/* Calculate projections of reflections to normal */
			for ( k = 0; k < N_refl_max; k++ ) {
				gsl_blas_ddot(normal, refl_sample[k], &projections[k]);
			}

			/* Find ds - period in 1d lattice of projections */
			ds = find_ds_fft(projections, N_refl_max, d_max, fftw);
			if ( ds < 0.0 ) {
				ERROR("find_ds_fft() failed.\n");
				continue;
			}

			/* Refine ds, write 1 to fits[i] if reflections[i]
			 * fits ds */
			ds = refine_ds(projections, N_refl_max, ds, LevelFit, fits);

			/* n - number of reflections fitting ds */
			n = check_refl_fitting_ds(projections, N_refl_max, ds, LevelFit);

			/* normal/ds - possible direct vector */
			gsl_vector_scale(normal, 1/ds);

			if ( n > N_refl_max / 3 && n > 6 ) {

				tvectors[N_tvectors] = tvector_new(N_refl_max);

				gsl_vector_memcpy(tvectors[N_tvectors].t, normal);
				memcpy(tvectors[N_tvectors].fits, fits,
				       N_refl_max * sizeof(int));

				tvectors[N_tvectors].n = n;

				N_tvectors++;

				if (n > n_max) n_max = n;
			}
		}

		if ( (i != 0 && i % (N_triplets_max/2) == 0) || i == N_triplets - 1 ) {
			/* Sort tvectors by length */
			qsort(tvectors, N_tvectors, sizeof(struct tvector),
			      compare_tvectors);

			/* Three shortest independent tvectors with t.n > acl
			 * determine the final cell. acl is selected for the
			 * solution with the maximum number of fitting
			 * reflections */
			profile_start("asdf-findcell");
			find_cell(tvectors, N_tvectors, IndexFit, volume_min,
				  volume_max, n_max, refl_sample, N_refl_max, c);
			profile_end("asdf-findcell");

			if ( c->n > 4 * n_max / 5 ) {
				break;
			}
		}
	}
	profile_end("asdf-search");
	cffree(fits);

	for ( i = 0; i < N_tvectors; i++ ) {
		tvector_free(tvectors[i]);
	}
	cffree(tvectors);

	for ( i = 0; i < N_triplets; i++ ) {
		cffree(triplets[i]);
	}
	cffree(triplets);

	gsl_vector_free(normal);
	if ( N_reflections > N_refl_max ) cffree(refl_sample);

	if ( c->n ) return 1;

	return 0;
}


int run_asdf(struct image *image, void *ipriv)
{
	int i, j;

	double LevelFit = 1./1000;
	double IndexFit = 1./500;
	double d_max = 1000.; // thrice the maximum expected axis length
	double volume_min = 100.;
	double volume_max = 100000000.;

	int N_refl_max;
	int N_triplets_max;

	struct asdf_private *dp = (struct asdf_private *)ipriv;

	if ( dp->fast_execution ) {
		N_refl_max = 120; // maximum number of peaks to be used for indexing
		N_triplets_max = 10000; // maximum number of triplets
	} else {
		N_refl_max = 2000;
		N_triplets_max = 20000;
	}

	if ( dp->indm & INDEXING_USE_CELL_PARAMETERS ) {

		double a, b, c, gamma, beta, alpha;
		cell_get_parameters(dp->template, &a, &b, &c,
                                                  &alpha, &beta, &gamma);

		d_max = max(a, b, c) * 3 * 1e10;

		double volume = cell_get_volume(dp->template) * 1e30;

		/* Divide volume constraints by number of lattice points per
		 * unit cell since asdf always finds primitive cell */
		int latt_points_per_uc = 1;
		char centering = cell_get_centering(dp->template);
		if ( centering == 'A' ||
		     centering == 'B' ||
		     centering == 'C' ||
		     centering == 'I' ) latt_points_per_uc = 2;
		else if ( centering == 'F' ) latt_points_per_uc = 4;

		volume_min = volume * 0.9/latt_points_per_uc;
		volume_max = volume * 1.1/latt_points_per_uc;
	}

	int n = image_feature_count(image->features);
	int N_reflections = 0;
	gsl_vector *reflections[n];

	for ( i=0; i<n; i++ ) {
		struct imagefeature *f;
		double r[3];

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		detgeom_transform_coords(&image->detgeom->panels[f->pn],
		                         f->fs, f->ss, image->lambda,
		                         0.0, 0.0, r);

		reflections[N_reflections] = gsl_vector_alloc(3);
		gsl_vector_set(reflections[N_reflections], 0, r[0]/1e10);
		gsl_vector_set(reflections[N_reflections], 1, r[1]/1e10);
		gsl_vector_set(reflections[N_reflections], 2, r[2]/1e10);
		N_reflections++;
	}

	struct asdf_cell *c = asdf_cell_new(N_reflections);
	if (c == NULL) {
		ERROR("Failed to allocate asdf_cell in run_asdf!\n");
		return 0;
	}

	if ( N_reflections == 0 ) return 0;

	j = index_refls(reflections, N_reflections, N_refl_max, d_max, volume_min,
	                volume_max, LevelFit, IndexFit, N_triplets_max, c,
	                dp->fftw);

	for ( i = 0; i < N_reflections; i++ ) {
		gsl_vector_free(reflections[i]);
	}

	if ( j ) {

		UnitCell *uc;
		Crystal *cr;
		uc = cell_new();

		cell_set_cartesian(uc, gsl_vector_get(c->axes[0], 0) * 1e-10,
				       gsl_vector_get(c->axes[0], 1) * 1e-10,
				       gsl_vector_get(c->axes[0], 2) * 1e-10,
				       gsl_vector_get(c->axes[1], 0) * 1e-10,
				       gsl_vector_get(c->axes[1], 1) * 1e-10,
				       gsl_vector_get(c->axes[1], 2) * 1e-10,
				       gsl_vector_get(c->axes[2], 0) * 1e-10,
				       gsl_vector_get(c->axes[2], 1) * 1e-10,
				       gsl_vector_get(c->axes[2], 2) * 1e-10);

		cr = crystal_new();
		if ( cr == NULL ) {
			ERROR("Failed to allocate crystal.\n");
			return 0;
		}
		crystal_set_cell(cr, uc);
		image_add_crystal(image, cr);
		asdf_cell_free(c);
		return 1;

	}

	asdf_cell_free(c);
	return 0;
}


/**
 * Prepare the ASDF indexing algorithm
 */
void *asdf_prepare(IndexingMethod *indm, UnitCell *cell, struct asdf_options *asdf_opts)
{
	struct asdf_private *dp;

	/* Flags that asdf knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_USE_CELL_PARAMETERS;

	dp = cfmalloc(sizeof(struct asdf_private));
	if ( dp == NULL ) return NULL;

	dp->template = cell;
	dp->indm = *indm;
	dp->fast_execution = asdf_opts->fast_execution;
	dp->fftw = fftw_vars_new();

	return (void *)dp;
}


void asdf_cleanup(void *pp)
{
	struct asdf_private *p;
	p = (struct asdf_private *)pp;
	fftw_vars_free(p->fftw);
	cffree(p);
}


const char *asdf_probe(UnitCell *cell)
{
	return "asdf";
}

#else /* HAVE_FFTW */

int run_asdf(struct image *image, void *ipriv)
{
       ERROR("This copy of CrystFEL was compiled without FFTW support.\n");
       return 0;
}


void *asdf_prepare(IndexingMethod *indm, UnitCell *cell)
{
       ERROR("This copy of CrystFEL was compiled without FFTW support.\n");
       ERROR("To use asdf indexing, recompile with FFTW.\n");
       return NULL;
}


const char *asdf_probe(UnitCell *cell)
{
       return NULL;
}


void asdf_cleanup(void *pp)
{
}

#endif /* HAVE_FFTW */

static void asdf_show_help()
{
	printf("Parameters for the asdf indexing algorithm:\n"
"     --asdf-fast\n"
"                           Speed up execution by limiting maximum number of peaks\n"
"                            used for indexing and the number of unit cell search\n"
"                            iterations\n"
);
}


int asdf_default_options(struct asdf_options **opts_ptr)
{
	struct asdf_options *opts;

	opts = cfmalloc(sizeof(struct asdf_options));
	if ( opts == NULL ) return ENOMEM;

	opts->fast_execution = 0;

	*opts_ptr = opts;
	return 0;
}


static error_t asdf_parse_arg(int key, char *arg,
                                  struct argp_state *state)
{
	struct asdf_options **opts_ptr = state->input;
	int r;

	switch ( key ) {

		case ARGP_KEY_INIT :
		r = asdf_default_options(opts_ptr);
		if ( r ) return r;
		break;

		case 1 :
		asdf_show_help();
		return EINVAL;

		case 2 :
		(*opts_ptr)->fast_execution = 1;
		break;

	}

	return 0;
}


static struct argp_option asdf_options[] = {

	{"help-asdf", 1, NULL, OPTION_NO_USAGE,
	 "Show options for asdf indexing algorithm", 99},

	{"asdf-fast", 2, NULL, OPTION_HIDDEN, NULL},

	{0}
};


struct argp asdf_argp = { asdf_options, asdf_parse_arg,
                              NULL, NULL, NULL, NULL, NULL };
