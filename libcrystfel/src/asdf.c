#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <fftw3.h>

#include "image.h"
#include "dirax.h"
#include "utils.h"
#include "peaks.h"
#include "cell-utils.h"
#include "asdf.h"


struct asdf_private {
	IndexingMethod          indm;
	float                   *ltl;
	UnitCell                *template;
};


/* Possible direct vector */
struct tvector {
	gsl_vector *t;
	int n; // number of fitting reflections
	int *fits; 
};

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

struct tvector tvector_new(int n) {
	struct tvector t;
	
	t.t = gsl_vector_alloc(3);
	t.n = 0;
	t.fits = malloc(sizeof(int) * n);
	
	return t;
}

static int tvector_free(struct tvector t) {
	gsl_vector_free(t.t);
	free(t.fits);
	
	return 1;
}

static int tvector_memcpy(struct tvector *dest, struct tvector *src, int n) {
	gsl_vector_memcpy(dest->t, src->t);
	dest->n = src->n;
	memcpy(dest->fits, src->fits, sizeof(int) * n);
	
	return 1;
}

static int asdf_cell_free(struct asdf_cell c) {
	int i;
	for ( i = 0; i < 3; i++ ) {
		gsl_vector_free(c.axes[i]);
		gsl_vector_free(c.reciprocal[i]);
	}
	
	free(c.reflections);
	for ( i = 0; i < c.N_refls; i++ ) {
		free(c.indices[i]);
	}
	free(c.indices);
	
	return 1;
}

static struct asdf_cell asdf_cell_new(int n) {
	
	struct asdf_cell c;

	int i;
	for ( i = 0; i < 3; i++ ) {
		c.axes[i] = gsl_vector_alloc(3);
		c.reciprocal[i] = gsl_vector_alloc(3);
	}
	
	c.N_refls = n;
	c.reflections = malloc(sizeof(int) * n);
	
	c.indices = malloc(sizeof(int *) * n);
	for ( i = 0; i < n; i++ ) {
		c.indices[i] = malloc(sizeof(int) * 3);
	}
	
	c.n = 0;
	
	c.acl = 0;
	c.n_max = 0;
	
	return c;
}

static int asdf_cell_memcpy(struct asdf_cell *dest, struct asdf_cell *src) {
	int i;
	for ( i = 0; i < 3; i++ ) {
		gsl_vector_memcpy(dest->axes[i], src->axes[i]);
		gsl_vector_memcpy(dest->reciprocal[i], src->reciprocal[i]);
	}
	
	dest->volume = src->volume;
	
	int n = src->N_refls;
	dest->N_refls = n;
	
	dest->n = src->n;
	memcpy(dest->reflections, src->reflections, sizeof(int) * n);
	
	memcpy(dest->indices, src->indices, sizeof(int *) * n);
	
	for (i  = 0; i < n; i++ ) {
		memcpy(dest->indices[i], src->indices[i], sizeof(int) * 3);
	}
	
	dest->acl = src->acl;
	dest->n_max = src->n_max;
	return 1;
}

/* result = vec1 cross vec2 */
static int cross_product(gsl_vector *vec1, gsl_vector *vec2, 
                         gsl_vector **result) {
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

static int calc_reciprocal(gsl_vector **direct, gsl_vector **reciprocal) {
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

static int check_cell(struct asdf_private *dp, struct image *image,
                      UnitCell *cell)
{
	UnitCell *out;
	Crystal *cr;

	if ( dp->indm & INDEXING_CHECK_CELL_COMBINATIONS ) {

		out = match_cell(cell, dp->template, 0, dp->ltl, 1);
		if ( out == NULL ) return 0;

	} else if ( dp->indm & INDEXING_CHECK_CELL_AXES ) {

		out = match_cell(cell, dp->template, 0, dp->ltl, 0);
		if ( out == NULL ) return 0;

	} else {
		out = cell_new_from_cell(cell);
	}

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 0;
	}

	crystal_set_cell(cr, out);

	if ( dp->indm & INDEXING_CHECK_PEAKS ) {
		if ( !peak_sanity_check(image, &cr, 1) ) {
			crystal_free(cr);  /* Frees the cell as well */
			cell_free(out);
			return 0;
		}
	}

	image_add_crystal(image, cr);

	return 1;
}

static int compare_doubles (const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	
	return (*da > *db) - (*da < *db);
}

/* Compares tvectors by length */
static int compare_tvectors (const void *a, const void *b)
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
                       gsl_vector *normal) {
    
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
                         fftw_plan p) {    
 
	int n = N_projections; 
	double projections_sorted[n];
	memcpy(projections_sorted, projections, sizeof(double) * n);    
	qsort(projections_sorted, n, sizeof(double), compare_doubles);
	
	int i;
	
	int N = 1024; // number of points in fft calculation
	
	double in[N];
	fftw_complex *out = fftw_malloc(sizeof (fftw_complex) * N);
	
	for ( i = 0; i < N; i++ ) {
		in[i] = 0;
	}
	
	for ( i = 0; i < n; i++ ) {
		in[(int)((projections_sorted[i] - projections_sorted[0]) / 
		   (projections_sorted[n - 1] - projections_sorted[0]) * N)] ++;
	}
	
	if ( !p ) {
		p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);
	}
	
	fftw_execute_dft_r2c(p, in, out);
	
	int i_max = (int)(d_max * (projections_sorted[n - 1] - 
	                           projections_sorted[0]));
	
	int d = 1;
	double max = 0;
	double a;
	for ( i = 1; i <= i_max; i++ ) {
		a = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
		if (a > max) {
			max = a;
			d = i;
		}
	}
	
	double ds = (projections_sorted[n - 1] - projections_sorted[0]) / d;
	
	fftw_free(out);
	return ds;
}

/* Returns number of reflections fitting ds.
 * A projected reflection fits a one-dimensional lattice with elementary 
 * lattice vector d* if its absolute distance to the nearest lattice 
 * point is less than LevelFit. */
static int check_refl_fitting_ds(double *projections, int N_projections, 
                                 double ds, double LevelFit) {
	
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
                       double LevelFit, int *fits) {
	
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
                                   int N_reflections, double IndexFit) {
	      
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

static void print_asdf_cell(struct asdf_cell cc) {
	double a, b, c, alpha, beta, gamma, ab, bc, ca;
	
	a = gsl_blas_dnrm2(cc.axes[0]);
	b = gsl_blas_dnrm2(cc.axes[1]);
	c = gsl_blas_dnrm2(cc.axes[2]);
	
	gsl_blas_ddot(cc.axes[0], cc.axes[1], &ab);
	gsl_blas_ddot(cc.axes[1], cc.axes[2], &bc);
	gsl_blas_ddot(cc.axes[0], cc.axes[2], &ca);
	
	alpha = acos(bc/b/c)/M_PI*180;
	beta = acos(ca/a/c)/M_PI*180;
	gamma = acos(ab/a/b)/M_PI*180;
	
	//~ int i, j;
	//~ for (i = 0; i < 3; i ++) {
		//~ for (j = 0; j < 3; j ++) {
			//~ printf("%f ", gsl_vector_get(cc.axes[i], j));
		//~ }
		//~ printf("\n");
	//~ }
	
	printf("%.2f %.2f %.2f %.2f %.2f %.2f %.0f %d \n", a, b, c, 
						           alpha, beta, gamma, 
						           cc.volume, cc.n);

}

/* Returns 0 when refinement doesn't converge (i.e. all fitting reflections  
 * lie in the same plane) */
static int refine_asdf_cell(struct asdf_cell c, gsl_vector **reflections, 
                            int N_reflections, double IndexFit) {
    
	gsl_matrix *X = gsl_matrix_alloc(c.n, 3);
	
	gsl_vector *r[] = {gsl_vector_alloc(c.n), 
	gsl_vector_alloc(c.n), 
	gsl_vector_alloc(c.n)};
	
	gsl_vector *res = gsl_vector_alloc(3);
	gsl_matrix *cov = gsl_matrix_alloc (3, 3);
	double chisq;
	
	int i, j;
	int n = 0;
	for ( i = 0; i < N_reflections; i++ ) if ( c.reflections[i] == 1 ) 
	{
		for ( j = 0; j < 3; j++ ) {
			gsl_matrix_set(X, n, j, round(c.indices[i][j]));
			gsl_vector_set(r[j], n, 
			               gsl_vector_get(reflections[i], j));
		}
		n++;
	}
	
	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(c.n, 3);
	
	for ( i = 0; i < 3; i++ ) {
		gsl_multifit_linear (X, r[i], res, cov, &chisq, work);
		
		for (j = 0; j < 3; j++ ) {
			gsl_vector_set(c.reciprocal[j], i, 
			               gsl_vector_get(res, j));
		}
	}
	
	calc_reciprocal(c.reciprocal, c.axes);
	
	double a[3];
	for ( i = 0; i < 3; i++ ) {
		a[i] = gsl_blas_dnrm2(c.axes[i]);
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

static int reduce_asdf_cell(struct asdf_cell *cl) {
	double a, b, c, alpha, beta, gamma, ab, bc, ca, bb, cc;
   
	gsl_vector *va = gsl_vector_alloc(3);	
	gsl_vector *vb = gsl_vector_alloc(3);	
	gsl_vector *vc = gsl_vector_alloc(3);
	
	int changed = 1;
	while ( changed ) {
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
	
		if ( changed == 0 ) {

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
		}
		
		if ( changed == 0 ) {
			
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
static int find_acl(struct tvector t1, struct tvector t2, struct tvector t3) {
	int i = t1.n, j = t2.n, k = t3.n;
	if ( i <= j && i <= k ) return i;
	if ( j <= i && j <= k ) return j;
	if ( k <= i && k <= j ) return k;
}

static int find_cell(struct tvector *tvectors, int N_tvectors, double IndexFit, 
                     double volume_min, double volume_max, int n_max, 
                     gsl_vector **reflections, int N_reflections, 
                     struct asdf_cell *result) {
				  
	int i, j, k, n, cell_correct;
	
	double volume;
	
	/* Only tvectors with the number of fitting reflections > acl are 
	 * considered */
	int acl = N_reflections < 18 ? 6 : N_reflections/3;
	
	struct asdf_cell c = asdf_cell_new(N_reflections);
    
	/* Traversing a 3d array in slices perpendicular to the main diagonal */
	int sl;
	for ( sl = 0; sl < 3 * N_tvectors - 1; sl++ ) {
        
	int i_min = sl < 2 * N_tvectors ? 0 : sl - 2 * N_tvectors;
	int i_max = sl < N_tvectors ? sl : N_tvectors;
	
	for ( i = i_min; i < i_max; i++) if (tvectors[i].n > acl ) {
	
	int j_min = sl - N_tvectors - 2 * i - 1 < 0 ? i + 1 : sl - N_tvectors - i;
	int j_max = sl - N_tvectors - i < 0 ? sl - i : N_tvectors;
		
	for ( j = j_min; j < j_max; j++) if (tvectors[j].n > acl ) {
	
		k = sl - i - j - 1;
		
		if ( k > j && tvectors[k].n > acl &&
		     check_cell_angles(tvectors[i].t, 
				       tvectors[j].t, 
				       tvectors[k].t, 0.99) ) 
		{
		
			volume = calc_volume(tvectors[i].t, 
			                     tvectors[j].t, 
			                     tvectors[k].t);
			
			if ( fabs(volume) > volume_min && 
			     fabs(volume) < volume_max )
			{
				
				gsl_vector_memcpy(c.axes[0], tvectors[i].t);
				gsl_vector_memcpy(c.axes[1], tvectors[j].t);
				gsl_vector_memcpy(c.axes[2], tvectors[k].t);
				
				c.volume = volume;
				check_refl_fitting_cell(&c, reflections, 
				                        N_reflections, 
				                        IndexFit);
				
				if ( c.n < 6 ) break;
				
				reduce_asdf_cell(&c);
				
				/* if one of the cell angles > 135 or < 45 
				 * do not continue */
				if ( !check_cell_angles(c.axes[0], c.axes[1], 
				     c.axes[2], 0.71) ) break;
				
				/* index reflections with new cell axes */
				check_refl_fitting_cell(&c, reflections, 
				                        N_reflections, 
				                        IndexFit);
				
				acl = find_acl(tvectors[i], 
				               tvectors[j], 
				               tvectors[k]);
				
				c.acl = acl;
				c.n_max = n_max;
				
				/* refine cell until the number of fitting 
				 * reflections stops increasing */
				n = 0;
				cell_correct = 1;
				while ( c.n - n && cell_correct ) {
					n = c.n;
					cell_correct = refine_asdf_cell(c, 
					                          reflections, 
					                          N_reflections, 
					                          IndexFit);
					
					check_refl_fitting_cell(&c, reflections, 
					                        N_reflections, 
					                        IndexFit);
				}
				
				if ( cell_correct ) {
					reduce_asdf_cell(&c);
					if ( result->n < c.n ) {
						asdf_cell_memcpy(result, &c);      
					}
					acl++;                             
				
					if (acl > n_max) break;
					if (tvectors[j].n <= acl || 
					    tvectors[i].n <= acl) break;
				}
			}
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

static void shuffle_triplets(int **triplets, int n) {
	int i, j;
	int t[3]; 
	for ( i = 0; i < n - 1; i++ ) {
		j = i + rand() / (RAND_MAX / (n - i) + 1);
		memcpy(t, triplets[j], 3 * sizeof(int));
		memcpy(triplets[j], triplets[i], 3 * sizeof(int));
		memcpy(triplets[i], t,  3 * sizeof(int));
	}
}

static double angle_between_gsl(gsl_vector *a, gsl_vector *b) {
	double ab;
	gsl_blas_ddot(a, b, &ab);
	return acos(ab/gsl_blas_dnrm2(a)/gsl_blas_dnrm2(b)) * 180 / M_PI;
}

static int index_refls(gsl_vector **reflections, int N_reflections, 
                       double d_max, double volume_min, double volume_max,
                       double LevelFit, double IndexFit, int i_max, 
                       struct asdf_cell *c) {

	int i, j, k, l, n;
	
	/* Number of triplets = c_n^3 if n - number of reflections */
	int N_triplets = N_reflections * (N_reflections - 1) * 
	                                 (N_reflections - 2) / 6;

	int **triplets = malloc(N_triplets * sizeof(int *));
	l = 0;
	for ( i = 0; i < N_reflections; i++ ) {
		for ( j = i + 1; j < N_reflections; j++ ) {
			for ( k = j + 1; k < N_reflections; k++ ) {
				triplets[l] = malloc(3 * sizeof(int));
				
				triplets[l][0] = i;
				triplets[l][1] = j;
				triplets[l][2] = k;
				l++;
			}
		}
	}

	/* Triplets are processed in a random sequence if N_triplets > 10000 */
	if ( N_reflections > 40 ) shuffle_triplets(triplets, N_triplets);
	
	gsl_vector *normal = gsl_vector_alloc(3);
		
	double projections[N_reflections];
	double ds;
	
	int *fits = malloc(N_reflections * sizeof(int));
	
	if ( i_max > N_triplets ) i_max = N_triplets;
   
	struct tvector *tvectors = malloc(i_max * sizeof(struct tvector));
	int N_tvectors = 0;
	
	int n_max = 0; // maximum number of reflections fitting one of tvectors 
	
	fftw_plan p;
	
	for ( i = 0; i < i_max; i++ ) {
		if ( calc_normal(reflections[triplets[i][0]],
		                 reflections[triplets[i][1]],
		                 reflections[triplets[i][2]],
		                 normal) ) 
		{
			
			/* Calculate projections of reflections to normal */
			for ( k = 0; k < N_reflections; k++ ) {
				gsl_blas_ddot(normal, reflections[k], 
				              &projections[k]);
			}
			
			/* Find ds - period in 1d lattice of projections */
			ds = find_ds_fft(projections, N_reflections, d_max, p);
			
			/* Refine ds, write 1 to fits[i] if reflections[i] 
			 * fits ds */
			ds = refine_ds(projections, N_reflections, ds, LevelFit, 
			               fits);
			
			/* n - number of reflections fitting ds */
			n = check_refl_fitting_ds(projections, N_reflections, 
			                          ds, LevelFit);
			
			/* normal/ds - possible direct vector */
			gsl_vector_scale(normal, 1/ds);
			
			if ( n > N_reflections/3 && n > 6 ) {
                
				tvectors[N_tvectors] = tvector_new(N_reflections);
				
				gsl_vector_memcpy(tvectors[N_tvectors].t, 
				                  normal);
				memcpy(tvectors[N_tvectors].fits, fits, 
				       N_reflections * sizeof(int));
				
				tvectors[N_tvectors].n = n;

				N_tvectors++;
				
				if (n > n_max) n_max = n;
			}
		}
		
		if ( (i != 0 && i % 10000 == 0) || i == i_max - 1 ) {
			
			/* Sort tvectors by length */
			qsort(tvectors, N_tvectors, sizeof(struct tvector), 
			      compare_tvectors);
            
			/* Three shortest independent tvectors with t.n > acl 
			 * determine the final cell. acl is selected for the 
			 * solution with the maximum number of fitting 
			 * reflections */
		             
			find_cell(tvectors, N_tvectors, IndexFit, volume_min, 
			          volume_max, n_max, reflections, 
			          N_reflections, c);
                        
			if ( c->n > 4 * n_max / 5 ) {
				break;
			}
		}	
	}	
	free(fits);
	
	for ( i = 0; i < N_tvectors; i++ ) {
		tvector_free(tvectors[i]);
	}
	free(tvectors);
	
	for ( i = 0; i < N_triplets; i++ ) {
		free(triplets[i]);
	}
	free(triplets);

	fftw_destroy_plan(p);

	if ( c->n ) return 1;
	
	return 0;
}

double cell_get_volume(UnitCell *cell)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct rvec aCb;
	double volume;

	if ( cell_get_reciprocal(cell, &asx, &asy, &asz,
	                               &bsx, &bsy, &bsz,
	                               &csx, &csy, &csz) ) {
		ERROR("Couldn't get cell cartesian.\n");
		return 0;
	}
	
	/* "a" cross "b" */
	aCb.u = asy*bsz - asz*bsy;
	aCb.v = - (asx*bsz - asz*bsx);
	aCb.w = asx*bsy - asy*bsx;

	/* "a cross b" dot "c" */
	volume = (aCb.u*csx + aCb.v*csy + aCb.w*csz)/1e30;

	return 1/volume;
}

int run_asdf(struct image *image, IndexingPrivate *ipriv) {
	int i;
	
	double LevelFit = 1./1000;
	double IndexFit = 1./500;
	double d_max = 220.; // thrice the maximum expected axis length
	double volume_min = 100.;
	double volume_max = 1000000.;
	
	int i_max = 10000; // maximum number of triplets 
	
	struct asdf_private *dp = (struct asdf_private *)ipriv;
	
	if ( dp->indm & INDEXING_CHECK_CELL_AXES ) {
		double volume = cell_get_volume(dp->template);
		volume_min = volume * 0.95;
		volume_max = volume * 1.05;
	}
	
	int N_reflections = image_feature_count(image->features);
	gsl_vector *reflections[N_reflections];
	
	for ( i = 0; i < N_reflections; i++ ) {
		struct imagefeature *f;
	
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		
		reflections[i] = gsl_vector_alloc(3);
		gsl_vector_set(reflections[i], 0, f->rx/1e10);
		gsl_vector_set(reflections[i], 1, f->ry/1e10);
		gsl_vector_set(reflections[i], 2, f->rz/1e10);	
	}
	
	struct asdf_cell c = asdf_cell_new(N_reflections);
	
	if ( N_reflections == 0 ) return 0;
	
	i = index_refls(reflections, N_reflections, d_max, volume_min, volume_max, 
		        LevelFit, IndexFit, i_max, &c);
	
	for ( i = 0; i < N_reflections; i++ ) {
		gsl_vector_free(reflections[i]);
	}
	
	if ( i ) {
		UnitCell *uc;
		uc = cell_new();
		 
		cell_set_cartesian(uc, gsl_vector_get(c.axes[0], 0) * 1e-10, 
				       gsl_vector_get(c.axes[0], 1) * 1e-10,
				       gsl_vector_get(c.axes[0], 2) * 1e-10,
				       gsl_vector_get(c.axes[1], 0) * 1e-10,
				       gsl_vector_get(c.axes[1], 1) * 1e-10,
				       gsl_vector_get(c.axes[1], 2) * 1e-10,
				       gsl_vector_get(c.axes[2], 0) * 1e-10,
				       gsl_vector_get(c.axes[2], 1) * 1e-10,
				       gsl_vector_get(c.axes[2], 2) * 1e-10);
	
		if ( check_cell(dp, image, uc) ) {
			cell_free(uc);
			return 1;
		}
		
	cell_free(uc);
	}
	
	return 0;
}

IndexingPrivate *asdf_prepare(IndexingMethod *indm, UnitCell *cell,
                               struct detector *det, float *ltl)
{
	struct asdf_private *dp;
	int need_cell = 0;

	if ( *indm & INDEXING_CHECK_CELL_COMBINATIONS ) need_cell = 1;
	if ( *indm & INDEXING_CHECK_CELL_AXES ) need_cell = 1;

	if ( need_cell && !cell_has_parameters(cell) ) {
		ERROR("Altering your asdf flags because cell parameters were"
		      " not provided.\n");
		*indm &= ~INDEXING_CHECK_CELL_COMBINATIONS;
		*indm &= ~INDEXING_CHECK_CELL_AXES;
	}

	/* Flags that asdf knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_CELL_COMBINATIONS
	       | INDEXING_CHECK_CELL_AXES | INDEXING_CHECK_PEAKS;

	dp = malloc(sizeof(struct asdf_private));
	if ( dp == NULL ) return NULL;

	dp->ltl = ltl;
	dp->template = cell;
	dp->indm = *indm;

	return (IndexingPrivate *)dp;
}


void asdf_cleanup(IndexingPrivate *pp)
{
	struct asdf_private *p;
	p = (struct asdf_private *)pp;
	free(p);
}
