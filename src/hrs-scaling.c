/*
 * hrs-scaling.c
 *
 * Intensity scaling using generalised HRS target function
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "utils.h"
#include "reflist.h"


/* Maximum number of iterations of NLSq scaling per macrocycle. */
#define MAX_CYCLES (30)


char *find_common_reflections(struct image *images, int n)
{
	int i;
	char *cref;

	cref = calloc(n*n, sizeof(char));

	for ( i=0; i<n; i++ ) {

		Reflection *refl;
		RefListIterator *iter;

		for ( refl = first_refl(images[i].reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			signed int h, k, l;
			int j;

			get_indices(refl, &h, &k, &l);

			for ( j=0; j<i; j++ ) {
				Reflection *r2;
				r2 = find_refl(images[j].reflections, h, k, l);
				if ( r2 != NULL ) {
					cref[i+n*j] = 1;
					cref[j+n*i] = 1;
					break;
				}
			}

		}

	}

	return cref;
}


static void s_uhavha(signed int hat, signed int kat, signed int lat,
                     struct image *image, double *uha, double *vha)
{
	double uha_val = 0.0;
	double vha_val = 0.0;
	RefList *reflections = image->reflections;
	Reflection *refl;

	for ( refl = find_refl(reflections, hat, kat, lat);
	      refl != NULL;
	      refl = next_found_refl(refl) ) {

		double ic, sigi;

		if ( !get_scalable(refl) ) continue;

		ic = get_intensity(refl) / get_partiality(refl);
		sigi = sqrt(fabs(ic));

		uha_val += 1.0 / pow(sigi, 2.0);
		vha_val += ic / pow(sigi, 2.0);

	}

	if ( uha != NULL ) *uha = uha_val;
	if ( vha != NULL ) *vha = vha_val;
}


static double s_uh(struct image *images, int n,
                   signed int h, signed int k, signed int l, const char *sym)
{
	int a;
	double val = 0.0;

	for ( a=0; a<n; a++ ) {
		double uha;
		s_uhavha(h, k, l, &images[a], &uha, NULL);
		val += pow(images[a].osf, 2.0) * uha;
	}

	return val;
}


static double s_vh(struct image *images, int n,
                   signed int h, signed int k, signed int l, const char *sym)
{
	int a;
	double val = 0.0;

	for ( a=0; a<n; a++ ) {
		double vha;
		s_uhavha(h, k, l, &images[a], NULL, &vha);
		val += images[a].osf * vha;
	}

	return val;
}


static double iterate_scale(struct image *images, int n,
                            ReflItemList *obs, const char *sym, char *cref)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	double max_shift;
	int n_ref;
	double *uh_arr;
	double *vh_arr;
	double *uha_arr;
	double *vha_arr;
	int h;  /* Reflection index */
	int frame;
	int refidx;

	M = gsl_matrix_calloc(n, n);
	v = gsl_vector_calloc(n);
	n_ref = num_items(obs);

	uh_arr = new_list_intensity();
	vh_arr = new_list_intensity();
	for ( h=0; h<n_ref; h++ ) {

		double uh, vh;
		struct refl_item *it = get_item(obs, h);

		uh = s_uh(images, n, it->h, it->k, it->l, sym);
		vh = s_vh(images, n, it->h, it->k, it->l, sym);

		set_intensity(uh_arr, it->h, it->k, it->l, uh);
		set_intensity(vh_arr, it->h, it->k, it->l, vh);

	}

	uha_arr = malloc(n*sizeof(double));
	vha_arr = malloc(n*sizeof(double));

	for ( refidx=0; refidx<n_ref; refidx++ ) {

		int a;
		struct refl_item *it = get_item(obs, refidx);
		const signed int h = it->h;
		const signed int k = it->k;
		const signed int l = it->l;

		/* For this reflection, calculate all the possible
		 * values of uha and vha */
		for ( a=0; a<n; a++ ) {

			double uha, vha;

			s_uhavha(h, k, l, &images[a], &uha, &vha);
			uha_arr[a] = uha;
			vha_arr[a] = vha;

		}

		for ( a=0; a<n; a++ ) {

			int b;  /* Frame (scale factor) number */
			struct image *image_a = &images[a];
			double vc, Ih, uh, vh, rha, vha, uha;
			double vval;

			/* Determine the "solution" vector component */
			uha = uha_arr[a];
			vha = vha_arr[a];
			uh = lookup_intensity(uh_arr, h, k, l);
			vh = lookup_intensity(vh_arr, h, k, l);
			Ih = vh / uh;
			if ( isnan(Ih) ) Ih = 0.0;  /* 0 / 0 = 0, not NaN */
			rha = vha - image_a->osf * uha * Ih;
			vc = Ih * rha;

			/* Determine the matrix component */
			for ( b=0; b<n; b++ ) {

				double mc = 0.0;
				double tval, rhb, vhb, uhb;
				struct image *image_b = &images[b];

				/* Matrix is symmetric */
				if ( b > a ) continue;

				/* Zero if no common reflections */
				if ( cref[a+n*b] != 0 ) continue;

				uhb = uha_arr[b];
				vhb = vha_arr[b];
				rhb = vhb - image_b->osf * uhb * Ih;

				mc = (rha*vhb + vha*rhb - vha*vhb) / uh;
				if ( isnan(mc) ) mc = 0.0; /* 0 / 0 = 0 */

				if ( a == b ) {
					mc += pow(Ih, 2.0) * uha;
				}

				tval = gsl_matrix_get(M, a, b);
				gsl_matrix_set(M, a, b, tval+mc);
				gsl_matrix_set(M, b, a, tval+mc);

			}

			vval = gsl_vector_get(v, a);
			gsl_vector_set(v, a, vval+vc);

		}

	}

	free(uh_arr);
	free(vh_arr);
	free(uha_arr);
	free(vha_arr);

	/* Fox and Holmes method */
	gsl_eigen_symmv_workspace *work;
	gsl_vector *e_val;
	gsl_matrix *e_vec;
	int val;

	/* Diagonalise */
	work = gsl_eigen_symmv_alloc(n);
	e_val = gsl_vector_alloc(n);
	e_vec = gsl_matrix_alloc(n, n);
	val = gsl_eigen_symmv(M, e_val, e_vec, work);
	gsl_eigen_symmv_free(work);
	val = gsl_eigen_symmv_sort(e_val, e_vec, GSL_EIGEN_SORT_ABS_DESC);

	/* Rotate the "solution vector" */
	gsl_vector *rprime;
	rprime = gsl_vector_alloc(n);
	val = gsl_blas_dgemv(CblasTrans, 1.0, e_vec, v, 0.0, rprime);

	/* Solve (now easy) */
	gsl_vector *sprime;
	sprime = gsl_vector_alloc(n);
	for ( frame=0; frame<n-1; frame++ ) {
		double num, den;
		num = gsl_vector_get(rprime, frame);
		den = gsl_vector_get(e_val, frame);
		gsl_vector_set(sprime, frame, num/den);
	}
	gsl_vector_set(sprime, n-1, 0.0);  /* Set last shift to zero */

	/* Rotate back */
	shifts = gsl_vector_alloc(n);
	val = gsl_blas_dgemv(CblasNoTrans, 1.0, e_vec, sprime, 0.0, shifts);

	/* Apply shifts */
	max_shift = 0.0;
	for ( frame=0; frame<n; frame++ ) {

		double shift = gsl_vector_get(shifts, frame);

		images[frame].osf += shift;

		if ( fabs(shift) > fabs(max_shift) ) {
			max_shift = fabs(shift);
		}

	}

	gsl_vector_free(shifts);
	gsl_vector_free(sprime);
	gsl_matrix_free(e_vec);
	gsl_vector_free(e_val);
	gsl_matrix_free(M);
	gsl_vector_free(v);

	return max_shift;
}


static double *lsq_intensities(struct image *images, int n,
                               ReflItemList *obs, const char *sym)
{
	double *I_full;
	int i;

	I_full = new_list_intensity();
	for ( i=0; i<num_items(obs); i++ ) {

		struct refl_item *it = get_item(obs, i);
		double num = 0.0;
		double den = 0.0;
		int m;

		/* For each frame */
		for ( m=0; m<n; m++ ) {

			double G;
			Reflection *refl;

			G = images[m].osf;

			for ( refl = find_refl(images[m].reflections,
			                       it->h, it->k, it->l);
			      refl != NULL;
			      refl = next_found_refl(refl) ) {

				double p;

				if ( !get_scalable(refl) ) continue;

				p = get_partiality(refl);

				num += get_intensity(refl) * p * G;
				den += pow(p, 2.0) * pow(G, 2.0);

			}

		}

		set_intensity(I_full, it->h, it->k, it->l, num/den);

	}

	return I_full;
}


static void normalise_osfs(struct image *images, int n)
{
	int m;
	double tot = 0.0;
	double norm;

	for ( m=0; m<n; m++ ) {
		tot += images[m].osf;
	}
	norm = n / tot;

	for ( m=0; m<n; m++ ) {
		images[m].osf *= norm;
	}
}


/* Scale the stack of images */
double *scale_intensities(struct image *images, int n, const char *sym,
                          ReflItemList *obs, char *cref)
{
	int m;
	double *I_full;
	int i;
	double max_shift;

	/* Start with all scale factors equal */
	for ( m=0; m<n; m++ ) images[m].osf = 1.0;

	/* Iterate */
	i = 0;
	do {

		max_shift = iterate_scale(images, n, obs, sym, cref);
		STATUS("Iteration %2i: max shift = %5.2f\n", i, max_shift);
		i++;
		normalise_osfs(images, n);

	} while ( (max_shift > 0.01) && (i < MAX_CYCLES) );

	I_full = lsq_intensities(images, n, obs, sym);
	return I_full;
}
