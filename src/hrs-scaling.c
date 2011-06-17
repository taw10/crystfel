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
#define MAX_CYCLES (10)


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
		sigi = ic / 10.0; /* FIXME */

		uha_val += 1.0 / pow(sigi, 2.0);
		vha_val += ic / pow(sigi, 2.0);

	}

	*uha = uha_val;
	*vha = vha_val;
}


static void s_uhvh(struct image *images, int n,
                   signed int h, signed int k, signed int l, const char *sym,
                   double *uhp, double *vhp)
{
	int a;
	double uh = 0.0;
	double vh = 0.0;

	for ( a=0; a<n; a++ ) {

		double uha, vha;

		s_uhavha(h, k, l, &images[a], &uha, &vha);

		uh += pow(images[a].osf, 2.0) * uha;
		vh += images[a].osf * vha;

	}

	*uhp = uh;
	*vhp = vh;
}


static gsl_vector *solve_by_eigenvalue_filtration(gsl_vector *v, gsl_matrix *M)
{
	gsl_eigen_symmv_workspace *work;
	gsl_vector *e_val;
	gsl_matrix *e_vec;
	int val, n, frame;
	gsl_vector *shifts;

	n = v->size;
	if ( v->size != M->size1 ) return NULL;
	if ( v->size != M->size2 ) return NULL;

	/* Diagonalise */
	work = gsl_eigen_symmv_alloc(n);
	e_val = gsl_vector_alloc(n);
	e_vec = gsl_matrix_alloc(n, n);
	val = gsl_eigen_symmv(M, e_val, e_vec, work);
	if ( val ) {
		ERROR("Couldn't diagonalise matrix.\n");
		return NULL;
	}
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

	gsl_vector_free(sprime);
	gsl_vector_free(rprime);
	gsl_matrix_free(e_vec);
	gsl_vector_free(e_val);

	return shifts;
}


static gsl_vector *solve_diagonal(gsl_vector *v, gsl_matrix *M)
{
	gsl_vector *shifts;
	int n, frame;

	n = v->size;
	if ( v->size != M->size1 ) return NULL;
	if ( v->size != M->size2 ) return NULL;

	shifts = gsl_vector_alloc(n);
	if ( shifts == NULL ) return NULL;

	for ( frame=0; frame<n; frame++ ) {

		double num, den, sh;

		num = gsl_vector_get(v, frame);
		den = gsl_matrix_get(M, frame, frame);
		sh = num/den;

		if ( isnan(sh) ) {
			gsl_vector_set(shifts, frame, 0.0);
		} else {
			gsl_vector_set(shifts, frame, sh);
		}

	}

	return shifts;
}


static double iterate_scale(struct image *images, int n,
                            ReflItemList *obs, const char *sym, char *cref,
                            double *reference)
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

		s_uhvh(images, n, it->h, it->k, it->l, sym, &uh, &vh);

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
			double vc, Ih, uh, rha, vha, uha;
			double vval;

			/* Determine the "solution" vector component */
			uha = uha_arr[a];
			vha = vha_arr[a];
			uh = lookup_intensity(uh_arr, h, k, l);

			if ( !reference ) {
				double vh;
				vh = lookup_intensity(vh_arr, h, k, l);
				Ih = vh / uh;
				/* 0 / 0 = 0, not NaN */
				if ( isnan(Ih) ) Ih = 0.0;
			} else {
				/* Look up by asymmetric indices */
				Ih = lookup_intensity(reference, h, k, l);
			}

			rha = vha - image_a->osf * uha * Ih;
			vc = Ih * rha;

			/* Determine the matrix component */
			for ( b=0; b<n; b++ ) {

				double mc = 0.0;
				double tval, rhb, vhb, uhb;
				struct image *image_b = &images[b];

				/* Matrix is symmetric */
				if ( b > a ) continue;

				/* Off-diagonals zero if reference available */
				if ( (reference != NULL) && (a!=b) ) continue;

				/* Zero if no common reflections */
				if ( cref[a+n*b] != 0 ) continue;

				uhb = uha_arr[b];
				vhb = vha_arr[b];
				rhb = vhb - image_b->osf * uhb * Ih;

				mc = (rha*vhb + vha*rhb - vha*vhb) / uh;
				if ( isnan(mc) ) mc = 0.0; /* 0 / 0 = 0 */

				if ( reference != NULL ) mc = 0.0;

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

	//show_matrix_eqn(M, v, n);

	if ( reference == NULL ) {
		shifts = solve_by_eigenvalue_filtration(v, M);
	} else {
		shifts = solve_diagonal(v, M);
	}

	gsl_matrix_free(M);
	gsl_vector_free(v);

	if ( shifts == NULL ) return 0.0;

	/* Apply shifts */
	max_shift = 0.0;
	for ( frame=0; frame<n; frame++ ) {

		double shift = gsl_vector_get(shifts, frame);

		images[frame].osf += shift;

		//STATUS("Shift %i: %5.2f: -> %5.2f\n",
		//       frame, shift, images[frame].osf);

		if ( fabs(shift) > fabs(max_shift) ) {
			max_shift = fabs(shift);
		}

	}

	gsl_vector_free(shifts);

	return max_shift;
}


static RefList *lsq_intensities(struct image *images, int n,
                                ReflItemList *obs, const char *sym)
{
	RefList *full;
	int i;

	full = reflist_new();
	for ( i=0; i<num_items(obs); i++ ) {

		struct refl_item *it = get_item(obs, i);
		double num = 0.0;
		double den = 0.0;
		int m;
		int redundancy = 0;
		Reflection *new;

		/* For each frame */
		for ( m=0; m<n; m++ ) {

			double G;
			Reflection *refl;

			G = images[m].osf;

			for ( refl = find_refl(images[m].reflections,
			                       it->h, it->k, it->l);
			      refl != NULL;
			      refl = next_found_refl(refl) )
			{

				double p;

				if ( !get_scalable(refl) ) continue;
				p = get_partiality(refl);

				num += get_intensity(refl) * p * G;
				den += pow(p, 2.0) * pow(G, 2.0);
				redundancy++;

			}

		}

		if ( !isnan(num/den) ) {
			new = add_refl(full, it->h, it->k, it->l);
			set_int(new, num/den);
			set_redundancy(new, redundancy);
		} else {
			ERROR("Couldn't calculate LSQ full intensity for"
			      "%3i %3i %3i\n", it->h, it->k, it->l);
			/* Doom is probably impending */
		}

	}

	return full;
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
RefList *scale_intensities(struct image *images, int n, const char *sym,
                           ReflItemList *obs, char *cref, double *reference)
{
	int m;
	RefList *full;
	int i;
	double max_shift;

	/* Start with all scale factors equal */
	for ( m=0; m<n; m++ ) images[m].osf = 1.0;

	/* Iterate */
	i = 0;
	do {

		max_shift = iterate_scale(images, n, obs, sym, cref, reference);
		STATUS("Scaling iteration %2i: max shift = %5.2f\n",
		       i, max_shift);
		i++;
		//normalise_osfs(images, n);

	} while ( (max_shift > 0.01) && (i < MAX_CYCLES) );

	full = lsq_intensities(images, n, obs, sym);
	return full;
}
