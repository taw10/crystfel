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
#define MAX_CYCLES (50)


static char *find_common_reflections(struct image *images, int n)
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

			if ( !get_scalable(refl) ) continue;

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
	      refl = next_found_refl(refl) )
	{
		double ic, sigi;

		if ( !get_scalable(refl) ) continue;

		ic = get_intensity(refl) / get_partiality(refl);

		/* Get the error in the estimated full intensity */
		sigi = get_esd_intensity(refl) / get_partiality(refl);

		uha_val += 1.0 / pow(sigi, 2.0);
		vha_val += ic / pow(sigi, 2.0);
	}

	*uha = uha_val;
	*vha = vha_val;
}


static void s_uhvh(struct image *images, int n,
                   signed int h, signed int k, signed int l,
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


static double iterate_scale(struct image *images, int n, RefList *scalable,
                            char *cref, RefList *reference)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	double max_shift;
	double *uha_arr;
	double *vha_arr;
	int frame;
	Reflection *refl;
	RefListIterator *iter;

	M = gsl_matrix_calloc(n, n);
	v = gsl_vector_calloc(n);

	uha_arr = malloc(n*sizeof(double));
	vha_arr = malloc(n*sizeof(double));

	for ( refl = first_refl(scalable, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		int a;
		signed int h, k, l;
		double uh, vh, Ih;

		get_indices(refl, &h, &k, &l);

		s_uhvh(images, n, h, k, l, &uh, &vh);

		if ( !reference ) {
			Ih = vh / uh;
			/* 0 / 0 = 0, not NaN */
			if ( isnan(Ih) ) Ih = 0.0;
		} else {
			/* Look up by asymmetric indices */
			Reflection *r = find_refl(reference, h, k, l);
			if ( r == NULL ) {
				ERROR("%3i %3i %3i isn't in the "
				      "reference list, so why is it "
				      "marked as scalable?\n", h, k, l);
				Ih = 0.0;
			} else {
				Ih = get_intensity(r);
			}
		}

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
			double vc, rha, vha, uha;
			double vval;

			/* Determine the "solution" vector component */
			uha = uha_arr[a];
			vha = vha_arr[a];
			rha = vha - image_a->osf * uha * Ih;
			vc = Ih * rha;
			vval = gsl_vector_get(v, a);
			gsl_vector_set(v, a, vval+vc);

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
				if ( (reference == NULL) && cref[a+n*b] != 0 ) {
					continue;
				}

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

		}
	}

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


static void lsq_intensities(struct image *images, int n, RefList *full)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double num = 0.0;
		double den = 0.0;
		int m;
		int redundancy = 0;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		/* For each frame */
		for ( m=0; m<n; m++ ) {

			double G;
			Reflection *f;

			/* Don't scale intensities from this image if
			 * post refinement failed on the last step. */
			if ( images[m].pr_dud ) continue;

			G = images[m].osf;

			for ( f = find_refl(images[m].reflections, h, k, l);
			      f != NULL;  f = next_found_refl(refl) )
			{

				double p;

				if ( !get_scalable(f) ) continue;
				p = get_partiality(f);

				num += get_intensity(f) * p * G;
				den += pow(p, 2.0) * pow(G, 2.0);
				redundancy++;

			}

		}

		set_int(refl, num/den);

		//STATUS("%3i %3i %3i %i %i\n", h, k, l,
		//        redundancy, get_redundancy(refl));

		if ( get_redundancy(refl) != redundancy ) {
			ERROR("Didn't find all the expected parts for"
			      " %3i %3i %3i\n", h, k, l);
		}

	}
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


static RefList *guess_scaled_reflections(struct image *images, int n)
{
	RefList *scalable;
	int i;

	scalable = reflist_new();
	for ( i=0; i<n; i++ ) {

		Reflection *refl;
		RefListIterator *iter;

		/* Don't scale intensities from this image if
		 * post refinement failed on the last step. */
		if ( images[i].pr_dud ) continue;

		for ( refl = first_refl(images[i].reflections, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			if ( get_scalable(refl) ) {

				signed int h, k, l;
				Reflection *f;
				int red;

				get_indices(refl, &h, &k, &l);
				f = find_refl(scalable, h, k, l);
				if ( f == NULL ) {
					f = add_refl(scalable, h, k, l);
				}

				red = get_redundancy(f);
				set_redundancy(f, red+1);

			}
		}

	}

	return scalable;
}


static UNUSED void show_scale_factors(struct image *images, int n)
{
	int i;
	for ( i=0; i<n; i++ ) {
		STATUS("Image %4i: scale factor %5.2f\n", i, images[i].osf);
	}
}


static UNUSED double total_dev(struct image *image, const RefList *full)
{
	double dev = 0.0;

	/* For each reflection */
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double G, p;
		signed int h, k, l;
		Reflection *full_version;
		double I_full, I_partial;

		if ( !get_scalable(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		assert((h!=0) || (k!=0) || (l!=0));

		full_version = find_refl(full, h, k, l);
		if ( full_version == NULL ) continue;
		/* Some reflections may have recently become scalable, but
		 * scale_intensities() might not yet have been called, so the
		 * full version may not have been calculated yet. */

		G = image->osf;
		p = get_partiality(refl);
		I_partial = get_intensity(refl);
		I_full = get_intensity(full_version);
		//STATUS("%3i %3i %3i  %5.2f  %5.2f  %5.2f  %5.2f  %5.2f\n",
		//       h, k, l, G, p, I_partial, I_full,
		//       I_partial - p*G*I_full);

		dev += pow(I_partial - p*G*I_full, 2.0);

	}

	return dev;
}


static UNUSED void plot_graph(struct image *image, RefList *reference)
{
	double sc;

	for ( sc=0.0; sc<3.0; sc+=0.1 ) {

		image->osf = sc;
		STATUS("%5.2f: %e\n", sc, total_dev(image, reference));

	}
}


/* Scale the stack of images */
RefList *scale_intensities(struct image *images, int n, RefList *reference)
{
	int m;
	RefList *full;
	int i;
	double max_shift;
	char *cref;

	/* Start with all scale factors equal */
	for ( m=0; m<n; m++ ) images[m].osf = 1.0;

	//plot_graph(images, reference);

	cref = find_common_reflections(images, n);
	full = guess_scaled_reflections(images, n);

	/* Iterate */
	i = 0;
	do {

		max_shift = iterate_scale(images, n, full, cref, reference);
		STATUS("Scaling iteration %2i: max shift = %5.2f\n",
		       i+1, max_shift);
		i++;
		if ( reference == NULL ) normalise_osfs(images, n);

	} while ( (max_shift > 0.01) && (i < MAX_CYCLES) );

	free(cref);

	//show_scale_factors(images, n);

	lsq_intensities(images, n, full);
	return full;
}
