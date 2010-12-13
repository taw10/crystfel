/*
 * hrs-scaling.c
 *
 * Intensity scaling using generalised HRS target function
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"


/* Maximum number of iterations of NLSq scaling per macrocycle. */
#define MAX_CYCLES (10)


static void show_matrix_eqn(gsl_matrix *M, gsl_vector *v, int r)
{
	int i, j;

	for ( i=0; i<r; i++ ) {
		STATUS("[ ");
		for ( j=0; j<r; j++ ) {
			STATUS("%+9.3e ", gsl_matrix_get(M, i, j));
		}
		STATUS("][ a%2i ] = [ %+9.3e ]\n", i, gsl_vector_get(v, i));
	}
}


static double gradient(int j, signed int h, signed int k, signed int l,
                       struct image *images, int n, const char *sym)
{
	struct image *image;
	struct cpeak *spots;
	double num1 = 0.0;
	double num2 = 0.0;
	double den = 0.0;
	double num1_this = 0.0;
	double num2_this = 0.0;
	int m, i;

	/* "Everything" terms */
	for ( m=0; m<n; m++ ) {

		image = &images[m];
		spots = image->cpeaks;

		for ( i=0; i<image->n_cpeaks; i++ ) {

			signed int ha, ka, la;

			if ( !spots[i].valid ) continue;
			if ( spots[i].p < 0.1 ) continue;

			get_asymm(spots[i].h, spots[i].k, spots[i].l,
			          &ha, &ka, &la, sym);

			if ( ha != h ) continue;
			if ( ka != k ) continue;
			if ( la != l ) continue;

			num1 += pow(image->osf, 2.0) * pow(spots[i].p, 2.0);
			num2 += image->osf * spots[i].intensity * spots[i].p;
			den += pow(image->osf, 2.0) * pow(spots[i].p, 2.0);

		}

	}

	/* "This frame" terms */
	image = &images[j];
	spots = image->cpeaks;
	for ( i=0; i<image->n_cpeaks; i++ ) {

		signed int ha, ka, la;

		if ( !spots[i].valid ) continue;
		if ( spots[i].p < 0.1 ) continue;

		get_asymm(spots[i].h, spots[i].k, spots[i].l,
		          &ha, &ka, &la, sym);

		if ( ha != h ) continue;
		if ( ka != k ) continue;
		if ( la != l ) continue;

		num1_this += spots[i].intensity * spots[i].p;
		num2_this += pow(spots[i].p, 2.0);

	}

	return (num1*num1_this - num2*num2_this) / den;
}


static double iterate_scale(struct image *images, int n, const char *sym,
                            double *I_full_old)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int m;
	double max_shift;
	int crossed = 0;

	M = gsl_matrix_calloc(n-1, n-1);
	v = gsl_vector_calloc(n-1);

	for ( m=0; m<n; m++ ) {  /* "Equation number": one equation per frame */

		int k;  /* Frame (scale factor) number */
		int h;
		int mcomp;
		double vc_tot = 0.0;
		struct image *image = &images[m];
		struct cpeak *spots = image->cpeaks;

		if ( m == crossed ) continue;

		/* Determine the "solution" vector component */
		for ( h=0; h<image->n_cpeaks; h++ ) {

			double v_c;
			double ic, ip;
			signed int ha, ka, la;

			if ( !spots[h].valid ) continue;
			if ( spots[h].p < 0.1 ) continue;

			get_asymm(spots[h].h, spots[h].k, spots[h].l,
			          &ha, &ka, &la, sym);
			ic = lookup_intensity(I_full_old, ha, ka, la);

			v_c = ip - (spots[h].p * image->osf * ic);
			v_c *= pow(spots[h].p, 2.0);
			v_c *= pow(image->osf, 2.0);
			v_c *= ic;
			v_c *= gradient(m, ha, ka, la, images, n, sym);
			vc_tot += v_c;

		}

		mcomp = m;
		if ( m > crossed ) mcomp--;
		gsl_vector_set(v, mcomp, vc_tot);

		/* Now fill in the matrix components */
		for ( k=0; k<n; k++ ) {

			double val = 0.0;
			int kcomp;

			if ( k == crossed ) continue;

			for ( h=0; h<image->n_cpeaks; h++ ) {

				signed int ha, ka, la;
				double con;
				double ic;

				if ( !spots[h].valid ) continue;
				if ( spots[h].p < 0.1 ) continue;

				get_asymm(spots[h].h, spots[h].k, spots[h].l,
				          &ha, &ka, &la, sym);
				ic = lookup_intensity(I_full_old, ha, ka, la);

				con = -pow(spots[h].p, 3.0);
				con *= pow(image->osf, 3.0);
				con *= ic;
				con *= gradient(m, ha, ka, la, images, n, sym);
				con *= gradient(k, ha, ka, la, images, n, sym);
				val += con;

			}

			kcomp = k;
			if ( k > crossed ) kcomp--;
			gsl_matrix_set(M, mcomp, kcomp, val);

		}


	}
	show_matrix_eqn(M, v, n-1);

	shifts = gsl_vector_alloc(n-1);
	gsl_linalg_HH_solve(M, v, shifts);
	max_shift = 0.0;
	for ( m=0; m<n-1; m++ ) {

		double shift = gsl_vector_get(shifts, m);
		int mimg;

		mimg = m;
		if ( mimg >= crossed ) mimg++;

		images[mimg].osf += shift;

		if ( shift > max_shift ) max_shift = shift;

	}
	gsl_vector_free(shifts);

	gsl_matrix_free(M);
	gsl_vector_free(v);

	return max_shift;
}


static double *lsq_intensities(struct image *images, int n, const char *sym,
                               ReflItemList *obs)
{
	double *I_full;
	int i;

	I_full = new_list_intensity();
	for ( i=0; i<num_items(obs); i++ ) {

		signed int h, k, l;
		struct refl_item *it = get_item(obs, i);
		double num = 0.0;
		double den = 0.0;
		int m;

		get_asymm(it->h, it->k, it->l, &h, &k, &l, sym);

		/* For each frame */
		for ( m=0; m<n; m++ ) {

			double G;
			int a;

			G = images[m].osf;

			/* For each peak */
			for ( a=0; a<images[m].n_cpeaks; a++ ) {

				signed int ha, ka, la;

				if ( !images[m].cpeaks[a].valid ) continue;
				if ( images[m].cpeaks[a].p < 0.1 ) continue;

				/* Correct reflection? */
				get_asymm(images[m].cpeaks[a].h,
				          images[m].cpeaks[a].k,
				          images[m].cpeaks[a].l,
				          &ha, &ka, &la, sym);
				if ( ha != h ) continue;
				if ( ka != k ) continue;
				if ( la != l ) continue;

				num += images[m].cpeaks[a].intensity
				     * images[m].cpeaks[a].p * G;

				den += pow(images[m].cpeaks[a].p, 2.0)
				     * pow(G, 2.0);

			}

		}

		set_intensity(I_full, h, k, l, num/den);

	}

	return I_full;
}


/* Scale the stack of images */
double *scale_intensities(struct image *images, int n, const char *sym,
                          ReflItemList *obs)
{
	int m;
	double *I_full;
	int i;
	double max_shift;

	/* Start with all scale factors equal */
	for ( m=0; m<n; m++ ) {
		images[m].osf = 1.0;
	}

	/* Calculate LSQ intensities using these scale factors */
	I_full = lsq_intensities(images, n, sym, obs);

	/* Iterate */
	i = 0;
	do {

		max_shift = iterate_scale(images, n, sym, I_full);
		STATUS("Iteration %2i: max shift = %5.2f\n", i, max_shift);
		free(I_full);
		I_full = lsq_intensities(images, n, sym, obs);
		i++;

	} while ( (max_shift > 0.1) && (i < MAX_CYCLES) );

	return I_full;
}
