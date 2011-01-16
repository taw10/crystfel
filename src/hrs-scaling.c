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
#include <gsl/gsl_eigen.h>

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"


/* Maximum number of iterations of NLSq scaling per macrocycle. */
#define MAX_CYCLES (30)


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


static void show_eigen(gsl_matrix *e_vec, gsl_vector *e_val, int r)
{
	int i, j;

	for ( i=0; i<r; i++ ) {
		STATUS("[ ");
		for ( j=0; j<r; j++ ) {
			STATUS("%+5.2f ", gsl_matrix_get(e_vec, i, j));
		}
		STATUS("]   [ %+9.3e ]\n", gsl_vector_get(e_val, i));
	}
}


static void sum_GI(struct image *images, int n, const char *sym,
                   signed int hat, signed int kat, signed int lat,
                   double *sigma_GI, double *sigma_Gsq)
{
	int k;

	*sigma_GI = 0.0;
	*sigma_Gsq = 0.0;
	for ( k=0; k<n; k++ ) {

		int hi;
		struct image *image = &images[k];
		struct cpeak *spots = images->cpeaks;
		int found = 0;

		for ( hi=0; hi<image->n_cpeaks; hi++ ) {

			double ic;
			signed int ha, ka, la;

			if ( !spots[hi].valid ) continue;
			if ( spots[hi].p < 0.1 ) continue;
			get_asymm(spots[hi].h, spots[hi].k, spots[hi].l,
				  &ha, &ka, &la, sym);
			if ( ha != hat ) continue;
			if ( ka != kat ) continue;
			if ( la != lat ) continue;

			ic = spots[hi].intensity / spots[hi].p;

			*sigma_GI += ic * image->osf;
			found = 1;

		}

		if ( found ) {
			*sigma_Gsq += pow(image->osf, 2.0);
		}
	}
}


static double find_occurrances(struct image *image, const char *sym,
                               signed int h, signed int k, signed int l)
{
	double Ihl = 0.0;
	int find;
	struct cpeak *spots = image->cpeaks;

	for ( find=0; find<image->n_cpeaks; find++ ) {

		signed int ha, ka, la;

		if ( !spots[find].valid ) continue;
		if ( spots[find].p < 0.1 ) continue;
		get_asymm(spots[find].h, spots[find].k,
		          spots[find].l, &ha, &ka, &la, sym);
		if ( ha != h ) continue;
		if ( ka != k ) continue;
		if ( la != l ) continue;

		Ihl += spots[find].intensity / spots[find].p;

	}

	return Ihl;
}


static double iterate_scale(struct image *images, int n,
                            ReflItemList *obs, const char *sym)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int l;
	double max_shift;
	int n_ref;

	M = gsl_matrix_calloc(n, n);
	v = gsl_vector_calloc(n);
	n_ref = num_items(obs);

	for ( l=0; l<n; l++ ) {  /* "Equation number": one equation per frame */

		int m;  /* Frame (scale factor) number */
		int h;
		double vc_tot = 0.0;
		struct image *imagel = &images[l];

		/* Determine the "solution" vector component */
		for ( h=0; h<n_ref; h++ ) {

			double sigma_GI, sigma_Gsq;
			double vc;
			double Ihl;
			struct refl_item *it = get_item(obs, h);

			sum_GI(images, n, sym, it->h, it->k, it->l,
			       &sigma_GI, &sigma_Gsq);

			/* Add up symmetric equivalents within the pattern */
			Ihl = find_occurrances(imagel, sym,
			                       it->h, it->k, it->l);

			vc = Ihl * sigma_GI / sigma_Gsq;
			vc -= imagel->osf * pow(sigma_GI, 2.0) / sigma_Gsq;

			vc_tot += vc;

		}

		gsl_vector_set(v, l, vc_tot);

		/* Now fill in the matrix components */
		for ( m=0; m<n; m++ ) {

			double mc_tot = 0.0;
			struct image *imagem = &images[m];

			if ( m > l ) continue;  /* Matrix is symmetric */

			for ( h=0; h<n_ref; h++ ) {

				double mc = 0.0;
				double Ihl, Ihm;
				struct refl_item *it = get_item(obs, h);
				double sigma_GI, sigma_Gsq;

				sum_GI(images, n, sym, it->h, it->k, it->l,
				       &sigma_GI, &sigma_Gsq);

				if ( l == m ) {
					mc += pow(sigma_GI, 2.0)
					                  / pow(sigma_Gsq, 2.0);
				}

				Ihl = find_occurrances(imagel, sym,
			                       it->h, it->k, it->l);
			        Ihm = find_occurrances(imagem, sym,
			                       it->h, it->k, it->l);

				mc += Ihl * Ihm / sigma_Gsq;

				mc += (sigma_GI / pow(sigma_Gsq, 2.0) )
				       * ( imagel->osf*Ihm + imagem->osf * Ihl);

				mc_tot += mc;

			}

			gsl_matrix_set(M, l, m, mc_tot);
			gsl_matrix_set(M, m, l, mc_tot);

		}


	}
	show_matrix_eqn(M, v, n);

	gsl_eigen_symmv_workspace *work;
	gsl_vector *e_val;
	gsl_matrix *e_vec;
	int val;

	work = gsl_eigen_symmv_alloc(n);
	e_val = gsl_vector_alloc(n);
	e_vec = gsl_matrix_alloc(n, n);
	val = gsl_eigen_symmv(M, e_val, e_vec, work);
	STATUS("gsl_eigen_symmv said %i (%s)\n", val, gsl_strerror(val));
	gsl_eigen_symmv_free(work);

	show_eigen(e_vec, e_val, n);

#if 0  /* HRS method */
	shifts = gsl_vector_alloc(n);
	gsl_linalg_HH_solve(M, v, shifts);
	max_shift = 0.0;
	for ( l=0; l<n-1; l++ ) {

		double shift = gsl_vector_get(shifts, l);

		images[l].osf += shift;

		if ( fabs(shift) > fabs(max_shift) ) {
			max_shift = fabs(shift);
		}

	}
	gsl_vector_free(shifts);
#endif

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


static void normalise_osfs(struct image *images, int n)
{
	int m;
	double tot = 0.0;
	double mean;

	for ( m=0; m<n; m++ ) {
		tot += images[m].osf;
	}
	mean = tot / (double)n;

	for ( m=0; m<n; m++ ) {
		images[m].osf /= mean;
	}
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

		double tot = 0.0;
		int j;

		for ( j=0; j<images[m].n_cpeaks; j++ ) {
			tot += images[m].cpeaks[j].intensity
			        / images[m].cpeaks[j].p;
		}

		images[m].osf = tot;

	}
	normalise_osfs(images, n);

	/* Iterate */
	i = 0;
	do {

		max_shift = iterate_scale(images, n, obs, sym);
		STATUS("Iteration %2i: max shift = %5.2f\n", i, max_shift);
		i++;
		normalise_osfs(images, n);

	} while ( (max_shift > 0.01) && (i < MAX_CYCLES) );

	for ( m=0; m<n; m++ ) {
		images[m].osf /= images[0].osf;
	}

	I_full = lsq_intensities(images, n, obs, sym);
	return I_full;
}
