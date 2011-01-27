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


static double s_uha(signed int hat, signed int kat, signed int lat,
                    struct image *images, int n, const char *sym)
{
	int k;
	double val = 0.0;

	for ( k=0; k<n; k++ ) {

		int hi;
		struct image *image = &images[k];
		struct cpeak *spots = images->cpeaks;

		for ( hi=0; hi<image->n_cpeaks; hi++ ) {

			double ic, sigi;
			signed int ha, ka, la;

			if ( !spots[hi].valid ) continue;
			if ( spots[hi].p < 0.1 ) continue;
			get_asymm(spots[hi].h, spots[hi].k, spots[hi].l,
				  &ha, &ka, &la, sym);
			if ( ha != hat ) continue;
			if ( ka != kat ) continue;
			if ( la != lat ) continue;

			ic = spots[hi].intensity / spots[hi].p;
			sigi = sqrt(fabs(ic));

			val += 1.0 / pow(sigi, 2.0);

		}

	}

	return val;
}


static double s_vha(signed int hat, signed int kat, signed int lat,
                    struct image *images, int n, const char *sym)
{
	int k;
	double val = 0.0;

	for ( k=0; k<n; k++ ) {

		int hi;
		struct image *image = &images[k];
		struct cpeak *spots = images->cpeaks;

		for ( hi=0; hi<image->n_cpeaks; hi++ ) {

			double ic, sigi;
			signed int ha, ka, la;

			if ( !spots[hi].valid ) continue;
			if ( spots[hi].p < 0.1 ) continue;
			get_asymm(spots[hi].h, spots[hi].k, spots[hi].l,
				  &ha, &ka, &la, sym);
			if ( ha != hat ) continue;
			if ( ka != kat ) continue;
			if ( la != lat ) continue;

			ic = spots[hi].intensity / spots[hi].p;
			sigi = sqrt(fabs(ic));

			val += ic / pow(sigi, 2.0);  /* Yes, I know this = 1 */

		}

	}

	return val;
}


static double s_uh(struct image *images, int n, double uha)
{
	int k;
	double val = 0.0;

	for ( k=0; k<n; k++ ) {
		val += pow(images[k].osf, 2.0) * uha;
	}

	return val;
}


static double s_vh(struct image *images, int n, double vha)
{
	int k;
	double val = 0.0;

	for ( k=0; k<n; k++ ) {
		val += images[k].osf * vha;
	}

	return val;
}


static double iterate_scale(struct image *images, int n,
                            ReflItemList *obs, const char *sym)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int a;
	double max_shift;
	int crossed = 0;
	int n_ref;

	M = gsl_matrix_calloc(n-1, n-1);
	v = gsl_vector_calloc(n-1);
	n_ref = num_items(obs);

	for ( a=0; a<n; a++ ) {  /* "Equation number": one equation per frame */

		int b;  /* Frame (scale factor) number */
		int h;  /* Reflection index */
		int acomp;
		double vc_tot = 0.0;
		struct image *image_a = &images[a];

		if ( a == crossed ) continue;
		acomp = a;
		if ( a > crossed ) acomp--;

		for ( h=0; h<n_ref; h++ ) {

			double vc, Ih, uh, rha, vha, uha;
			struct refl_item *it = get_item(obs, h);

			/* Determine the "solution" vector component */
			vha = s_vha(it->h, it->k, it->l, images, n, sym);
			uha = s_uha(it->h, it->k, it->l, images, n, sym);
			uh = s_uh(images, n, uha);
			Ih = s_vh(images, n, vha) / uh;
			rha = vha - image_a->osf * uha * Ih;
			vc = Ih * rha;
			vc_tot += vc;

			/* Determine the matrix component */
			for ( b=0; b<n; b++ ) {

				double mc = 0.0;
				double tval, rhb, vhb, uhb;
				int bcomp;
				struct image *image_b = &images[b];

				/* Matrix is symmetric */
				if ( b > a ) continue;
				if ( b == crossed ) continue;
				bcomp = b;
				if ( b > crossed ) bcomp--;

				vhb = vha;  /* Because we ignore the weights */
				uhb = uha;  /* Same as above */
				rhb = vhb - image_b->osf * uhb * Ih;

				mc = (rha*vhb + vha*rhb - vha*vhb) / uh;

				if ( a == b ) {
					mc += pow(Ih, 2.0) * uha;
				}

				tval = gsl_matrix_get(M, acomp, bcomp);
				gsl_matrix_set(M, acomp, bcomp, tval+mc);
				gsl_matrix_set(M, bcomp, acomp, tval+mc);

			}

		}

		gsl_vector_set(v, acomp, vc_tot);

	}

	show_matrix_eqn(M, v, n-1);

#if 0  /* Fox and Holmes method */
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
#endif

	shifts = gsl_vector_alloc(n-1);

	gsl_linalg_HH_solve(M, v, shifts);
	max_shift = 0.0;
	for ( a=0; a<n-1; a++ ) {

		double shift = gsl_vector_get(shifts, a);
		int aimg;

		aimg = a;
		if ( aimg >= crossed ) aimg++;

		images[aimg].osf += shift;

		if ( fabs(shift) > fabs(max_shift) ) {
			max_shift = fabs(shift);
		}

	}
	gsl_vector_free(shifts);

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
