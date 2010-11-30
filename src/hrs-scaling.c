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


/* Scale the stack of images */
void scale_intensities(struct image *images, int n, const char *sym)
{
#if 0
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int j;

	M = gsl_matrix_calloc(n, n);
	v = gsl_vector_calloc(n);

	for ( j=0; j<n; j++ ) {

		signed int hind, kind, lind;
		signed int ha, ka, la;
		double I_full, delta_I;
		float I_partial;
		float xc, yc;
		int h;
		struct image *image = &images[j];
		struct cpeak *spots = image->cpeaks;

		for ( h=0; h<image->n_cpeaks; h++ ) {

			int g;
			double v_c, gr, I_full;

			hind = spots[h].h;
			kind = spots[h].k;
			lind = spots[h].l;

			/* Don't attempt to use spots with very small
			 * partialities, since it won't be accurate. */
			if ( spots[h].p < 0.1 ) continue;

			if ( integrate_peak(image, spots[h].x, spots[h].y,
				            &xc, &yc, &I_partial, NULL, NULL, 1, 1) ) {
				continue;
			}

			get_asymm(hind, kind, lind, &ha, &ka, &la, sym);
			I_full = lookup_intensity(i_full, ha, ka, la);
			delta_I = I_partial - (spots[h].p * I_full / image->osf);


			for ( g=0; g<NUM_PARAMS; g++ ) {

				double M_curr, M_c;

				M_curr = gsl_matrix_get(M, g, k);

				M_c = gradient(image, g, spots[h],
				               image->profile_radius)
				    * gradient(image, k, spots[h],
				               image->profile_radius);
				M_c *= pow(I_full, 2.0);

				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			gr = gradient(image, k, spots[h],
			              image->profile_radius);
			v_c = delta_I * I_full * gr;
			gsl_vector_set(v, k, v_c);

		}

	}
	show_matrix_eqn(M, v, NUM_PARAMS);

	shifts = gsl_vector_alloc(NUM_PARAMS);
	gsl_linalg_HH_solve(M, v, shifts);
	for ( param=0; param<NUM_PARAMS; param++ ) {
		double shift = gsl_vector_get(shifts, param);
		apply_shift(image, param, shift);
	}

	gsl_matrix_free(M);
	gsl_vector_free(v);
	gsl_vector_free(shifts);

	free(spots);
	spots = find_intersections(image, image->indexed_cell, n, 0);
	*pspots = spots;
	return mean_partial_dev(image, spots, *n, sym, i_full, NULL);
#endif
}
