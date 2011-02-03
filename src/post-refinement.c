/*
 * post-refinement.c
 *
 * Post refinement
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

#include "image.h"
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"


/* Returns dp/dr at "r" */
static double partiality_gradient(double r, double profile_radius)
{
	double q, dpdq, dqdr;

	/* Calculate degree of penetration */
	q = (r + profile_radius)/(2.0*profile_radius);

	/* dp/dq */
	dpdq = 6.0*(q-pow(q, 2.0));

	/* dq/dr */
	dqdr = 1.0 / (2.0*profile_radius);

	return dpdq * dqdr;
}


/* Returns dp/drad at "r" */
static double partiality_rgradient(double r, double profile_radius)
{
	double q, dpdq, dqdrad;

	/* Calculate degree of penetration */
	q = (r + profile_radius)/(2.0*profile_radius);

	/* dp/dq */
	dpdq = 6.0*(q-pow(q, 2.0));

	/* dq/drad */
	dqdrad = 0.5 * (1.0 - r * pow(profile_radius, -2.0));

	return dpdq * dqdrad;
}


/* Return the gradient of parameter 'k' given the current status of 'image'. */
double gradient(struct image *image, int k, struct cpeak spot, double r)
{
	double ds, tt, azi;
	double nom, den;
	double g = 0.0;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xl, yl, zl;

	cell_get_reciprocal(image->indexed_cell, &asx, &asy, &asz,
	                                         &bsx, &bsy, &bsz,
	                                         &csx, &csy, &csz);
	xl = spot.h*asx + spot.k*bsx + spot.l*csx;
	yl = spot.h*asy + spot.k*bsy + spot.l*csy;
	zl = spot.h*asz + spot.k*bsz + spot.l*csz;

	ds = 2.0 * resolution(image->indexed_cell, spot.h, spot.k, spot.l);
	tt = angle_between(0.0, 0.0, 1.0,  xl, yl, zl+k);
	azi = angle_between(1.0, 0.0, 0.0, xl, yl, 0.0);

	/* Calculate the gradient of partiality wrt excitation error. */
	if ( spot.clamp1 == 0 ) {
		g += partiality_gradient(spot.r1, r);
	}
	if ( spot.clamp2 == 0 ) {
		g += partiality_gradient(spot.r2, r);
	}

	/* For many gradients, just multiply the above number by the gradient
	 * of excitation error wrt whatever. */
	switch ( k ) {

	case REF_SCALE :
		return -spot.p*pow(image->osf, -2.0);

	case REF_DIV :
		nom = sqrt(2.0) * ds * sin(image->div/2.0);
		den = sqrt(1.0 - cos(image->div/2.0));
		return (nom/den) * g;

	case REF_R :
		if ( spot.clamp1 == 0 ) {
			g += partiality_rgradient(spot.r1, r);
		}
		if ( spot.clamp2 == 0 ) {
			g += partiality_rgradient(spot.r2, r);
		}
		return g;

	/* Cell parameters and orientation */
	case REF_ASX :
		return spot.h * pow(sin(tt), -1.0) * cos(azi) * g;
	case REF_BSX :
		return spot.k * pow(sin(tt), -1.0) * cos(azi) * g;
	case REF_CSX :
		return spot.l * pow(sin(tt), -1.0) * cos(azi) * g;
	case REF_ASY :
		return spot.h * pow(sin(tt), -1.0) * sin(azi) * g;
	case REF_BSY :
		return spot.k * pow(sin(tt), -1.0) * sin(azi) * g;
	case REF_CSY :
		return spot.l * pow(sin(tt), -1.0) * sin(azi) * g;
	case REF_ASZ :
		return spot.h * pow(cos(tt), -1.0) * g;
	case REF_BSZ :
		return spot.k * pow(cos(tt), -1.0) * g;
	case REF_CSZ :
		return spot.l * pow(cos(tt), -1.0) * g;

	}

	ERROR("No gradient defined for parameter %i\n", k);
	abort();
}


static void apply_cell_shift(UnitCell *cell, int k, double shift)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	switch ( k )
	{
		case REF_ASX :  asx += shift;  break;
		case REF_ASY :  asy += shift;  break;
		case REF_ASZ :  asz += shift;  break;
		case REF_BSX :  bsx += shift;  break;
		case REF_BSY :  bsy += shift;  break;
		case REF_BSZ :  bsz += shift;  break;
		case REF_CSX :  csx += shift;  break;
		case REF_CSY :  csy += shift;  break;
		case REF_CSZ :  csz += shift;  break;
	}

	cell_set_reciprocal(cell, asx, asy, asz,
	                          bsx, bsy, bsz,
	                          csx, csy, csz);

	if ( k == REF_CSZ ) {
		double a, b, c, al, be, ga;
		cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
		STATUS("New cell: %5.2f %5.2f %5.2f nm %5.2f %5.2f %5.2f deg\n",
		       a/1.0e-9, b/1.0e-9, c/1.0e-9,
		       rad2deg(al), rad2deg(be), rad2deg(ga));
	}
}


/* Apply the given shift to the 'k'th parameter of 'image'. */
void apply_shift(struct image *image, int k, double shift)
{
	switch ( k ) {

	case REF_SCALE :
		image->osf += shift;
		STATUS("New OSF = %f (shift %e)\n", image->osf, shift);
		break;

	case REF_DIV :
		STATUS("Shifting div by %e\n", shift);
		image->div += shift;
		break;

	case REF_R :
		STATUS("Shifting r by %e\n", shift);
		image->profile_radius += shift;
		break;

	case REF_ASX :
	case REF_ASY :
	case REF_ASZ :
	case REF_BSX :
	case REF_BSY :
	case REF_BSZ :
	case REF_CSX :
	case REF_CSY :
	case REF_CSZ :
		apply_cell_shift(image->indexed_cell, k, shift);
		break;

	default :
		ERROR("No shift defined for parameter %i\n", k);
		abort();

	}
}


double mean_partial_dev(struct image *image, struct cpeak *spots, int n,
                        const char *sym, double *i_full, FILE *graph)
{
	int h, n_used;
	double delta_I = 0.0;

	n_used = 0;
	for ( h=0; h<n; h++ ) {

		signed int hind, kind, lind;
		signed int ha, ka, la;
		double I_full;
		float I_partial;
		float xc, yc;

		hind = spots[h].h;
		kind = spots[h].k;
		lind = spots[h].l;

		/* Don't attempt to use spots with very small
		 * partialities, since it won't be accurate. */
		if ( spots[h].p < 0.1 ) continue;

		/* Actual measurement of this reflection from this
		 * pattern? */
		/* FIXME: Coordinates aren't whole numbers */
		if ( integrate_peak(image, spots[h].x, spots[h].y,
		                    &xc, &yc, &I_partial, NULL, NULL,
		                    1, 0) ) {
			continue;
		}

		get_asymm(hind, kind, lind, &ha, &ka, &la, sym);
		I_full = lookup_intensity(i_full, ha, ka, la);
		delta_I += fabs(I_partial - (spots[h].p * I_full / image->osf));
		n_used++;

		if ( graph != NULL ) {
			fprintf(graph, "%3i %3i %3i %5.2f (at %5.2f,%5.2f)"
			               " %5.2f %5.2f\n",
			       hind, kind, lind, I_partial/spots[h].p, xc, yc,
			       spots[h].p, I_partial / I_full);
		}

	}

	return delta_I / (double)n_used;
}


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


/* Perform one cycle of post refinement on 'image' against 'i_full' */
double pr_iterate(struct image *image, double *i_full, const char *sym,
                  struct cpeak **pspots, int *n)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int h, param;
	struct cpeak *spots = *pspots;

	M = gsl_matrix_calloc(NUM_PARAMS, NUM_PARAMS);
	v = gsl_vector_calloc(NUM_PARAMS);

	/* Construct the equations, one per reflection in this image */
	for ( h=0; h<*n; h++ ) {

		signed int hind, kind, lind;
		signed int ha, ka, la;
		double I_full, delta_I;
		float I_partial;
		float xc, yc;
		int k;

		hind = spots[h].h;
		kind = spots[h].k;
		lind = spots[h].l;

		/* Don't attempt to use spots with very small
		 * partialities, since it won't be accurate. */
		if ( spots[h].p < 0.1 ) continue;

		/* Actual measurement of this reflection from this
		 * pattern? */
		/* FIXME: Coordinates aren't whole numbers */
		if ( integrate_peak(image, spots[h].x, spots[h].y,
		                    &xc, &yc, &I_partial, NULL, NULL,
		                    1, 0) ) {
			continue;
		}

		get_asymm(hind, kind, lind, &ha, &ka, &la, sym);
		I_full = lookup_intensity(i_full, ha, ka, la);
		delta_I = I_partial - (spots[h].p * I_full / image->osf);

		for ( k=0; k<NUM_PARAMS; k++ ) {

			int g;
			double v_c, gr;

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
}
