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


/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (100)


/* Refineable parameters */
enum {
	REF_SCALE,
	REF_DIV,
	REF_R,
	REF_ASX,
	REF_BSX,
	REF_CSX,
	REF_ASY,
	REF_BSY,
	REF_CSY,
	REF_ASZ,
	REF_BSZ,
	REF_CSZ,
	NUM_PARAMS
};


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
static double gradient(struct image *image, int k, Reflection *refl, double r)
{
	double ds, tt, azi;
	double nom, den;
	double g = 0.0;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xl, yl, zl;
	signed int hi, ki, li;
	double r1, r2, p;
	int clamp_low, clamp_high;

	get_indices(refl, &hi, &ki, &li);

	cell_get_reciprocal(image->indexed_cell, &asx, &asy, &asz,
	                                         &bsx, &bsy, &bsz,
	                                         &csx, &csy, &csz);
	xl = hi*asx + ki*bsx + li*csx;
	yl = hi*asy + ki*bsy + li*csy;
	zl = hi*asz + ki*bsz + li*csz;


	ds = 2.0 * resolution(image->indexed_cell, hi, ki, li);
	tt = angle_between(0.0, 0.0, 1.0,  xl, yl, zl+k);
	azi = angle_between(1.0, 0.0, 0.0, xl, yl, 0.0);

	get_partial(refl, &r1, &r2, &p, &clamp_low, &clamp_high);

	/* Calculate the gradient of partiality wrt excitation error. */
	if ( clamp_low == 0 ) {
		g += partiality_gradient(r1, r);
	}
	if ( clamp_high == 0 ) {
		g += partiality_gradient(r2, r);
	}

	/* For many gradients, just multiply the above number by the gradient
	 * of excitation error wrt whatever. */
	switch ( k ) {

	case REF_SCALE :
		return -p*pow(image->osf, -2.0);

	case REF_DIV :
		nom = sqrt(2.0) * ds * sin(image->div/2.0);
		den = sqrt(1.0 - cos(image->div/2.0));
		return (nom/den) * g;

	case REF_R :
		if ( clamp_low == 0 ) {
			g += partiality_rgradient(r1, r);
		}
		if ( clamp_high == 0 ) {
			g += partiality_rgradient(r2, r);
		}
		return g;

	/* Cell parameters and orientation */
	case REF_ASX :
		return hi * pow(sin(tt), -1.0) * cos(azi) * g;
	case REF_BSX :
		return ki * pow(sin(tt), -1.0) * cos(azi) * g;
	case REF_CSX :
		return li * pow(sin(tt), -1.0) * cos(azi) * g;
	case REF_ASY :
		return hi * pow(sin(tt), -1.0) * sin(azi) * g;
	case REF_BSY :
		return ki * pow(sin(tt), -1.0) * sin(azi) * g;
	case REF_CSY :
		return li * pow(sin(tt), -1.0) * sin(azi) * g;
	case REF_ASZ :
		return hi * pow(cos(tt), -1.0) * g;
	case REF_BSZ :
		return ki * pow(cos(tt), -1.0) * g;
	case REF_CSZ :
		return li * pow(cos(tt), -1.0) * g;

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
static void apply_shift(struct image *image, int k, double shift)
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


/* Perform one cycle of post refinement on 'image' against 'full' */
static double pr_iterate(struct image *image, const RefList *full,
                         const char *sym)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int param;
	Reflection *refl;
	RefListIterator *iter;
	RefList *reflections;
	double max_shift;

	reflections = image->reflections;

	M = gsl_matrix_calloc(NUM_PARAMS, NUM_PARAMS);
	v = gsl_vector_calloc(NUM_PARAMS);

	/* Construct the equations, one per reflection in this image */
	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int hind, kind, lind;
		signed int ha, ka, la;
		double I_full, delta_I;
		double I_partial;
		int k;
		double p;
		Reflection *match;

		get_indices(refl, &hind, &kind, &lind);

		if ( !get_scalable(refl) ) continue;

		/* Actual measurement of this reflection from this pattern? */
		I_partial = get_intensity(refl);

		get_asymm(hind, kind, lind, &ha, &ka, &la, sym);
		match = find_refl(full, ha, ka, la);
		assert(match != NULL);
		I_full = get_intensity(match);
		p = get_partiality(refl);
		delta_I = I_partial - (p * I_full / image->osf);

		for ( k=0; k<NUM_PARAMS; k++ ) {

			int g;
			double v_c, gr;

			for ( g=0; g<NUM_PARAMS; g++ ) {

				double M_curr, M_c;

				M_curr = gsl_matrix_get(M, g, k);

				M_c = gradient(image, g, refl,
				               image->profile_radius)
				    * gradient(image, k, refl,
				               image->profile_radius);
				M_c *= pow(I_full, 2.0);

				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			gr = gradient(image, k, refl, image->profile_radius);
			v_c = delta_I * I_full * gr;
			gsl_vector_set(v, k, v_c);

		}

	}
	show_matrix_eqn(M, v, NUM_PARAMS);

	shifts = gsl_vector_alloc(NUM_PARAMS);
	gsl_linalg_HH_solve(M, v, shifts);

	max_shift = 0.0;
	for ( param=0; param<NUM_PARAMS; param++ ) {
		double shift = gsl_vector_get(shifts, param);
		apply_shift(image, param, shift);
		if ( fabs(shift) > max_shift ) max_shift = fabs(shift);
	}

	gsl_matrix_free(M);
	gsl_vector_free(v);
	gsl_vector_free(shifts);

	return max_shift;
}


void pr_refine(struct image *image, const RefList *full, const char *sym)
{
	double max_shift;
	int i;

	i = 0;
	do {
		max_shift = pr_iterate(image, full, sym);
		STATUS("Iteration %2i: max shift = %5.2f\n", i, max_shift);
		i++;
	} while ( (max_shift > 0.01) && (i < MAX_CYCLES) );
}
