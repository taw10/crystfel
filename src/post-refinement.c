/*
 * post-refinement.c
 *
 * Post refinement
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
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
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "cell-utils.h"


/* Minimum partiality of a reflection for it to be used for refinement */
#define MIN_PART_REFINE (0.1)

/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (10)

/* Returns dp(gauss)/dr at "r" */
static double gaussian_fraction_gradient(double r, double R)
{
	const double ng = 2.6;
	const double sig = R/ng;

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return exp(-pow(r/sig, 2.0)/2.0) / (sig*sqrt(2.0*M_PI));
}


/* Returns dp(sph)/dr at "r" */
static double sphere_fraction_gradient(double r, double pr)
{
	double q, dpdq, dqdr;

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -pr ) return 0.0;
	if ( r > +pr ) return 0.0;

	q = (r + pr)/(2.0*pr);
	dpdq = 6.0*(q - q*q);
	dqdr = 1.0 / (2.0*pr);
	return dpdq * dqdr;
}


/* Returns dp/dr at "r" */
static double partiality_gradient(double r, double pr,
                                  PartialityModel pmodel,
                                  double rlow, double rhigh)
{
	double A, D;

	D = rlow - rhigh;

	switch ( pmodel ) {

		default:
		case PMODEL_UNITY:
		return 0.0;

		case PMODEL_SCSPHERE:
		A = sphere_fraction_gradient(r, pr)/D;
		return 4.0*pr*A/3.0;

		case PMODEL_SCGAUSSIAN:
		A = gaussian_fraction_gradient(r, pr)/D;
		return 4.0*pr*A/3.0;

	}
}


static double sphere_fraction_rgradient(double r, double R)
{
	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return -(3.0*r/(4.0*R*R)) * (1.0 - r*r/(R*R));
}


static double gaussian_fraction_rgradient(double r, double R)
{
	const double ng = 2.6;
	const double sig = R/ng;

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return -(ng*r/(sqrt(2.0*M_PI)*R*R))*exp(-r*r/(2.0*sig*sig));
}


static double volume_fraction_rgradient(double r, double pr,
                                       PartialityModel pmodel)
{
	switch ( pmodel )
	{
		case PMODEL_UNITY :
		return 1.0;

		case PMODEL_SCSPHERE :
		return sphere_fraction_rgradient(r, pr);

		case PMODEL_SCGAUSSIAN :
		return gaussian_fraction_rgradient(r, pr);
	}

	ERROR("No pmodel in volume_fraction_rgradient!\n");
	return 1.0;
}


static double volume_fraction(double rlow, double rhigh, double pr,
                              PartialityModel pmodel)
{
	switch ( pmodel )
	{
		case PMODEL_UNITY :
		return 1.0;

		case PMODEL_SCSPHERE :
		return sphere_fraction(rlow, rhigh, pr);

		case PMODEL_SCGAUSSIAN :
		return gaussian_fraction(rlow, rhigh, pr);
	}

	ERROR("No pmodel in volume_fraction!\n");
	return 1.0;
}


/* Return the gradient of partiality wrt parameter 'k' given the current status
 * of 'image'. */
double p_gradient(Crystal *cr, int k, Reflection *refl, PartialityModel pmodel)
{
	double glow, ghigh;
	double rlow, rhigh, p;
	struct image *image = crystal_get_image(cr);
	double R = crystal_get_profile_radius(cr);

	get_partial(refl, &rlow, &rhigh, &p);

	if ( k == GPARAM_R ) {

		double Rglow, Rghigh;
		double D, psph;

		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);

		Rglow = volume_fraction_rgradient(rlow, R, pmodel);
		Rghigh = volume_fraction_rgradient(rhigh, R, pmodel);

		return 4.0*psph/(3.0*D) + (4.0*R/(3.0*D))*(Rglow - Rghigh);

	}

	/* Calculate the gradient of partiality wrt excitation error. */
	glow = partiality_gradient(rlow, R, pmodel, rlow, rhigh);
	ghigh = partiality_gradient(rhigh, R, pmodel, rlow, rhigh);

	if ( k == GPARAM_DIV ) {

		double D, psph, ds;
		signed int hs, ks, ls;

		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);
		get_symmetric_indices(refl, &hs, &ks, &ls);
		ds = 2.0 * resolution(crystal_get_cell(cr), hs, ks, ls);

		return (ds/2.0)*(glow+ghigh) - 4.0*R*psph*ds/(3.0*D*D);

	}

	return r_gradient(crystal_get_cell(cr), k, refl, image) * (glow-ghigh);
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
		case GPARAM_ASX :  asx += shift;  break;
		case GPARAM_ASY :  asy += shift;  break;
		case GPARAM_ASZ :  asz += shift;  break;
		case GPARAM_BSX :  bsx += shift;  break;
		case GPARAM_BSY :  bsy += shift;  break;
		case GPARAM_BSZ :  bsz += shift;  break;
		case GPARAM_CSX :  csx += shift;  break;
		case GPARAM_CSY :  csy += shift;  break;
		case GPARAM_CSZ :  csz += shift;  break;
	}

	cell_set_reciprocal(cell, asx, asy, asz,
	                          bsx, bsy, bsz,
	                          csx, csy, csz);
}


/* Apply the given shift to the 'k'th parameter of 'image'. */
static void apply_shift(Crystal *cr, int k, double shift)
{
	double t;
	struct image *image = crystal_get_image(cr);

	switch ( k ) {

		case GPARAM_DIV :
		if ( isnan(shift) ) {
			ERROR("NaN divergence shift\n");
		} else {
			image->div += shift;
			if ( image->div < 0.0 ) image->div = 0.0;
		}
		break;

		case GPARAM_R :
		t = crystal_get_profile_radius(cr);
		t += shift;
		crystal_set_profile_radius(cr, t);
		break;

		case GPARAM_ASX :
		case GPARAM_ASY :
		case GPARAM_ASZ :
		case GPARAM_BSX :
		case GPARAM_BSY :
		case GPARAM_BSZ :
		case GPARAM_CSX :
		case GPARAM_CSY :
		case GPARAM_CSZ :
		apply_cell_shift(crystal_get_cell(cr), k, shift);
		break;

		default :
		ERROR("No shift defined for parameter %i\n", k);
		abort();

	}
}


/* Perform one cycle of post refinement on 'image' against 'full' */
static double pr_iterate(Crystal *cr, const RefList *full,
                         PartialityModel pmodel, int *n_filtered)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int param;
	Reflection *refl;
	RefListIterator *iter;
	RefList *reflections;
	double max_shift;
	int nref = 0;
	const int verbose = 0;
	int num_params = 0;
	enum gparam rv[32];

	*n_filtered = 0;

	/* If partiality model is anything other than "unity", refine all the
	 * geometrical parameters */
	if ( pmodel != PMODEL_UNITY ) {
		rv[num_params++] = GPARAM_ASX;
		rv[num_params++] = GPARAM_ASY;
		rv[num_params++] = GPARAM_ASZ;
		rv[num_params++] = GPARAM_BSX;
		rv[num_params++] = GPARAM_BSY;
		rv[num_params++] = GPARAM_BSZ;
		rv[num_params++] = GPARAM_CSX;
		rv[num_params++] = GPARAM_CSY;
		rv[num_params++] = GPARAM_CSZ;
	}

	STATUS("Refining %i parameters.\n", num_params);

	reflections = crystal_get_reflections(cr);

	M = gsl_matrix_calloc(num_params, num_params);
	v = gsl_vector_calloc(num_params);

	/* Construct the equations, one per reflection in this image */
	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int ha, ka, la;
		double I_full, delta_I;
		double w;
		double I_partial;
		int k;
		double p, l;
		Reflection *match;
		double gradients[num_params];

		/* Find the full version */
		get_indices(refl, &ha, &ka, &la);
		match = find_refl(full, ha, ka, la);
		if ( match == NULL ) continue;

		if ( (get_intensity(refl) < 3.0*get_esd_intensity(refl))
		  || (get_partiality(refl) < MIN_PART_REFINE)
		  || (get_redundancy(match) < 2) ) continue;

		I_full = get_intensity(match);

		/* Actual measurement of this reflection from this pattern? */
		I_partial = get_intensity(refl) / crystal_get_osf(cr);
		p = get_partiality(refl);
		l = get_lorentz(refl);

		/* Calculate the weight for this reflection */
		w =  pow(get_esd_intensity(refl), 2.0);
		w += l * p * I_full * pow(get_esd_intensity(match), 2.0);
		w = pow(w, -1.0);

		/* Calculate all gradients for this reflection */
		for ( k=0; k<num_params; k++ ) {
			gradients[k] = p_gradient(cr, rv[k], refl, pmodel) * l;
		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = gradients[g] * gradients[k];
				M_c *= w * pow(I_full, 2.0);

				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			delta_I = I_partial - (l * p * I_full);
			v_c = w * delta_I * I_full * gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		nref++;
	}
	if ( verbose ) {
		STATUS("The original equation:\n");
		show_matrix_eqn(M, v);
	}

	//STATUS("%i reflections went into the equations.\n", nref);
	if ( nref == 0 ) {
		crystal_set_user_flag(cr, 2);
		gsl_matrix_free(M);
		gsl_vector_free(v);
		return 0.0;
	}

	max_shift = 0.0;
	shifts = solve_svd(v, M, n_filtered, verbose);
	if ( shifts != NULL ) {

		for ( param=0; param<num_params; param++ ) {
			double shift = gsl_vector_get(shifts, param);
			apply_shift(cr, rv[param], shift);
			//STATUS("Shift %i: %e\n", param, shift);
			if ( fabs(shift) > max_shift ) max_shift = fabs(shift);
		}

	} else {
		crystal_set_user_flag(cr, 3);
	}

	gsl_matrix_free(M);
	gsl_vector_free(v);
	gsl_vector_free(shifts);

	return max_shift;
}


static double guide_dev(Crystal *cr, const RefList *full)
{
	double dev = 0.0;

	/* For each reflection */
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double G, p;
		signed int h, k, l;
		Reflection *full_version;
		double I_full, I_partial;

		if ( (get_intensity(refl) < 3.0*get_esd_intensity(refl))
		  || (get_partiality(refl) < MIN_PART_REFINE) ) continue;

		get_indices(refl, &h, &k, &l);
		assert((h!=0) || (k!=0) || (l!=0));

		full_version = find_refl(full, h, k, l);
		if ( full_version == NULL ) continue;
		/* Some reflections may have recently become scalable, but
		 * scale_intensities() might not yet have been called, so the
		 * full version may not have been calculated yet. */

		G = crystal_get_osf(cr);
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


struct prdata pr_refine(Crystal *cr, const RefList *full,
                        PartialityModel pmodel)
{
	double dev;
	int i;
	const int verbose = 0;
	struct prdata prdata;
	double mean_p_change = 0.0;

	prdata.refined = 0;
	prdata.n_filtered = 0;

	/* Don't refine crystal if scaling was bad */
	if ( crystal_get_user_flag(cr) != 0 ) return prdata;

	if ( verbose ) {
		dev = guide_dev(cr, full);
		STATUS("\n");  /* Deal with progress bar */
		STATUS("Before iteration:                       dev = %10.5e\n",
		       dev);
	}

	i = 0;
	do {

		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		double dev;

		cell_get_reciprocal(crystal_get_cell(cr), &asx, &asy, &asz,
			               &bsx, &bsy, &bsz, &csx, &csy, &csz);

		pr_iterate(cr, full, pmodel, &prdata.n_filtered);

		update_partialities(cr, pmodel);

		if ( verbose ) {
			dev = guide_dev(cr, full);
			STATUS("PR Iteration %2i: dev = %10.5e\n", i+1, dev);
		}

		i++;

	} while ( (mean_p_change > 0.01) && (i < MAX_CYCLES) );

	if ( crystal_get_user_flag(cr) == 0 ) {
		prdata.refined = 1;
	}

	return prdata;
}
