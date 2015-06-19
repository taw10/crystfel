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


/* Return the gradient of "fx" wrt parameter 'k' given the current
 * status of the crystal. */
double gradient(Crystal *cr, int k, Reflection *refl, PartialityModel pmodel)
{
	double glow, ghigh;
	double rlow, rhigh, p;
	struct image *image = crystal_get_image(cr);
	double R = crystal_get_profile_radius(cr);
	double gr;
	signed int hi, ki, li;
	double s;

	get_indices(refl, &hi, &ki, &li);
	s = resolution(crystal_get_cell(cr), hi, ki, li);
	get_partial(refl, &rlow, &rhigh, &p);

	if ( k == GPARAM_OSF ) {
		return 1.0;
	}

	if ( k == GPARAM_BFAC ) {
		return -s*s;
	}

	if ( k == GPARAM_R ) {

		double Rglow, Rghigh;
		double D, psph;

		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);

		Rglow = volume_fraction_rgradient(rlow, R, pmodel);
		Rghigh = volume_fraction_rgradient(rhigh, R, pmodel);

		gr = 4.0*psph/(3.0*D) + (4.0*R/(3.0*D))*(Rglow - Rghigh);
		return gr/p;

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

		gr = (ds/2.0)*(glow+ghigh) - 4.0*R*psph*ds/(3.0*D*D);
		return gr/p;

	}

	gr = r_gradient(crystal_get_cell(cr), k, refl, image) * (glow-ghigh);
	return gr/p;
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

	if ( isnan(shift) ) {
		ERROR("Refusing NaN shift for parameter %i\n", k);
		ERROR("Image serial %i\n", image->serial);
		return;
	}

	switch ( k ) {

		case GPARAM_DIV :
		image->div += shift;
		if ( image->div < 0.0 ) image->div = 0.0;
		break;

		case GPARAM_R :
		t = crystal_get_profile_radius(cr);
		t += shift;
		crystal_set_profile_radius(cr, t);
		break;

		case GPARAM_BFAC :
		t = crystal_get_Bfac(cr);
		t += shift;
		crystal_set_Bfac(cr, t);
		break;

		case GPARAM_OSF :
		t = crystal_get_osf(cr);
		t += shift;
		crystal_set_osf(cr, t);
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
                         PartialityModel pmodel, int no_scale, int *n_filtered,
                         int verbose)
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
	int num_params = 0;
	enum gparam rv[32];
	double G, B;

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
		//rv[num_params++] = GPARAM_R;
	}

	/* If we are scaling, refine scale factors (duh) */
	if ( !no_scale ) {
		rv[num_params++] = GPARAM_OSF;
		rv[num_params++] = GPARAM_BFAC;
	}

	if ( num_params == 0 ) return 0.0;

	reflections = crystal_get_reflections(cr);

	M = gsl_matrix_calloc(num_params, num_params);
	v = gsl_vector_calloc(num_params);

	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);

	/* Construct the equations, one per reflection in this image */
	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int ha, ka, la;
		double I_full, delta_I, esd;
		double w;
		double I_partial;
		int k;
		double p, L, s;
		double fx;
		Reflection *match;
		double gradients[num_params];

		/* If reflection is free-flagged, don't use it here */
		if ( get_flag(refl) ) continue;

		/* Find the full version */
		get_indices(refl, &ha, &ka, &la);
		match = find_refl(full, ha, ka, la);
		if ( match == NULL ) continue;

		/* Merged intensitty */
		I_full = get_intensity(match);

		/* Actual measurement of this reflection from this pattern */
		I_partial = get_intensity(refl);
		esd = get_esd_intensity(refl);

		if ( (get_partiality(refl) < MIN_PART_REFINE)
		  || (get_redundancy(match) < 2)
		  || (I_full <= 0) || (I_partial < 3*esd) ) continue;

		p = get_partiality(refl);
		L = get_lorentz(refl);
		s = resolution(crystal_get_cell(cr), ha, ka, la);

		/* Calculate the weight for this reflection */
		w = 1.0;

		/* Calculate all gradients for this reflection */
		for ( k=0; k<num_params; k++ ) {
			gradients[k] = gradient(cr, rv[k], refl, pmodel);
		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = w * gradients[g] * gradients[k];

				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			fx = G + log(p) - log(L) - B*s*s + log(I_full);
			delta_I = log(I_partial) - fx;
			v_c = w * delta_I * gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		nref++;
	}
	if ( verbose ) {
		STATUS("The original equation:\n");
		show_matrix_eqn(M, v);
		STATUS("%i reflections went into the equations.\n", nref);
	}

	if ( nref < num_params ) {
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
			if ( verbose ) STATUS("Shift %i: %e\n", param, shift);
			apply_shift(cr, rv[param], shift);
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


static double residual(Crystal *cr, const RefList *full, int verbose, int free)
{
	double dev = 0.0;
	double G, B;
	Reflection *refl;
	RefListIterator *iter;
	FILE *fh = NULL;

	if ( verbose ) {
		fh = fopen("residual.dat", "w");
	}

	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double p, L, s, w;
		signed int h, k, l;
		Reflection *match;
		double esd, I_full, I_partial;
		double fx, dc;

		if ( free && !get_flag(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;

		p = get_partiality(refl);
		L = get_lorentz(refl);
		I_partial = get_intensity(refl);
		I_full = get_intensity(match);
		esd = get_esd_intensity(refl);
		s = resolution(crystal_get_cell(cr), h, k, l);

		if ( (get_partiality(refl) < MIN_PART_REFINE)
		  || (get_redundancy(match) < 2)
		  || (I_full <= 0) || (I_partial < 3*esd) ) continue;

		fx = G + log(p) - log(L) - B*s*s + log(I_full);
		dc = log(I_partial) - fx;
		w = 1.0;
		dev += w*dc*dc;

		if ( fh != NULL ) {
			fprintf(fh, "%4i %4i %4i %e %.2f %e %f %f\n",
			        h, k, l, s, G, B, I_partial, I_full);
		}
	}

	if ( fh != NULL ) fclose(fh);

	return dev;
}


struct prdata pr_refine(Crystal *cr, const RefList *full,
		PartialityModel pmodel, int no_scale)
{
	int i;
	int verbose = 0;
	struct prdata prdata;
	int done = 0;
	double old_dev;

	prdata.refined = 0;
	prdata.n_filtered = 0;

	/* Don't refine crystal if scaling was bad */
	if ( crystal_get_user_flag(cr) != 0 ) return prdata;

	old_dev = residual(cr, full, 0, 0);
	prdata.initial_free_residual = residual(cr, full, 0, 1);
	prdata.initial_residual = old_dev;

	if ( verbose ) {
		STATUS("\n");  /* Deal with progress bar */
		STATUS("Initial G=%.2f, B=%e\n",
		       crystal_get_osf(cr), crystal_get_Bfac(cr));
		STATUS("Initial  dev =  %10.5e, free dev = %10.5e\n",
		       old_dev, residual(cr, full, 0, 1));
	}

	i = 0;
	do {

		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		double dev;

		cell_get_reciprocal(crystal_get_cell(cr), &asx, &asy, &asz,
			               &bsx, &bsy, &bsz, &csx, &csy, &csz);

		pr_iterate(cr, full, pmodel, no_scale, &prdata.n_filtered, 0);

		update_partialities(cr, pmodel);

		dev = residual(cr, full, 0, 0);
		if ( fabs(dev - old_dev) < dev*0.01 ) done = 1;

		if ( verbose ) {
			STATUS("Iter %2i: dev = %10.5e, free dev = %10.5e\n",
			       i+1, dev, residual(cr, full, 0, 1));
		}

		i++;
		old_dev = dev;

	} while ( i < MAX_CYCLES && !done );

	if ( crystal_get_user_flag(cr) == 0 ) {
		prdata.refined = 1;
	}
	prdata.final_free_residual = residual(cr, full, 0, 1);
	prdata.final_residual = old_dev;

	if ( verbose ) {
		STATUS("Final G=%.2f, B=%e\n", crystal_get_osf(cr),
		       crystal_get_Bfac(cr));
	}

	return prdata;
}
