/*
 * post-refinement.c
 *
 * Post refinement
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return ng/(R*sqrt(2.0*M_PI)) * exp(-pow(r*ng/R, 2.0)/2.0);
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
	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return -(3.0*r/(4.0*R*R)) * (1.0 - r*r/(R*R));
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
	double azi;
	double glow, ghigh;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xl, yl, zl;
	double ds;
	signed int hs, ks, ls;
	double rlow, rhigh, p;
	double philow, phihigh, phi;
	double khigh, klow;
	double tl, cet, cez;
	struct image *image = crystal_get_image(cr);
	double R = crystal_get_profile_radius(cr);
	double D, psph;

	get_partial(refl, &rlow, &rhigh, &p);

	if ( k == REF_R ) {

		double Rglow, Rghigh;

		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);

		Rglow = volume_fraction_rgradient(rlow, R, pmodel);
		Rghigh = volume_fraction_rgradient(rhigh, R, pmodel);

		return 4.0*psph/(3.0*D) + (4.0*R/(3.0*D))*(Rglow - Rghigh);

	}

	/* Calculate the gradient of partiality wrt excitation error. */
	glow = partiality_gradient(rlow, R, pmodel, rlow, rhigh);
	ghigh = partiality_gradient(rhigh, R, pmodel, rlow, rhigh);

	get_symmetric_indices(refl, &hs, &ks, &ls);
	ds = 2.0 * resolution(crystal_get_cell(cr), hs, ks, ls);

	cell_get_reciprocal(crystal_get_cell(cr), &asx, &asy, &asz,
	                                          &bsx, &bsy, &bsz,
	                                          &csx, &csy, &csz);
	xl = hs*asx + ks*bsx + ls*csx;
	yl = hs*asy + ks*bsy + ls*csy;
	zl = hs*asz + ks*bsz + ls*csz;

	/* "low" gives the largest Ewald sphere (wavelength short => k large)
	 * "high" gives the smallest Ewald sphere (wavelength long => k small)
	 */
	klow = 1.0/(image->lambda - image->lambda*image->bw/2.0);
	khigh = 1.0/(image->lambda + image->lambda*image->bw/2.0);

	tl = sqrt(xl*xl + yl*yl);

	cet = -sin(image->div/2.0) * klow;
	cez = -cos(image->div/2.0) * klow;
	philow = angle_between_2d(tl-cet, zl-cez, 0.0, 1.0);

	cet = -sin(image->div/2.0) * khigh;
	cez = -cos(image->div/2.0) * khigh;
	phihigh = angle_between_2d(tl-cet, zl-cez, 0.0, 1.0);

	/* Approximation: philow and phihigh are very similar */
	phi = (philow + phihigh) / 2.0;

	azi = atan2(yl, xl);

	/* For many gradients, just multiply the above number by the gradient
	 * of excitation error wrt whatever. */
	switch ( k ) {

		/* Cell parameters and orientation */
		case REF_ASX :
		return - hs * sin(phi) * cos(azi) * (glow-ghigh);

		case REF_BSX :
		return - ks * sin(phi) * cos(azi) * (glow-ghigh);

		case REF_CSX :
		return - ls * sin(phi) * cos(azi) * (glow-ghigh);

		case REF_ASY :
		return - hs * sin(phi) * sin(azi) * (glow-ghigh);

		case REF_BSY :
		return - ks * sin(phi) * sin(azi) * (glow-ghigh);

		case REF_CSY :
		return - ls * sin(phi) * sin(azi) * (glow-ghigh);

		case REF_ASZ :
		return - hs * cos(phi) * (glow-ghigh);

		case REF_BSZ :
		return - ks * cos(phi) * (glow-ghigh);

		case REF_CSZ :
		return - ls * cos(phi) * (glow-ghigh);

		case REF_DIV :
		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);
		return (ds/2.0)*(glow+ghigh) - 4.0*R*psph*ds/(3.0*D*D);

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
}


/* Apply the given shift to the 'k'th parameter of 'image'. */
static void apply_shift(Crystal *cr, int k, double shift)
{
	double t;
	struct image *image = crystal_get_image(cr);

	switch ( k ) {

		case REF_DIV :
		if ( isnan(shift) ) {
			ERROR("NaN divergence shift\n");
		} else {
			image->div += shift;
			if ( image->div < 0.0 ) image->div = 0.0;
		}
		break;

		case REF_R :
		t = crystal_get_profile_radius(cr);
		t += shift;
		crystal_set_profile_radius(cr, t);
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
		apply_cell_shift(crystal_get_cell(cr), k, shift);
		break;

		default :
		ERROR("No shift defined for parameter %i\n", k);
		abort();

	}
}


static int check_eigen(gsl_vector *e_val, int verbose)
{
	int i;
	double vmax, vmin;
	const int n = e_val->size;
	const double max_condition = 1e6;
	int n_filt = 0;

	if ( verbose ) STATUS("Eigenvalues:\n");
	vmin = +INFINITY;
	vmax = 0.0;
	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( verbose ) STATUS("%i: %e\n", i, val);
		if ( val > vmax ) vmax = val;
		if ( val < vmin ) vmin = val;
	}

	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( val < vmax/max_condition ) {
			gsl_vector_set(e_val, i, 0.0);
			n_filt++;
		}
	}

	vmin = +INFINITY;
	vmax = 0.0;
	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( val == 0.0 ) continue;
		if ( val > vmax ) vmax = val;
		if ( val < vmin ) vmin = val;
	}
	if ( verbose ) {
		STATUS("Condition number: %e / %e = %5.2f\n",
		       vmax, vmin, vmax/vmin);
		STATUS("%i out of %i eigenvalues filtered.\n", n_filt, n);
	}

	return n_filt;
}


static gsl_vector *solve_svd(gsl_vector *v, gsl_matrix *M, int *n_filt,
                             int verbose)
{
	gsl_matrix *s_vec;
	gsl_vector *s_val;
	int err, n;
	gsl_vector *shifts;
	gsl_vector *SB;
	gsl_vector *SinvX;
	gsl_matrix *S;  /* rescaling matrix due to Bricogne */
	gsl_matrix *AS;
	gsl_matrix *SAS;
	int i;
	gsl_matrix *SAS_copy;

	n = v->size;
	if ( v->size != M->size1 ) return NULL;
	if ( v->size != M->size2 ) return NULL;

	/* Calculate the rescaling matrix S */
	S = gsl_matrix_calloc(n, n);
	for ( i=0; i<n; i++ ) {
		double sii = pow(gsl_matrix_get(M, i, i), -0.5);
		gsl_matrix_set(S, i, i, sii);
	}

	/* Calculate the matrix SAS, which we will be (not) inverting */
	AS = gsl_matrix_calloc(n, n);
	SAS = gsl_matrix_calloc(n, n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, M, S, 0.0, AS);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S, AS, 0.0, SAS);
	gsl_matrix_free(AS);

	/* Calculate the vector SB, which is the RHS of the equation */
	SB = gsl_vector_calloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1.0, S, v, 0.0, SB);

	if ( verbose ) {
		STATUS("The equation after rescaling:\n");
		show_matrix_eqn(SAS, SB);
	}

	SAS_copy = gsl_matrix_alloc(SAS->size1, SAS->size2);
	gsl_matrix_memcpy(SAS_copy, SAS);

	/* Do the SVD */
	s_val = gsl_vector_calloc(n);
	s_vec = gsl_matrix_calloc(n, n);
	err = gsl_linalg_SV_decomp_jacobi(SAS, s_vec, s_val);
	if ( err ) {
		if ( verbose ) ERROR("SVD failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		gsl_matrix_free(SAS);
		gsl_matrix_free(S);
		return NULL;
	}
	/* "SAS" is now "U" */

	/* Filter the eigenvalues */
	*n_filt = check_eigen(s_val, verbose);

	gsl_matrix_free(SAS_copy);

	/* Solve the equation SAS.SinvX = SB */
	SinvX = gsl_vector_calloc(n);
	err = gsl_linalg_SV_solve(SAS, s_vec, s_val, SB, SinvX);
	gsl_vector_free(SB);
	gsl_matrix_free(SAS);
	gsl_matrix_free(s_vec);
	gsl_vector_free(s_val);

	if ( err ) {
		ERROR("Matrix solution failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(S);
		gsl_vector_free(SinvX);
		return NULL;
	}

	/* Calculate S.SinvX to get X, the shifts */
	shifts = gsl_vector_calloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1.0, S, SinvX, 0.0, shifts);

	gsl_matrix_free(S);
	gsl_vector_free(SinvX);

	return shifts;
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

	*n_filtered = 0;

	reflections = crystal_get_reflections(cr);

	M = gsl_matrix_calloc(NUM_PARAMS, NUM_PARAMS);
	v = gsl_vector_calloc(NUM_PARAMS);

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
		double gradients[NUM_PARAMS];

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
		for ( k=0; k<NUM_PARAMS; k++ ) {
			gradients[k] = p_gradient(cr, k, refl, pmodel) * l;
		}

		for ( k=0; k<NUM_PARAMS; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<NUM_PARAMS; g++ ) {

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

		for ( param=0; param<NUM_PARAMS; param++ ) {
			double shift = gsl_vector_get(shifts, param);
			apply_shift(cr, param, shift);
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


struct param_backup
{
	UnitCell *cell;
	double profile_radius;
	double div;
};


static struct param_backup backup_crystal(Crystal *cr)
{
	struct param_backup b;
	struct image *image = crystal_get_image(cr);

	b.cell = cell_new_from_cell(crystal_get_cell(cr));
	b.profile_radius = crystal_get_profile_radius(cr);
	b.div = image->div;

	return b;
}


static void revert_crystal(Crystal *cr, struct param_backup b)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct image *image = crystal_get_image(cr);

	cell_get_reciprocal(b.cell, &asx, &asy, &asz,
	                            &bsx, &bsy, &bsz,
	                            &csx, &csy, &csz);

	cell_set_reciprocal(crystal_get_cell(cr), asx, asy, asz,
	                                          bsx, bsy, bsz,
	                                          csx, csy, csz);

	crystal_set_profile_radius(cr, b.profile_radius);
	image->div = b.div;
}


static void free_backup_crystal(struct param_backup b)
{
	cell_free(b.cell);
}


struct prdata pr_refine(Crystal *cr, const RefList *full,
                        PartialityModel pmodel)
{
	double dev;
	int i;
	struct param_backup backup;
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

	backup = backup_crystal(cr);

	i = 0;
	do {

		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		double dev;
		int n_total;
		int n_gained = 0;
		int n_lost = 0;

		n_total = num_reflections(crystal_get_reflections(cr));
		cell_get_reciprocal(crystal_get_cell(cr), &asx, &asy, &asz,
			               &bsx, &bsy, &bsz, &csx, &csy, &csz);

		pr_iterate(cr, full, pmodel, &prdata.n_filtered);

		update_partialities_2(cr, pmodel, &n_gained, &n_lost,
		                      &mean_p_change);

		if ( verbose ) {
			dev = guide_dev(cr, full);
			STATUS("PR Iteration %2i: mean p change = %10.2f"
			       " dev = %10.5e, %i gained, %i lost, %i total\n",
			       i+1, mean_p_change, dev,
			       n_gained, n_lost, n_total);
		}

		if ( 3*n_lost > n_total ) {
			revert_crystal(cr, backup);
			update_partialities_2(cr, pmodel, &n_gained, &n_lost,
			                      &mean_p_change);
			crystal_set_user_flag(cr, 4);
			break;
		}

		i++;

	} while ( (mean_p_change > 0.01) && (i < MAX_CYCLES) );

	free_backup_crystal(backup);

	if ( crystal_get_user_flag(cr) == 0 ) {
		prdata.refined = 1;
	}

	return prdata;
}
