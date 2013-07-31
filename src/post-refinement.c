/*
 * post-refinement.c
 *
 * Post refinement
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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


/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (10)


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
	dqdrad = -0.5 * r * pow(profile_radius, -2.0);

	return dpdq * dqdrad;
}


/* Return the gradient of partiality wrt parameter 'k' given the current status
 * of 'image'. */
double p_gradient(Crystal *cr, int k, Reflection *refl, PartialityModel pmodel)
{
	double ds, azi;
	double glow, ghigh;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xl, yl, zl;
	signed int hs, ks, ls;
	double rlow, rhigh, p;
	int clamp_low, clamp_high;
	double philow, phihigh, phi;
	double khigh, klow;
	double tl, cet, cez;
	double gr;
	struct image *image = crystal_get_image(cr);
	double r = crystal_get_profile_radius(cr);

	get_symmetric_indices(refl, &hs, &ks, &ls);

	cell_get_reciprocal(crystal_get_cell(cr), &asx, &asy, &asz,
	                                          &bsx, &bsy, &bsz,
	                                          &csx, &csy, &csz);
	xl = hs*asx + ks*bsx + ls*csx;
	yl = hs*asy + ks*bsy + ls*csy;
	zl = hs*asz + ks*bsz + ls*csz;

	ds = 2.0 * resolution(crystal_get_cell(cr), hs, ks, ls);
	get_partial(refl, &rlow, &rhigh, &p, &clamp_low, &clamp_high);

	/* "low" gives the largest Ewald sphere (wavelength short => k large)
	 * "high" gives the smallest Ewald sphere (wavelength long => k small)
	 */
	klow = 1.0/(image->lambda - image->lambda*image->bw/2.0);
	khigh = 1.0/(image->lambda + image->lambda*image->bw/2.0);

	tl = sqrt(xl*xl + yl*yl);
	ds = modulus(xl, yl, zl);

	cet = -sin(image->div/2.0) * klow;
	cez = -cos(image->div/2.0) * klow;
	philow = M_PI_2 - angle_between_2d(tl-cet, zl-cez, 1.0, 0.0);

	cet = -sin(image->div/2.0) * khigh;
	cez = -cos(image->div/2.0) * khigh;
	phihigh = M_PI_2 - angle_between_2d(tl-cet, zl-cez, 1.0, 0.0);

	/* Approximation: philow and phihigh are very similar */
	phi = (philow + phihigh) / 2.0;

	azi = atan2(yl, xl);

	/* Calculate the gradient of partiality wrt excitation error. */
	if ( clamp_low == 0 ) {
		glow = partiality_gradient(rlow, r);
	} else {
		glow = 0.0;
	}
	if ( clamp_high == 0 ) {
		ghigh = partiality_gradient(rhigh, r);
	} else {
		ghigh = 0.0;
	}

	/* For many gradients, just multiply the above number by the gradient
	 * of excitation error wrt whatever. */
	switch ( k ) {

		case REF_DIV :
		/* Small angle approximation */
		return (ds*glow + ds*ghigh) / 2.0;

		case REF_R :
		gr  = partiality_rgradient(rlow, r);
		gr -= partiality_rgradient(rhigh, r);
		return gr;

		/* Cell parameters and orientation */
		case REF_ASX :
		return hs * sin(phi) * cos(azi) * (ghigh-glow);

		case REF_BSX :
		return ks * sin(phi) * cos(azi) * (ghigh-glow);

		case REF_CSX :
		return ls * sin(phi) * cos(azi) * (ghigh-glow);

		case REF_ASY :
		return hs * sin(phi) * sin(azi) * (ghigh-glow);

		case REF_BSY :
		return ks * sin(phi) * sin(azi) * (ghigh-glow);

		case REF_CSY :
		return ls * sin(phi) * sin(azi) * (ghigh-glow);

		case REF_ASZ :
		return hs * cos(phi) * (ghigh-glow);

		case REF_BSZ :
		return ks * cos(phi) * (ghigh-glow);

		case REF_CSZ :
		return ls * cos(phi) * (ghigh-glow);

	}

	ERROR("No gradient defined for parameter %i\n", k);
	abort();
}


/* Return the gradient of Lorentz factor wrt parameter 'k' given the current
 * status of 'image'. */
double l_gradient(Crystal *cr, int k, Reflection *refl, PartialityModel pmodel)
{
	double ds;
	signed int hs, ks, ls;

	switch ( k ) {

		/* Cell parameters do not affect Lorentz factor */
		case REF_ASX :
		case REF_BSX :
		case REF_CSX :
		case REF_ASY :
		case REF_BSY :
		case REF_CSY :
		case REF_ASZ :
		case REF_BSZ :
		case REF_CSZ :
		return 0.0;

		/* Nor does change of radius */
		case REF_R :
		return 0.0;

		default:
		break;

	}

	assert(k == REF_DIV);

	get_symmetric_indices(refl, &hs, &ks, &ls);

	ds = 2.0 * resolution(crystal_get_cell(cr), hs, ks, ls);

	return 2.0*crystal_get_profile_radius(cr)*ds;
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


static void check_eigen(gsl_vector *e_val)
{
	int i;
	double vmax, vmin;
	const int n = e_val->size;
	const double max_condition = 1e6;
	const int verbose = 0;

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
	}
}


static gsl_vector *solve_svd(gsl_vector *v, gsl_matrix *M)
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

	/* Do the SVD */
	s_val = gsl_vector_calloc(n);
	s_vec = gsl_matrix_calloc(n, n);
	err = gsl_linalg_SV_decomp_jacobi(SAS, s_vec, s_val);
	if ( err ) {
		ERROR("SVD failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		gsl_matrix_free(SAS);
		gsl_matrix_free(S);
		return NULL;
	}
	/* "SAS" is now "U" */

	/* Filter the eigenvalues */
	check_eigen(s_val);

	/* Calculate the vector SB, which is the RHS of the equation */
	SB = gsl_vector_calloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1.0, S, v, 0.0, SB);

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
                         PartialityModel pmodel)
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
		const double osf = crystal_get_osf(cr);

		if ( !get_refinable(refl) ) continue;

		/* Find the full version */
		get_indices(refl, &ha, &ka, &la);
		match = find_refl(full, ha, ka, la);
		if ( match == NULL ) {
			ERROR("%3i %3i %3i isn't in the reference list, so why "
			      " is it marked as refinable?\n", ha, ka, la);
			continue;
		}
		I_full = get_intensity(match);

		/* Actual measurement of this reflection from this pattern? */
		I_partial = osf * get_intensity(refl);
		p = get_partiality(refl);
		l = get_lorentz(refl);

		/* Calculate the weight for this reflection */
		w =  pow(get_esd_intensity(refl), 2.0);
		w += l * p * I_full * pow(get_esd_intensity(match), 2.0);
		w = pow(w, -1.0);

		/* Calculate all gradients for this reflection */
		for ( k=0; k<NUM_PARAMS; k++ ) {
			double gr;
			gr  = p_gradient(cr, k, refl, pmodel) * l;
			gr += l_gradient(cr, k, refl, pmodel) * p;
			gradients[k] = gr;
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
	//show_matrix_eqn(M, v, NUM_PARAMS);

	//STATUS("%i reflections went into the equations.\n", nref);
	if ( nref == 0 ) {
		crystal_set_user_flag(cr, 1);
		gsl_matrix_free(M);
		gsl_vector_free(v);
		return 0.0;
	}

	max_shift = 0.0;
	shifts = solve_svd(v, M);
	if ( shifts != NULL ) {

		for ( param=0; param<NUM_PARAMS; param++ ) {
			double shift = gsl_vector_get(shifts, param);
			apply_shift(cr, param, shift);
			//STATUS("Shift %i: %e\n", param, shift);
			if ( fabs(shift) > max_shift ) max_shift = fabs(shift);
		}

	} else {
		ERROR("Problem solving equations.\n");
		/* Leave things as they were */
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

		if ( !get_refinable(refl) ) continue;

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


static Crystal *backup_crystal(Crystal *cr)
{
	Crystal *b;

	b = crystal_new();
	crystal_set_cell(b, cell_new_from_cell(crystal_get_cell(cr)));

	return b;
}


static void revert_crystal(Crystal *cr, Crystal *backup)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(crystal_get_cell(backup),
	                    &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz,
	                    &csx, &csy, &csz);

	cell_set_reciprocal(crystal_get_cell(cr),
	                    asx, asy, asz,
	                    bsx, bsy, bsz,
	                    csx, csy, csz);
}


static void free_backup_crystal(Crystal *cr)
{
	cell_free(crystal_get_cell(cr));
	crystal_free(cr);
}


void pr_refine(Crystal *cr, const RefList *full, PartialityModel pmodel)
{
	double max_shift, dev;
	int i;
	Crystal *backup;
	const int verbose = 0;

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

		max_shift = pr_iterate(cr, full, pmodel);

		update_partialities_2(cr, pmodel, &n_gained, &n_lost);

		if ( verbose ) {
			dev = guide_dev(cr, full);
			STATUS("PR Iteration %2i: max shift = %10.2f"
			       " dev = %10.5e, %i gained, %i lost, %i total\n",
			       i+1, max_shift, dev, n_gained, n_lost, n_total);
		}

		if ( 3*n_lost > n_total ) {
			revert_crystal(cr, backup);
			crystal_set_user_flag(cr, 1);
		}

		i++;

	} while ( (max_shift > 50.0) && (i < MAX_CYCLES) );

	free_backup_crystal(backup);
}
