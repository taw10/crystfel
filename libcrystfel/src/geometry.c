/*
 * geometry.c
 *
 * Geometry of diffraction
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <assert.h>
#include <fenv.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>
#include <ctype.h>
#include <math.h>

#include "utils.h"
#include "cell.h"
#include "cell-utils.h"
#include "image.h"
#include "peaks.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "geometry.h"

/** \file geometry.h */

static int locate_peak_on_panel(double x, double y, double z, double k,
                                struct detgeom_panel *p,
                                double det_shift_x, double det_shift_y,
                                double *pfs, double *pss)
{
	gsl_vector *v;
	gsl_vector *t;
	gsl_matrix *M;
	double fs, ss, one_over_mu;

	/* Set up matrix equation */
	M = gsl_matrix_alloc(3, 3);
	v = gsl_vector_alloc(3);
	t = gsl_vector_alloc(3);
	if ( (M==NULL) || (v==NULL) || (t==NULL) ) {
		ERROR("Failed to allocate vectors for prediction\n");
		return 0;
	}

	/* Scattered ray vector, whose direction is the same in reciproal and
	 * real space.   In ray_vector_gradient(), we will calculate the
	 * gradient of this vector w.r.t. everything. */
	gsl_vector_set(t, 0, x);
	gsl_vector_set(t, 1, y);
	gsl_vector_set(t, 2, k+z);

	gsl_matrix_set(M, 0, 0, p->cnx+(det_shift_x/p->pixel_pitch));
	gsl_matrix_set(M, 0, 1, p->fsx);
	gsl_matrix_set(M, 0, 2, p->ssx);
	gsl_matrix_set(M, 1, 0, p->cny+(det_shift_y/p->pixel_pitch));
	gsl_matrix_set(M, 1, 1, p->fsy);
	gsl_matrix_set(M, 1, 2, p->ssy);
	gsl_matrix_set(M, 2, 0, p->cnz);
	gsl_matrix_set(M, 2, 1, p->fsz);
	gsl_matrix_set(M, 2, 2, p->ssz);

	if ( gsl_linalg_HH_solve(M, t, v) ) {
		ERROR("Failed to solve prediction equation\n");
		return 0;
	}

	one_over_mu = gsl_vector_get(v, 0);
	fs = gsl_vector_get(v, 1) / one_over_mu;
	ss = gsl_vector_get(v, 2) / one_over_mu;
	gsl_vector_free(v);
	gsl_vector_free(t);
	gsl_matrix_free(M);

	*pfs = fs;  *pss = ss;

	/* If "mu" is negative, then the reflection is in the
	 * wrong direction */
	if ( one_over_mu < 0.0 ) return 0;

	/* Now, is this on this panel? */
	if ( fs < 0.0 ) return 0;
	if ( fs >= p->w ) return 0;
	if ( ss < 0.0 ) return 0;
	if ( ss >= p->h ) return 0;

	return 1;
}

static signed int locate_peak(double x, double y, double z, double k,
                              struct detgeom *det,
                              double det_shift_x, double det_shift_y,
                              double *pfs, double *pss)
{
	int i;

	*pfs = -1;  *pss = -1;

	for ( i=0; i<det->n_panels; i++ ) {

		struct detgeom_panel *p;

		p = &det->panels[i];

		if ( locate_peak_on_panel(x, y, z, k, p,
			                  det_shift_x, det_shift_y,
			                  pfs, pss) ) return i; /* Woohoo! */

	}

	return -1;
}


double sphere_fraction(double rlow, double rhigh, double pr)
{
	double qlow, qhigh;
	double plow, phigh;

	/* If the "lower" Ewald sphere is a long way away, use the
	 * position at which the Ewald sphere would just touch the
	 * reflection.
	 *
	 * The six possible combinations of clamp_{low,high} (including
	 * zero) correspond to the six situations in Table 3 of Rossmann
	 * et al. (1979).
	 */
	if ( rlow < -pr ) rlow = -pr;
	if ( rlow > +pr ) rlow = +pr;
	if ( rhigh < -pr ) rhigh = -pr;
	if ( rhigh > +pr ) rhigh = +pr;

	/* Calculate degrees of penetration */
	qlow  = (rlow + pr)/(2.0*pr);
	qhigh = (rhigh + pr)/(2.0*pr);

	plow  = 3.0*qlow*qlow - 2.0*qlow*qlow*qlow;
	phigh = 3.0*qhigh*qhigh - 2.0*qhigh*qhigh*qhigh;

	return plow - phigh;
}


double gaussian_fraction(double rlow, double rhigh, double R)
{
	double plow, phigh;
	const double ng = 2.6;
	const double sig = R/ng;

	/* If the "lower" Ewald sphere is a long way away, use the
	 * position at which the Ewald sphere would just touch the
	 * reflection.
	 *
	 * The six possible combinations of clamp_{low,high} (including
	 * zero) correspond to the six situations in Table 3 of Rossmann
	 * et al. (1979).
	 */
	if ( rlow < -R ) rlow = -R;
	if ( rlow > +R ) rlow = +R;
	if ( rhigh < -R ) rhigh = -R;
	if ( rhigh > +R ) rhigh = +R;

	plow =  0.5*(1.0 + gsl_sf_erf(rlow/(sig*sqrt(2.0))));
	phigh =  0.5*(1.0 + gsl_sf_erf(rhigh/(sig*sqrt(2.0))));

	return plow - phigh;
}


static void reseed(gsl_rng *rng)
{
	unsigned long int seed = gsl_rng_get(rng);
	gsl_rng_set(rng, seed);
}


static double random_partiality(signed int h, signed int k, signed int l,
                                int serial)
{
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	double p;
	int i;

	gsl_rng_set(rng, serial);
	reseed(rng);

	for ( i=0; i<abs(h)+1; i++ ) {
		gsl_rng_get(rng);
	}
	reseed(rng);
	if ( h >= 0 ) {
		reseed(rng);
	}

	for ( i=0; i<abs(k)+1; i++ ) {
		gsl_rng_get(rng);
	}
	reseed(rng);
	if ( k >= 0 ) {
		reseed(rng);
	}

	for ( i=0; i<abs(l)+1; i++ ) {
		gsl_rng_get(rng);
	}
	reseed(rng);
	if ( l >= 0 ) {
		reseed(rng);
	}

	p = gsl_rng_uniform(rng);
	gsl_rng_free(rng);
	return p;
}


static inline double safe_khalf(const double xl,
                                const double yl,
                                const double zl)
{
	if ( zl > 0.0 ) return NAN;
	return -(xl*xl+yl*yl+zl*zl) / (2.0*zl);
}


static Reflection *check_reflection(struct image *image, Crystal *cryst,
                                    signed int h, signed int k, signed int l,
                                    double xl, double yl, double zl,
                                    Reflection *updateme)
{
	Reflection *refl;
	double R;
	int i, n;
	double knom, khalf;
	double dcs, exerr;
	double partiality, mean_kpred, M2_kpred;
	double sumw_k, mean_k, M2_k;

	/* This arbitrary value is there to mimic previous behaviour */
	const double min_partiality = exp(-0.5*1.7*1.7);

	/* Don't predict 000 */
	if ( (updateme == NULL) && (abs(h)+abs(k)+abs(l) == 0) ) return NULL;

	R = fabs(crystal_get_profile_radius(cryst));
	n = spectrum_get_num_gaussians(image->spectrum);
	assert(n > 0);

	partiality = 0.0;
	mean_kpred = 0.0;
	M2_kpred = 0.0;
	sumw_k = 0.0;
	mean_k = 0.0;
	M2_k = 0.0;
	for ( i=0; i<n; i++ ) {

		struct gaussian g;
		double kpred;
		double exerr2, x, y, z, norm;
		double sigma_proj, w0, w1;
		double sigma2, exponent, overlap_integral;

		g = spectrum_get_gaussian(image->spectrum, i);

		/* Project lattice point onto Ewald sphere */
		x = xl;
		y = yl;
		z = zl + g.kcen;
		norm = 1.0/sqrt(x*x+y*y+z*z);
		x *= norm;
		y *= norm;
		z *= norm;

		/* Width of Ewald sphere in the direction of the projection */
		sigma_proj = (1-z)*g.sigma;

		mean_variance(g.kcen, g.area, &sumw_k, &mean_k, &M2_k);
		M2_k += g.area * g.sigma * g.sigma;
		w0 = 1.0/(R*R);
		w1 = 1.0/(sigma_proj*sigma_proj);

		x *= g.kcen;
		y *= g.kcen;
		z *= g.kcen;
		z -= g.kcen;

		/* Three because the general case fails in extreme cases */
		if ( isinf(w0) || (w0 / w1 <= DBL_MIN) ) {

			/* 'Laue' corner case */
			kpred = g.kcen;
			exerr2 = g.kcen - safe_khalf(xl, yl, zl);
			exerr2 *= exerr2;

		} else if ( w1 / w0 <= DBL_MIN ) {

			/* 'Monochromatic' corner case */
			kpred = safe_khalf(xl,yl,zl);
			exerr2 = g.kcen - kpred;
			exerr2*= exerr2;

		} else {

			/* General case */

			/* Closest point on Ewald sphere.
			 * Project zl to 0, bit of a hack... */
			const double zlp0 = zl<0?zl:0;
			exerr2 = (x-xl)*(x-xl) + (y-yl)*(y-yl) + (z-zl)*(z-zl);

			/* Weighted average between projected lattice point
			 * and Ewald sphere */
			x = ( xl  *w0 + x*w1 ) / ( w0 + w1 );
			y = ( yl  *w0 + y*w1 ) / ( w0 + w1 );
			z = ( zlp0*w0 + z*w1 ) / ( w0 + w1 );
			kpred = safe_khalf(x,y,z);

		}
		sigma2 = R*R + sigma_proj*sigma_proj;
		exponent = - 0.5 * exerr2 / sigma2;
		if ( exponent > -700.0 ) {
			overlap_integral = exp(exponent) * sqrt(2*M_PI*R*R) / sqrt(2*M_PI*sigma2);
		} else {
			overlap_integral = 0.0;
		}

		mean_variance(kpred, g.area*overlap_integral,
		              &partiality, &mean_kpred, &M2_kpred);

	}

	/* Revert the 'Lorentz' factor */
	partiality *= sqrt( ( R*R + M2_k/sumw_k) / ( R*R ) );
	if ( isnan(partiality) ) partiality = 0.0;

	if ( (updateme == NULL) && ( partiality < min_partiality ) ) return NULL;

	/* Calculate excitation error */
	knom = 1.0/image->lambda;
	dcs = distance3d(0.0, 0.0, -knom, xl, yl, zl);
	exerr = 1.0/image->lambda - dcs;  /* Loss of precision */

	/* Could also estimate excitation error like this, but it loses the sign
	 * (which is needed elsewhere):
	 * exerr = R * sqrt(-2*log(partiality)); */

	if ( updateme == NULL ) {
		refl = reflection_new(h, k, l);
	} else {
		refl = updateme;
	}

	/* If we are updating a previous reflection, assume it stays
	 * on the same panel and calculate the new position even if it's
	 * fallen off the edge of the panel. */
	if ( (image->detgeom != NULL) && (updateme != NULL) ) {

		double fs, ss;
		double det_shift_x, det_shift_y;

		assert(get_panel_number(updateme) <= image->detgeom->n_panels);
		crystal_get_det_shift(cryst, &det_shift_x, &det_shift_y);
		locate_peak_on_panel(xl, yl, zl, mean_kpred,
		                     &image->detgeom->panels[get_panel_number(updateme)],
		                     det_shift_x, det_shift_y,
		                     &fs, &ss);
		set_detector_pos(refl, fs, ss);

	}

	/* otherwise, calculate position if we have a detector structure, and
	 * if we don't then just make do with partiality calculation */
	if ( (image->detgeom != NULL) && (updateme == NULL) ) {

		double fs, ss;        /* position on detector */
		signed int p;         /* panel number */
		double det_shift_x, det_shift_y;

		crystal_get_det_shift(cryst, &det_shift_x, &det_shift_y);
		p = locate_peak(xl, yl, zl, mean_kpred,
		                image->detgeom,
		                det_shift_x, det_shift_y,
		                &fs, &ss);
		if ( p == -1 ) {
			reflection_free(refl);
			return NULL;
		}
		set_detector_pos(refl, fs, ss);
		set_panel_number(refl, p);

	}

	khalf = (- xl*xl - yl*yl - zl*zl) / (2.0*zl);
	set_kpred(refl, mean_kpred);
	set_khalf(refl, khalf);
	set_exerr(refl, exerr);
	set_lorentz(refl, 1.0);
	set_symmetric_indices(refl, h, k, l);
	set_redundancy(refl, 1);
	set_partiality(refl, partiality);

	return refl;
}


/**
 * \param cryst: A \ref Crystal
 * \param image: An image structure
 * \param max_res: Maximum resolution to predict to (m^-1)
 *
 * Calculates reflection positions for \p crys, as seen in \p image,
 * up to maximum 1/d value \p max_res
 *
 * \returns A list of predicted reflections
 */
RefList *predict_to_res(Crystal *cryst, struct image *image, double max_res)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	RefList *reflections;
	int hmax, kmax, lmax;
	double mres;
	signed int h, k, l;
	UnitCell *cell;

	cell = crystal_get_cell(cryst);
	if ( cell == NULL ) return NULL;

	reflections = reflist_new();

	/* Cell angle check from Foadi and Evans (2011) */
	if ( !cell_is_sensible(cell) ) {
		ERROR("Invalid unit cell parameters given to"
		      " predict_to_res()\n");
		cell_print(cell);
		return NULL;
	}

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	mres = detgeom_max_resolution(image->detgeom, image->lambda);
	if ( mres > max_res ) mres = max_res;

	hmax = mres * modulus(ax, ay, az);
	kmax = mres * modulus(bx, by, bz);
	lmax = mres * modulus(cx, cy, cz);

	if ( (hmax >= 512) || (kmax >= 512) || (lmax >= 512) ) {
		ERROR("Unit cell is too large - will only integrate reflections"
		      " up to 511th order.\n");
		cell_print(cell);
		if ( hmax >= 512 ) hmax = 511;
		if ( kmax >= 512 ) kmax = 511;
		if ( lmax >= 512 ) lmax = 511;
	}

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	for ( h=-hmax; h<=hmax; h++ ) {
	for ( k=-kmax; k<=kmax; k++ ) {
	for ( l=-lmax; l<=lmax; l++ ) {

		Reflection *refl;
		double xl, yl, zl;

		if ( forbidden_reflection(cell, h, k, l) ) continue;
		if ( 2.0*resolution(cell, h, k, l) > max_res ) continue;

		/* Get the coordinates of the reciprocal lattice point */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		refl = check_reflection(image, cryst,
		                        h, k, l, xl, yl, zl, NULL);

		if ( refl != NULL ) {
			add_refl_to_list(refl, reflections);
		}

	}
	}
	}

	return reflections;
}


static void set_unity_partialities(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;

	if ( list == NULL ) {
		ERROR("No reflections for partiality calculation!\n");
		return;
	}
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		set_partiality(refl, 1.0);
		set_lorentz(refl, 1.0);
	}
}


static void set_random_partialities(RefList *list, int image_serial)
{
	Reflection *refl;
	RefListIterator *iter;

	if ( list == NULL ) {
		ERROR("No reflections for partiality calculation!\n");
		return;
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		get_symmetric_indices(refl, &h, &k, &l);
		set_partiality(refl, random_partiality(h, k, l, image_serial));
		set_lorentz(refl, 1.0);
	}
}


static double do_integral(double q2, double zl, double R,
                          Spectrum *spectrum, char *verbose)
{
	int i;
	double kmin, kmax, kstart, kfinis;
	double inc;
	double total = 0.0;
	double k0, k1;
	const int SAMPLES = 50;  /* Number of samples for integration */
	FILE *fh = NULL;

	assert(R*R < q2);
	assert(R > 0.0);

	/* Range over which P is different from zero */
	k0 = (R*R - q2)/(2.0*(zl+R));
	k1 = (R*R - q2)/(2.0*(zl-R));

	/* Check for reflections partially "behind the beam" */
	if ( k0 < 0.0 ) k0 = +INFINITY;
	if ( k1 < 0.0 ) k1 = +INFINITY;

	/* Range over which E is significantly different from zero */
	spectrum_get_range(spectrum, &kmin, &kmax);

	/* Calculate range over which E*P is different from zero */
	if ( k0 < k1 ) {
		STATUS("%e %e\n", k0, k1);
		STATUS("%e %e %e\n", q2, zl, R);
		STATUS("%e %e\n", kmin, kmax);
	}
	assert((isinf(k0) && isinf(k1)) || (k0 > k1));
	assert(kmax > kmin);
	if ( kmin < k1 ) {
		if ( kmax < k1 ) {
			/* Case 1 */
			return 0.0;
		} else if ( kmax < k0 ) {
			/* Case 2 (kmax > k1)*/
			kstart = k1;   kfinis = kmax;
		} else {
			/* Case 3 (kmax > k1 and kmax > k0) */
			kstart = k1;   kfinis = k0;
		}
	} else if ( kmin < k0 ) { /* kmin > k1 */
		if ( kmax < k0 ) {
			/* Case 4 */
			kstart = kmin;   kfinis = kmax;
		} else {
			/* Case 5  (kmax > k0) */
			kstart = kmin;   kfinis = k0;
		}
	} else {
		/* Case 6 (kmin > k1 and (kmax>)kmin > k0) */
		return 0.0;
	}

	if ( kstart < 0.0 ) kstart = kmin;
	if ( kfinis < 0.0 ) kfinis = kmax;

	inc = (kfinis - kstart) / SAMPLES;

	if ( verbose ) {
		char fn[64];
		snprintf(fn, 63, "partial%s.graph", verbose);
		fh = fopen(fn, "w");
		fprintf(fh, "  n    p      wavelength   E           P\n");
		STATUS("k1/2 = %e m^-1\n", -q2/(2.0*zl));
		STATUS(" (wavelength %e m)\n", 1.0/(-q2/(2.0*zl)));
		STATUS("Reflection k goes from %e to %e m^-1\n", k1, k0);
		STATUS(" (wavelengths from %e to %e m\n", 1.0/k1, 1.0/k0);
		STATUS("Beam goes from %e to %e m^-1\n", kmin, kmax);
		STATUS(" (wavelengths from %e to %e m\n", 1.0/kmin, 1.0/kmax);
		STATUS("Integration goes from %e to %e m^-1\n", kstart, kfinis);
		STATUS(" (wavelengths from %e to %e m\n", 1.0/kstart, 1.0/kfinis);
	}

	for ( i=0; i<SAMPLES; i++ ) {

		double p, kp;
		double E, P;

		kp = kstart + i*inc;
		double pref = sqrt(q2 + kp*kp + 2.0*zl*kp)/(2.0*R);
		p = pref + 0.5 - kp/(2.0*R);

		/* Spectral energy term */
		E = spectrum_get_density_at_k(spectrum, kp);

		/* RLP profile term */
		P = 4.0*p * (1.0 - p);

		total += E*P*inc;

		if ( fh != NULL ) {
			fprintf(fh, "%3i %f %e %e %e\n", i, p, 1.0/kp, E, P);
		}
	}

	if ( fh != NULL ) fclose(fh);

	if ( isnan(total) ) {
		ERROR("NaN total!\n");
		STATUS("k1/2 = %e m^-1\n", -q2/(2.0*zl));
		STATUS(" (wavelength %e m)\n", 1.0/(-q2/(2.0*zl)));
		STATUS("Reflection k goes from %e to %e m^-1\n", k1, k0);
		STATUS(" (wavelengths from %e to %e m\n", 1.0/k1, 1.0/k0);
		STATUS("Beam goes from %e to %e m^-1\n", kmin, kmax);
		STATUS(" (wavelengths from %e to %e m\n", 1.0/kmin, 1.0/kmax);
		STATUS("Integration goes from %e to %e m^-1\n", kstart, kfinis);
		STATUS(" (wavelengths from %e to %e m\n", 1.0/kstart, 1.0/kfinis);
	}

	return total;
}


static void ginn_spectrum_partialities(RefList *list, Crystal *cryst, struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;
	double r0, m;
	UnitCell *cell;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	if ( list == NULL ) {
		ERROR("No reflections for partiality calculation!\n");
		return;
	}

	cell = crystal_get_cell(cryst);
	if ( cell == NULL ) {
		ERROR("No unit cell for partiality calculation!\n");
		return;
	}
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	r0 = fabs(crystal_get_profile_radius(cryst));
	m = crystal_get_mosaicity(cryst);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double R;
		signed int h, k, l;
		double xl, yl, zl;
		double q2;
		double total, norm;

		get_symmetric_indices(refl, &h, &k, &l);
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		/* Radius of rlp profile */
		q2 = xl*xl + yl*yl + zl*zl;

		R = r0 + m * sqrt(q2);

		//char tmp[256];
		//snprintf(tmp, 255, "-%i,%i,%i", h, k, l);
		char *tmp = NULL;

		total = do_integral(q2, zl, R, image->spectrum, tmp);
		norm = do_integral(q2, -0.5*q2*image->lambda, R,
		                   image->spectrum, NULL);

	        set_partiality(refl, total/norm);
		set_lorentz(refl, 1.0);

		assert(total <= 2.0*norm);

	}
}


static void ewald_offset_partialities(RefList *list, Crystal *cryst, struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;
	double r0, m;
	UnitCell *cell;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	if ( list == NULL ) {
		ERROR("No reflections for partiality calculation!\n");
		return;
	}

	cell = crystal_get_cell(cryst);
	if ( cell == NULL ) {
		ERROR("No unit cell for partiality calculation!\n");
		return;
	}
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	r0 = fabs(crystal_get_profile_radius(cryst));
	m = crystal_get_mosaicity(cryst);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double xl, yl, zl;
		double q2, R, t;

		get_symmetric_indices(refl, &h, &k, &l);
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		/* Radius of rlp profile */
		q2 = xl*xl + yl*yl + zl*zl;
		R = r0 + m * sqrt(q2);

		/* Excitation error */
		t = get_exerr(refl);

	        set_partiality(refl, exp(-(t*t)/(R*R)));
		set_lorentz(refl, 1.0);

	}
}


/**
 * \param list A \ref RefList
 * \param cryst A \ref Crystal
 * \param image An image structure
 * \param pmodel A \ref PartialityModel
 *
 * Calculates the partialities for the reflections in \p list, given the current
 * state of \p cryst and \p image.
 *
 * If \p pmodel is PMODEL_RANDOM or PMODEL_UNITY, then \p cryst can be NULL.
 * If \p pmodel is PMODEL_UNITY, then \p image can also be NULL.
 *
 * You must not have changed the crystal or image parameters since you last
 * called \ref predict_to_res or \ref update_predictions, because this function
 * relies on the limiting wavelength values calculated by those functions.
 */
void calculate_partialities(RefList *list, Crystal *cryst, struct image *image,
                            PartialityModel pmodel)
{
	switch ( pmodel ) {

		case PMODEL_UNITY :
		set_unity_partialities(list);
		break;

		case PMODEL_XSPHERE :
		ginn_spectrum_partialities(list, cryst, image);
		break;

		case PMODEL_OFFSET :
		ewald_offset_partialities(list, cryst, image);
		break;

		case PMODEL_RANDOM :
		set_random_partialities(list, image->serial);
		break;

		case PMODEL_GGPM :
		/* Do nothing, because this is the model used for prediction,
		 * so the partialities have been set by check_reflection */
		break;

		default :
		ERROR("Unknown partiality model %i\n", pmodel);
		break;

	}
}


/**
 * \param list A \ref RefList
 * \param cryst A \ref Crystal
 * \param image An image structure
 *
 * Updates the predicted reflections (positions and excitation errors, but not
 * the actual partialities) in \p list, to match the current statea of
 * \p crys as seen in \p image.
 *
 * If you need to update the partialities as well, call
 * \ref calculate_partialities afterwards.
 */
void update_predictions(RefList *list, Crystal *cryst, struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(crystal_get_cell(cryst), &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz, &csx, &csy, &csz);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double xl, yl, zl;
		signed int h, k, l;

		get_symmetric_indices(refl, &h, &k, &l);

		/* Get the coordinates of the reciprocal lattice point */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		check_reflection(image, cryst, h, k, l,
		                 xl, yl, zl, refl);

	}
}


struct polarisation parse_polarisation(const char *text)
{
	struct polarisation p;
	char angle[14];
	char deg[14];
	char frac[14];
	int i, n;

	if ( strlen(text) > 13 ) goto err;

	if ( strcmp(text, "none") == 0 ) {
		p.fraction = 0.5;
		p.angle = 0.0;
		p.disable = 1;
		return p;
	}

	p.fraction = 1.0;
	p.disable = 0;

	i = 0;
	n = 0;
	while ( isdigit(text[i]) ) {
		angle[n++] = text[i++];
	}
	angle[n] = '\0';

	n = 0;
	while ( isalpha(text[i]) ) {
		deg[n++] = text[i++];
	}
	deg[n] = '\0';

	n = 0;
	while ( isdigit(text[i]) ) {
		frac[n++] = text[i++];
	}
	frac[n] = '\0';

	if ( text[i] != '\0' ) goto err;

	if ( frac[0] != '\0' ) {
		int nfrac = atoi(frac);
		if ( nfrac > 100 ) goto err;
		if ( nfrac < 0 ) goto err;
		p.fraction = nfrac / 100.0;
	}

	if ( strcmp(deg, "deg") == 0 ) {
		p.angle = deg2rad(atoi(angle));
	} else {
		if ( angle[0] != '\0' ) goto err;
		if ( strncmp(deg, "horiz", 5) == 0 ) p.angle = 0.0;
		if ( strncmp(deg, "vert", 4) == 0 ) p.angle = M_PI_2;
	}

	return p;

err:
	ERROR("Invalid polarisation '%s'\n", text);
	p.fraction = 0.5;
	p.angle = 0.0;
	return p;
}


void polarisation_correction(RefList *list, UnitCell *cell,
                             struct polarisation p)
{
	Reflection *refl;
	RefListIterator *iter;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	if ( p.disable ) {
		return;
	}

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pol;
		double intensity, sigma;
		double xl, yl, zl;
		double phi, tt, kpred;
		signed int h, k, l;
		const double P = p.fraction;  /* degree of polarisation */

		get_symmetric_indices(refl, &h, &k, &l);

		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		kpred = get_kpred(refl);
		tt = angle_between(0.0, 0.0, 1.0, xl, yl, zl+kpred);
		phi = atan2(yl, xl) - p.angle;

		pol =         P*(1.0 - cos(phi)*cos(phi)*sin(tt)*sin(tt))
		     +  (1.0-P)*(1.0 - sin(phi)*sin(phi)*sin(tt)*sin(tt));

		intensity = get_intensity(refl);
		sigma = get_esd_intensity(refl);
		set_intensity(refl, intensity / pol);
		set_esd_intensity(refl, sigma / pol);
	}
}
