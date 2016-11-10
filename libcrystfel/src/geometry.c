/*
 * geometry.c
 *
 * Geometry of diffraction
 *
 * Copyright © 2012-2016 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2016 Thomas White <taw@physics.org>
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
#include <fenv.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>

#include "utils.h"
#include "cell.h"
#include "cell-utils.h"
#include "image.h"
#include "peaks.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "geometry.h"


static int locate_peak_on_panel(double x, double y, double z, double k,
                                struct panel *p,
                                double *pfs, double *pss)
{
	double ctt, tta, phi;
	gsl_vector *v;
	gsl_vector *t;
	gsl_matrix *M;
	double fs, ss, one_over_mu;

	/* Calculate 2theta (scattering angle) and azimuth (phi) */
	tta = atan2(sqrt(x*x+y*y), k+z);
	ctt = cos(tta);
	phi = atan2(y, x);

	/* Set up matrix equation */
	M = gsl_matrix_alloc(3, 3);
	v = gsl_vector_alloc(3);
	t = gsl_vector_alloc(3);
	if ( (M==NULL) || (v==NULL) || (t==NULL) ) {
		ERROR("Failed to allocate vectors for prediction\n");
		return 0;
	}

	gsl_vector_set(t, 0, sin(tta)*cos(phi));
	gsl_vector_set(t, 1, sin(tta)*sin(phi));
	gsl_vector_set(t, 2, ctt);

	gsl_matrix_set(M, 0, 0, p->cnx);
	gsl_matrix_set(M, 0, 1, p->fsx);
	gsl_matrix_set(M, 0, 2, p->ssx);
	gsl_matrix_set(M, 1, 0, p->cny);
	gsl_matrix_set(M, 1, 1, p->fsy);
	gsl_matrix_set(M, 1, 2, p->ssy);
	gsl_matrix_set(M, 2, 0, p->clen*p->res);
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

	/* Now, is this on this panel? */
	if ( fs < 0.0 ) return 0;
	if ( fs >= p->w ) return 0;
	if ( ss < 0.0 ) return 0;
	if ( ss >= p->h ) return 0;

	return 1;
}

static signed int locate_peak(double x, double y, double z, double k,
                              struct detector *det, double *pfs, double *pss)
{
	int i;

	*pfs = -1;  *pss = -1;

	for ( i=0; i<det->n_panels; i++ ) {

		struct panel *p;

		p = &det->panels[i];

		if ( locate_peak_on_panel(x, y, z, k, p, pfs, pss) ) {

			/* Woohoo! */
			return i;

		}

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


static double random_partiality(signed int h, signed int k, signed int l,
                                int serial)
{
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long int seed;
	double p;
	int i;

	gsl_rng_set(rng, serial);
	seed = gsl_rng_get(rng);
	gsl_rng_set(rng, seed);

	for ( i=0; i<abs(h)+1; i++ ) {
		seed = gsl_rng_get(rng);
	}
	gsl_rng_set(rng, seed);
	if ( h >= 0 ) {
		seed = gsl_rng_get(rng);
	}
	seed = gsl_rng_get(rng);
	gsl_rng_set(rng, seed);

	for ( i=0; i<abs(k)+1; i++ ) {
		seed = gsl_rng_get(rng);
	}
	gsl_rng_set(rng, seed);
	if ( k >= 0 ) {
		seed = gsl_rng_get(rng);
	}
	seed = gsl_rng_get(rng);
	gsl_rng_set(rng, seed);

	for ( i=0; i<abs(l)+1; i++ ) {
		seed = gsl_rng_get(rng);
	}
	gsl_rng_set(rng, seed);
	if ( l >= 0 ) {
		seed = gsl_rng_get(rng);
	}
	seed = gsl_rng_get(rng);
	gsl_rng_set(rng, seed);

	p = gsl_rng_uniform(rng);
	gsl_rng_free(rng);
	return p;
}


static Reflection *check_reflection(struct image *image, Crystal *cryst,
                                    signed int h, signed int k, signed int l,
                                    double xl, double yl, double zl,
                                    Reflection *updateme)
{
	Reflection *refl;
	double R, top;
	double kmin, kmax, k0, knom, k1;
	double dcs, exerr;

	/* Don't predict 000 */
	if ( (updateme == NULL) && (abs(h)+abs(k)+abs(l) == 0) ) return NULL;

	/* Calculate the limiting wavelengths, lambda0 and lambda1
	 * = 1/k0 and 1/k1 respectively */
	R = crystal_get_profile_radius(cryst);
	top = R*R - xl*xl - yl*yl - zl*zl;
	k0 = top/(2.0*(zl+R));
	k1 = top/(2.0*(zl-R));

	/* The reflection is excited if any of the reflection is within 2sigma
	 * of the nominal * wavelength of the X-ray beam
	 * (NB image->bw is full width) */
	kmin = 1.0/(image->lambda + image->lambda*image->bw);
	knom = 1.0/image->lambda;
	kmax = 1.0/(image->lambda - image->lambda*image->bw);
	if ( (k1>kmax) || (k0<kmin) ) return NULL;

	/* Calculate excitation error */
	dcs = distance3d(0.0, 0.0, -knom, xl, yl, zl);
	exerr = 1.0/image->lambda - dcs;  /* Loss of precision */

	if ( updateme == NULL ) {
		refl = reflection_new(h, k, l);
	} else {
		refl = updateme;
	}

	/* If we are updating a previous reflection, assume it stays
	 * on the same panel and calculate the new position even if it's
	 * fallen off the edge of the panel. */
	if ( (image->det != NULL) && (updateme != NULL) ) {

		double fs, ss;
		locate_peak_on_panel(xl, yl, zl, knom,
		                     get_panel(updateme), &fs, &ss);
		set_detector_pos(refl, fs, ss);

	}

	/* Otherwise, calculate position if we have a detector structure, and
	 * if we don't then just make do with partiality calculation */
	if ( (image->det != NULL) && (updateme == NULL) ) {

		double fs, ss;        /* Position on detector */
		signed int p;         /* Panel number */
		p = locate_peak(xl, yl, zl, knom,
		                image->det, &fs, &ss);
		if ( p == -1 ) {
			reflection_free(refl);
			return NULL;
		}
		set_detector_pos(refl, fs, ss);
		set_panel(refl, &image->det->panels[p]);

	}

	set_kpred(refl, knom);
	set_exerr(refl, exerr);
	set_lorentz(refl, 1.0);
	set_symmetric_indices(refl, h, k, l);
	set_redundancy(refl, 1);

	return refl;
}


double r_gradient(UnitCell *cell, int k, Reflection *refl, struct image *image)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xl, yl, zl;
	signed int hs, ks, ls;
	double tl, phi, azi;

	get_symmetric_indices(refl, &hs, &ks, &ls);

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	xl = hs*asx + ks*bsx + ls*csx;
	yl = hs*asy + ks*bsy + ls*csy;
	zl = hs*asz + ks*bsz + ls*csz;

	tl = sqrt(xl*xl + yl*yl);
	phi = angle_between_2d(tl, zl+1.0/image->lambda, 0.0, 1.0); /* 2theta */
	azi = atan2(yl, xl); /* azimuth */

	switch ( k ) {

		case GPARAM_ASX :
		return - hs * sin(phi) * cos(azi);

		case GPARAM_BSX :
		return - ks * sin(phi) * cos(azi);

		case GPARAM_CSX :
		return - ls * sin(phi) * cos(azi);

		case GPARAM_ASY :
		return - hs * sin(phi) * sin(azi);

		case GPARAM_BSY :
		return - ks * sin(phi) * sin(azi);

		case GPARAM_CSY :
		return - ls * sin(phi) * sin(azi);

		case GPARAM_ASZ :
		return - hs * cos(phi);

		case GPARAM_BSZ :
		return - ks * cos(phi);

		case GPARAM_CSZ :
		return - ls * cos(phi);

		case GPARAM_DETX :
		case GPARAM_DETY :
		case GPARAM_CLEN :
		return 0.0;

	}

	ERROR("No r gradient defined for parameter %i\n", k);
	abort();
}


RefList *predict_to_res(Crystal *cryst, double max_res)
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
		      " find_intersections()\n");
		cell_print(cell);
		return NULL;
	}

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	mres = largest_q(crystal_get_image(cryst));
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

		refl = check_reflection(crystal_get_image(cryst), cryst,
		                        h, k, l, xl, yl, zl, NULL);

		if ( refl != NULL ) {
			add_refl_to_list(refl, reflections);
		}

	}
	}
	}

	return reflections;
}


static void set_unity_partialities(Crystal *cryst)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;

	list = crystal_get_reflections(cryst);
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


static void set_random_partialities(Crystal *cryst)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	struct image *image;

	list = crystal_get_reflections(cryst);
	if ( list == NULL ) {
		ERROR("No reflections for partiality calculation!\n");
		return;
	}

	image = crystal_get_image(cryst);
	if ( image == NULL ) {
		ERROR("No image structure for partiality calculation!\n");
		return;
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		get_symmetric_indices(refl, &h, &k, &l);
		set_partiality(refl, random_partiality(h, k, l, image->serial));
		set_lorentz(refl, 1.0);
	}
}


static double do_integral(double q2, double zl, double R,
                          double lambda, double sig, int verbose)
{
	int i;
	double kmin, kmax, kstart, kfinis;
	double inc;
	double total = 0.0;
	double k0, k1;
	const int SAMPLES = 50;  /* Number of samples for integration */
	const double N = 1.5;  /* Pointiness of spectrum */
	FILE *fh = NULL;

	k0 = (R*R - q2)/(2.0*(zl+R));
	k1 = (R*R - q2)/(2.0*(zl-R));

	/* Range over which E is significantly different from zero */
	kmin = 1.0 / (lambda + 5.0*sig);
	kmax = 1.0 / (lambda - 5.0*sig);

	kstart = kmin > k1 ? kmin : k1;
	kfinis = (k0 < 0.0) || (kmax < k0) ? kmax : k0;
	inc = (kfinis - kstart) / SAMPLES;

	if ( verbose ) {
		char fn[64];
		snprintf(fn, 63, "partial%i.graph", verbose);
		fh = fopen(fn, "w");
		fprintf(fh, "  n    p      wavelength   E           P\n");
		STATUS("Nominal k = %e m^-1\n", 1.0/lambda);
		STATUS(" (wavelength %e m)\n", lambda);
		STATUS("Bandwidth %e m\n", sig);
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

		double p, kp, lrel;
		double E, P;

		kp = kstart + i*inc;
		double pref = sqrt(q2 + kp*kp + 2.0*zl*kp)/(2.0*R);
		p = pref + 0.5 - kp/(2.0*R);

		/* Spectral energy term */
		lrel = fabs(1.0/kp - lambda);
		E = exp(-0.5 * pow(lrel / sig, N));
		E /= sqrt(2.0 * M_PI * sig);

		/* RLP profile term */
		P = 4.0*p * (1.0 - p);

		total += E*P*inc;

		if ( fh != NULL ) {
			fprintf(fh, "%3i %f %e %e %e\n", i, p, 1.0/kp, E, P);
		}
	}

	if ( fh != NULL ) fclose(fh);

	return total;
}



static void ginn_spectrum_partialities(Crystal *cryst)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	double r0, m, lambda, sig;
	struct image *image;
	UnitCell *cell;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	list = crystal_get_reflections(cryst);
	if ( list == NULL ) {
		ERROR("No reflections for partiality calculation!\n");
		return;
	}

	image = crystal_get_image(cryst);
	if ( image == NULL ) {
		ERROR("No image for partiality calculation!\n");
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

	r0 = crystal_get_profile_radius(cryst);
	m = crystal_get_mosaicity(cryst);
	lambda = image->lambda;
	sig = image->bw * lambda;

	for ( refl = first_refl(crystal_get_reflections(cryst), &iter);
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

		total = do_integral(q2, zl, R, lambda, sig, 0);
		norm = do_integral(q2, -0.5*q2*lambda, R, lambda, sig, 0);

	        set_partiality(refl, total/norm);
		set_lorentz(refl, 1.0);

		if ( total > 2.0*norm ) {
			/* Error! */
			do_integral(q2, zl, R, lambda, sig, 1);
			do_integral(q2, -0.5*q2*lambda, R, lambda, sig, 2);
			abort();
		}

	}
}


/**
 * calculate_partialities:
 * @cryst: A %Crystal
 * @pmodel: A %PartialityModel
 *
 * Calculates the partialities for the reflections in @cryst, given the current
 * crystal and image parameters.  The crystal's image and reflection lists
 * must be set.  The specified %PartialityModel will be used.
 *
 * You must not have changed the crystal or image parameters since you last
 * called predict_to_res() or update_predictions(), because this function
 * relies on the limiting wavelength values calculated by those functions.
 */
void calculate_partialities(Crystal *cryst, PartialityModel pmodel)
{
	switch ( pmodel ) {

		case PMODEL_UNITY :
		set_unity_partialities(cryst);
		break;

		case PMODEL_XSPHERE :
		ginn_spectrum_partialities(cryst);
		break;

		case PMODEL_RANDOM :
		set_random_partialities(cryst);
		break;

		default :
		ERROR("Unknown partiality model %i\n", pmodel);
		break;

	}
}


/**
 * update_predictions:
 * @cryst: A %Crystal
 *
 * Updates the predicted reflections (positions and excitation errors, but not
 * the actual partialities) of @cryst's reflections according to
 * the current state of the crystal (e.g. its unit cell parameters).
 *
 * If you need to update the partialities as well, call calculate_partialities()
 * afterwards.
 */
void update_predictions(Crystal *cryst)
{
	Reflection *refl;
	RefListIterator *iter;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct image *image = crystal_get_image(cryst);

	cell_get_reciprocal(crystal_get_cell(cryst), &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz, &csx, &csy, &csz);

	for ( refl = first_refl(crystal_get_reflections(cryst), &iter);
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


void polarisation_correction(RefList *list, UnitCell *cell, struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pol;
		double intensity;
		double xl, yl;
		signed int h, k, l;
		const double P = 1.0;  /* degree of polarisation */

		get_indices(refl, &h, &k, &l);

		xl = (h*asx + k*bsx + l*csx)*image->lambda;
		yl = (h*asy + k*bsy + l*csy)*image->lambda;

		pol = P*(1.0 - xl*xl) + (1.0-P)*(1.0 - yl*yl);

		intensity = get_intensity(refl);
		set_intensity(refl, intensity / pol);
	}
}


/* Returns dx_h/dP, where P = any parameter */
double x_gradient(int param, Reflection *refl, UnitCell *cell, struct panel *p)
{
	signed int h, k, l;
	double xl, zl, kpred;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	get_indices(refl, &h, &k, &l);
	kpred = get_kpred(refl);
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	xl = h*asx + k*bsx + l*csx;
	zl = h*asz + k*bsz + l*csz;

	switch ( param ) {

		case GPARAM_ASX :
		return h * p->clen / (kpred + zl);

		case GPARAM_BSX :
		return k * p->clen / (kpred + zl);

		case GPARAM_CSX :
		return l * p->clen / (kpred + zl);

		case GPARAM_ASY :
		return 0.0;

		case GPARAM_BSY :
		return 0.0;

		case GPARAM_CSY :
		return 0.0;

		case GPARAM_ASZ :
		return -h * xl * p->clen / (kpred*kpred + 2.0*kpred*zl + zl*zl);

		case GPARAM_BSZ :
		return -k * xl * p->clen / (kpred*kpred + 2.0*kpred*zl + zl*zl);

		case GPARAM_CSZ :
		return -l * xl * p->clen / (kpred*kpred + 2.0*kpred*zl + zl*zl);

		case GPARAM_DETX :
		return -1;

		case GPARAM_DETY :
		return 0;

		case GPARAM_CLEN :
		return xl / (kpred+zl);

	}

	ERROR("Positional gradient requested for parameter %i?\n", param);
	abort();
}


/* Returns dy_h/dP, where P = any parameter */
double y_gradient(int param, Reflection *refl, UnitCell *cell, struct panel *p)
{
	signed int h, k, l;
	double yl, zl, kpred;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	get_indices(refl, &h, &k, &l);
	kpred = get_kpred(refl);
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	yl = h*asy + k*bsy + l*csy;
	zl = h*asz + k*bsz + l*csz;

	switch ( param ) {

		case GPARAM_ASX :
		return 0.0;

		case GPARAM_BSX :
		return 0.0;

		case GPARAM_CSX :
		return 0.0;

		case GPARAM_ASY :
		return h * p->clen / (kpred + zl);

		case GPARAM_BSY :
		return k * p->clen / (kpred + zl);

		case GPARAM_CSY :
		return l * p->clen / (kpred + zl);

		case GPARAM_ASZ :
		return -h * yl * p->clen / (kpred*kpred + 2.0*kpred*zl + zl*zl);

		case GPARAM_BSZ :
		return -k * yl * p->clen / (kpred*kpred + 2.0*kpred*zl + zl*zl);

		case GPARAM_CSZ :
		return -l * yl * p->clen / (kpred*kpred + 2.0*kpred*zl + zl*zl);

		case GPARAM_DETX :
		return 0;

		case GPARAM_DETY :
		return -1;

		case GPARAM_CLEN :
		return yl / (kpred+zl);

	}

	ERROR("Positional gradient requested for parameter %i?\n", param);
	abort();
}
