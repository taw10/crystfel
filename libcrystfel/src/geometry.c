/*
 * geometry.c
 *
 * Geometry of diffraction
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
#include <fenv.h>
#include <gsl/gsl_sf_erf.h>

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
	const double den = k + z;
	double fs, ss, plx, ply, xd, yd;

	/* Coordinates of peak relative to central beam, in m */
	xd = p->clen * x / den;
	yd = p->clen * y / den;

	/* Convert to pixels */
	xd *= p->res;
	yd *= p->res;

	/* Convert to relative to the panel corner */
	plx = xd - p->cnx;
	ply = yd - p->cny;

	fs = p->xfs*plx + p->yfs*ply;
	ss = p->xss*plx + p->yss*ply;

	fs += p->min_fs;
	ss += p->min_ss;

	*pfs = fs;  *pss = ss;

	/* Now, is this on this panel? */
	if ( fs < p->min_fs ) return 0;
	if ( fs > p->max_fs ) return 0;
	if ( ss < p->min_ss ) return 0;
	if ( ss > p->max_ss ) return 0;

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


static double partiality(PartialityModel pmodel,
                         double rlow, double rhigh, double pr)
{
	double D = rlow - rhigh;

	/* Convert to partiality */
	switch ( pmodel ) {

		default:
		case PMODEL_UNITY:
		return 1.0;

		case PMODEL_SCSPHERE:
		return 4.0*sphere_fraction(rlow, rhigh, pr)*pr / (3.0*D);

		case PMODEL_SCGAUSSIAN:
		return 4.0*gaussian_fraction(rlow, rhigh, pr)*pr / (3.0*D);

	}
}


static Reflection *check_reflection(struct image *image, Crystal *cryst,
                                    PartialityModel pmodel,
                                    signed int h, signed int k, signed int l,
                                    double xl, double yl, double zl,
                                    Reflection *updateme)
{
	const int output = 0;
	double tl;
	double rlow, rhigh;     /* "Excitation error" */
	double part;            /* Partiality */
	double klow, khigh;    /* Wavenumber */
	Reflection *refl;
	double cet, cez;  /* Centre of Ewald sphere */
	double pr;
	double del;

	/* Don't predict 000 */
	if ( (updateme == NULL) && (abs(h)+abs(k)+abs(l) == 0) ) return NULL;

	pr = crystal_get_profile_radius(cryst);
	del = image->div + crystal_get_mosaicity(cryst);

	/* "low" gives the largest Ewald sphere (wavelength short => k large)
	 * "high" gives the smallest Ewald sphere (wavelength long => k small)
	 */
	klow = 1.0/(image->lambda - image->lambda*image->bw/2.0);
	khigh = 1.0/(image->lambda + image->lambda*image->bw/2.0);

	/* If the point is looking "backscattery", reject it straight away */
	if ( (updateme == NULL) && (zl < -khigh/2.0) ) return NULL;

	tl = sqrt(xl*xl + yl*yl);

	cet = -sin(del/2.0) * khigh;
	cez = -cos(del/2.0) * khigh;
	rhigh = khigh - distance(cet, cez, tl, zl);  /* Loss of precision */

	cet =  sin(del/2.0) * klow;
	cez = -cos(del/2.0) * klow;
	rlow = klow - distance(cet, cez, tl, zl);  /* Loss of precision */

	/* Condition for reflection to be excited at all */
	if ( (updateme == NULL)
	     && (signbit(rlow) == signbit(rhigh))
	     && (fabs(rlow) > pr)
	     && (fabs(rhigh) > pr) ) return NULL;

	/* Calculate partiality */
	part = partiality(pmodel, rlow, rhigh, pr);

	if ( updateme == NULL ) {
		refl = reflection_new(h, k, l);
	} else {
		refl = updateme;
	}

	/* If we are updating a previous reflection, assume it stays
	 * on the same panel and calculate the new position even if it's
	 * fallen off the edge of the panel. */
	if ( updateme != NULL ) {

		double fs, ss;
		locate_peak_on_panel(xl, yl, zl, 1.0/image->lambda,
		                     get_panel(updateme), &fs, &ss);
		set_detector_pos(refl, fs, ss);

	}

	/* Otherwise, calculate position if we have a detector structure, and
	 * if we don't then just make do with partiality calculation */
	if ( (image->det != NULL) && (updateme == NULL) ) {

		double fs, ss;        /* Position on detector */
		signed int p;         /* Panel number */
		p = locate_peak(xl, yl, zl, 1.0/image->lambda,
		                image->det, &fs, &ss);
		if ( p == -1 ) {
			reflection_free(refl);
			return NULL;
		}
		set_detector_pos(refl, fs, ss);
		set_panel(refl, &image->det->panels[p]);

	}

	if ( unlikely(rlow < rhigh) ) {
		ERROR("Reflection with rlow < rhigh!\n");
		ERROR("%3i %3i %3i  rlow = %e, rhigh = %e\n",
		      h, k, l, rlow, rhigh);
		ERROR("div + m = %e, R = %e, bw = %e\n", del, pr, image->bw);
		/* If we are updating, this is (kind of) OK */
		if ( updateme == NULL ) {
			reflection_free(refl);
			return NULL;
		}
	}

	set_partial(refl, rlow, rhigh, part);
	set_lorentz(refl, 1.0);
	set_symmetric_indices(refl, h, k, l);
	set_redundancy(refl, 1);

	if ( output ) {
		printf("%3i %3i %3i %6f %5.2f\n", h, k, l, 0.0, part);
	}

	return refl;
}


double r_gradient(UnitCell *cell, int k, Reflection *refl, struct image *image)
{
	double azi;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xl, yl, zl;
	signed int hs, ks, ls;
	double rlow, rhigh, p;
	double philow, phihigh, phi;
	double khigh, klow;
	double tl, cet, cez;

	get_partial(refl, &rlow, &rhigh, &p);

	get_symmetric_indices(refl, &hs, &ks, &ls);

	cell_get_reciprocal(cell, &asx, &asy, &asz,
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


RefList *find_intersections(struct image *image, Crystal *cryst,
                            PartialityModel pmodel)
{
	return find_intersections_to_res(image, cryst, pmodel, INFINITY);
}


RefList *find_intersections_to_res(struct image *image, Crystal *cryst,
                                   PartialityModel pmodel, double max_res)
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

	mres = largest_q(image);
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

		refl = check_reflection(image, cryst, pmodel,
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
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(crystal_get_reflections(cryst), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		set_partiality(refl, 1.0);
		set_lorentz(refl, 1.0);
	}
}


/* Calculate partialities and apply them to the image's reflections */
void update_partialities(Crystal *cryst, PartialityModel pmodel)
{
	Reflection *refl;
	RefListIterator *iter;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct image *image = crystal_get_image(cryst);

	if ( pmodel == PMODEL_UNITY ) {
		set_unity_partialities(cryst);
		return;
	}

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

		check_reflection(image, cryst, pmodel,
		                 h, k, l, xl, yl, zl, refl);

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
		double pol, pa, pb, phi, tt, ool;
		double intensity;
		double xl, yl, zl;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		/* Polarisation correction assuming 100% polarisation
		 * along the x direction */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		ool = 1.0 / image->lambda;
		tt = angle_between(0.0, 0.0, 1.0,  xl, yl, zl+ool);
		phi = atan2(yl, xl);
		pa = pow(sin(phi)*sin(tt), 2.0);
		pb = pow(cos(tt), 2.0);
		pol = 1.0 - 2.0*(1.0-pa) + (1.0+pb);

		intensity = get_intensity(refl);
		set_intensity(refl, intensity / pol);
	}
}
