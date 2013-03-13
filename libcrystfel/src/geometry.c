/*
 * geometry.c
 *
 * Geometry of diffraction
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
#include <fenv.h>

#include "utils.h"
#include "cell.h"
#include "cell-utils.h"
#include "image.h"
#include "peaks.h"
#include "beam-parameters.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "geometry.h"


static signed int locate_peak(double x, double y, double z, double k,
                              struct detector *det, double *xdap, double *ydap)
{
	int i;
	signed int found = -1;
	const double den = k + z;

	*xdap = -1;  *ydap = -1;

	for ( i=0; i<det->n_panels; i++ ) {

		double xd, yd;
		double fs, ss, plx, ply;
		struct panel *p;

		p = &det->panels[i];

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

		/* Now, is this on this panel? */
		if ( fs < p->min_fs ) continue;
		if ( fs > p->max_fs ) continue;
		if ( ss < p->min_ss ) continue;
		if ( ss > p->max_ss ) continue;

		/* If peak appears on multiple panels, reject it */
		if ( found != -1 ) return -1;

		/* Woohoo! */
		found = i;
		*xdap = fs;
		*ydap = ss;

	}

	return found;
}


static double partiality(double rlow, double rhigh, double r)
{
	double qlow, qhigh;
	double plow, phigh;

	/* Calculate degrees of penetration */
	qlow  = (rlow + r)/(2.0*r);
	qhigh = (rhigh + r)/(2.0*r);

	/* Convert to partiality */
	plow  = 3.0*pow(qlow,2.0)  - 2.0*pow(qlow,3.0);
	phigh = 3.0*pow(qhigh,2.0) - 2.0*pow(qhigh,3.0);

	return plow - phigh;
}


static Reflection *check_reflection(struct image *image, Crystal *cryst,
                                    signed int h, signed int k, signed int l,
                                    double xl, double yl, double zl)
{
	const int output = 0;
	double tl;
	double rlow, rhigh;     /* "Excitation error" */
	signed int p;           /* Panel number */
	double xda, yda;        /* Position on detector */
	double part;            /* Partiality */
	int clamp_low, clamp_high;
	double klow, khigh;    /* Wavenumber */
	Reflection *refl;
	double cet, cez;
	double pr;
	double L;

	/* Don't predict 000 */
	if ( abs(h)+abs(k)+abs(l) == 0 ) return NULL;

	pr = crystal_get_profile_radius(cryst);

	/* "low" gives the largest Ewald sphere (wavelength short => k large)
	 * "high" gives the smallest Ewald sphere (wavelength long => k small)
	 */
	klow = 1.0/(image->lambda - image->lambda*image->bw/2.0);
	khigh = 1.0/(image->lambda + image->lambda*image->bw/2.0);

	/* If the point is looking "backscattery", reject it straight away */
	if ( zl < -khigh/2.0 ) return NULL;

	tl = sqrt(xl*xl + yl*yl);

	cet = -sin(image->div/2.0) * khigh;
	cez = -cos(image->div/2.0) * khigh;
	rhigh = khigh - distance(cet, cez, tl, zl);  /* Loss of precision */

	cet =  sin(image->div/2.0) * klow;
	cez = -cos(image->div/2.0) * klow;
	rlow = klow - distance(cet, cez, tl, zl);  /* Loss of precision */

	if ( (signbit(rlow) == signbit(rhigh))
	     && (fabs(rlow) > pr)
	     && (fabs(rhigh) > pr) ) return NULL;

	if ( rlow < rhigh ) {
		ERROR("Reflection with rlow < rhigh!\n");
		ERROR("%3i %3i %3i  rlow = %e, rhigh = %e\n",
		      h, k, l, rlow, rhigh);
		ERROR("div = %e\n", image->div);
		return NULL;
	}

	/* Lorentz factor is determined direction from the r values, before
	 * clamping.  The multiplication by the profile radius is to make the
	 * correction factor vaguely near 1. */
	L = pr / (rlow - rhigh);

	/* If the "lower" Ewald sphere is a long way away, use the
	 * position at which the Ewald sphere would just touch the
	 * reflection.
	 *
	 * The six possible combinations of clamp_{low,high} (including
	 * zero) correspond to the six situations in Table 3 of Rossmann
	 * et al. (1979).
	 */
	clamp_low = 0;  clamp_high = 0;
	if ( rlow < -pr ) {
		rlow = -pr;
		clamp_low = -1;
	}
	if ( rlow > +pr ) {
		rlow = +pr;
		clamp_low = +1;
	}
	if ( rhigh < -pr ) {
		rhigh = -pr;
		clamp_high = -1;
	}
	if ( rhigh > +pr ) {
		rhigh = +pr;
		clamp_high = +1;
	}

	/* Calculate partiality */
	part = partiality(rlow, rhigh, pr);

	/* Locate peak on detector. */
	p = locate_peak(xl, yl, zl, 1.0/image->lambda, image->det, &xda, &yda);
	if ( p == -1 ) return NULL;

	/* Add peak to list */
	refl = reflection_new(h, k, l);
	set_detector_pos(refl, 0.0, xda, yda);
	set_partial(refl, rlow, rhigh, part, clamp_low, clamp_high);
	set_lorentz(refl, L);
	set_symmetric_indices(refl, h, k, l);
	set_redundancy(refl, 1);

	if ( output ) {
		printf("%3i %3i %3i %6f (at %5.2f,%5.2f) %5.2f\n",
		       h, k, l, 0.0, xda, yda, part);
	}

	return refl;
}


RefList *find_intersections(struct image *image, Crystal *cryst)
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

	hmax = mres * modulus(ax, ay, az);
	kmax = mres * modulus(bx, by, bz);
	lmax = mres * modulus(cx, cy, cz);

	if ( (hmax >= 256) || (kmax >= 256) || (lmax >= 256) ) {
		ERROR("Unit cell is stupidly large.\n");
		cell_print(cell);
		if ( hmax >= 256 ) hmax = 255;
		if ( kmax >= 256 ) kmax = 255;
		if ( lmax >= 256 ) lmax = 255;
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

		/* Get the coordinates of the reciprocal lattice point */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		refl = check_reflection(image, cryst, h, k, l, xl, yl, zl);

		if ( refl != NULL ) {
			add_refl_to_list(refl, reflections);
		}

	}
	}
	}

	return reflections;
}


RefList *select_intersections(struct image *image, Crystal *cryst)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	const double min_dist = 0.25;
	RefList *list;
	int i;

	/* Round towards nearest */
	fesetround(1);

	/* Cell basis vectors for this image */
	cell_get_cartesian(crystal_get_cell(cryst), &ax, &ay, &az,
	                   &bx, &by, &bz, &cx, &cy, &cz);

	list = reflist_new();
	if ( list == NULL ) return NULL;

	/* Loop over peaks, checking proximity to nearest reflection */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec q;
		double h, k, l, hd, kd, ld;
		double dsq;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		/* Reciprocal space position of found peak */
		q = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);

		/* Decimal and fractional Miller indices of nearest
		 * reciprocal lattice point */
		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		/* Check distance */
		dsq = pow(h-hd, 2.0) + pow(k-kd, 2.0) + pow(l-ld, 2.0);

		if ( sqrt(dsq) < min_dist ) {

			Reflection *refl;

			refl = add_refl(list, h, k, l);
			set_detector_pos(refl, sqrt(dsq), f->fs, f->ss);

		}

	}

	return list;
}


/* Calculate partialities and apply them to the image's reflections */
void update_partialities(Crystal *cryst, PartialityModel pmodel)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *predicted;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct image *image = crystal_get_image(cryst);

	cell_get_reciprocal(crystal_get_cell(cryst), &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz, &csx, &csy, &csz);

	/* Scratch list to give check_reflection() something to add to */
	predicted = reflist_new();

	for ( refl = first_refl(crystal_get_reflections(cryst), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *vals;
		double r1, r2, p, x, y;
		double xl, yl, zl;
		signed int h, k, l;
		int clamp1, clamp2;

		get_symmetric_indices(refl, &h, &k, &l);

		/* Get the coordinates of the reciprocal lattice point */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		vals = check_reflection(image, cryst, h, k, l, xl, yl, zl);

		if ( vals == NULL ) {
			set_redundancy(refl, 0);
			continue;
		}
		set_redundancy(refl, 1);

		/* Transfer partiality stuff */
		get_partial(vals, &r1, &r2, &p, &clamp1, &clamp2);
		set_partial(refl, r1, r2, p, clamp1, clamp2);

		/* Transfer detector location */
		get_detector_pos(vals, &x, &y);
		set_detector_pos(refl, 0.0, x, y);

		reflection_free(vals);
	}

	reflist_free(predicted);
}
