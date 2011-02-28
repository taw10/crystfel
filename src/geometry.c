/*
 * geometry.c
 *
 * Geometry of diffraction
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <assert.h>

#include "utils.h"
#include "cell.h"
#include "image.h"
#include "peaks.h"
#include "beam-parameters.h"
#include "reflist.h"


#define MAX_CPEAKS (256 * 256)


static signed int locate_peak(double x, double y, double z, double k,
                              struct detector *det, double *xdap, double *ydap)
{
	int i;
	signed int found = -1;
	const double den = k + z;

	*xdap = -1;  *ydap = -1;

	for ( i=0; i<det->n_panels; i++ ) {

		double xd, yd, cl;
		double fs, ss, plx, ply;
		struct panel *p;

		p = &det->panels[i];

		/* Camera length for this panel */
		cl = p->clen;

		/* Coordinates of peak relative to central beam, in m */
		xd = cl * x / den;
		yd = cl * y / den;

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


static double excitation_error(double xl, double yl, double zl,
                               double ds, double k, double divergence)
{
	double tt, al;
	double r;
	double delta;

	tt = angle_between(0.0, 0.0, 1.0,  xl, yl, zl+k);
	al = M_PI_2 - asin(-zl/ds);

	r = ( ds * sin(al) / sin(tt) ) - k;

	delta = sqrt(2.0 * pow(ds, 2.0) * (1-cos(divergence)));
	if ( divergence > 0.0 ) {
		r += delta;
	} else {
		r -= delta;
	}

	return r;
}


static double partiality(double r1, double r2, double r)
{
	double q1, q2;
	double p1, p2;

	/* Calculate degrees of penetration */
	q1 = (r1 + r)/(2.0*r);
	q2 = (r2 + r)/(2.0*r);

	/* Convert to partiality */
	p1 = 3.0*pow(q1,2.0) - 2.0*pow(q1,3.0);
	p2 = 3.0*pow(q2,2.0) - 2.0*pow(q2,3.0);

	return p2 - p1;
}


static int check_reflection(struct image *image, double mres, int output,
                            RefList *reflections,
                            signed int h, signed int k, signed int l,
                            double asx, double asy, double asz,
                            double bsx, double bsy, double bsz,
                            double csx, double csy, double csz)
{
	double xl, yl, zl;
	double ds, ds_sq;
	double rlow, rhigh;     /* "Excitation error" */
	signed int p;           /* Panel number */
	double xda, yda;        /* Position on detector */
	int close, inside;
	double part;            /* Partiality */
	int clamp_low = 0;
	int clamp_high = 0;
	double bandwidth = image->bw;
	double divergence = image->div;
	double lambda = image->lambda;
	double klow, kcen, khigh;    /* Wavenumber */
	Reflection *refl;

	/* "low" gives the largest Ewald sphere,
	 * "high" gives the smallest Ewald sphere. */
	klow = 1.0/(lambda - lambda*bandwidth/2.0);
	kcen = 1.0/lambda;
	khigh = 1.0/(lambda + lambda*bandwidth/2.0);

	/* Get the coordinates of the reciprocal lattice point */
	zl = h*asz + k*bsz + l*csz;
	/* Throw out if it's "in front".  A tiny bit "in front" is OK. */
	if ( zl > image->profile_radius ) return 0;
	xl = h*asx + k*bsx + l*csx;
	yl = h*asy + k*bsy + l*csy;

	/* Calculate reciprocal lattice point modulus (and square) */
	ds_sq = modulus_squared(xl, yl, zl);  /* d*^2 */
	ds = sqrt(ds_sq);
	if ( ds > mres ) return 0;  /* Outside resolution range */

	/* Calculate excitation errors */
	rlow = excitation_error(xl, yl, zl, ds, klow, -divergence);
	rhigh = excitation_error(xl, yl, zl, ds, khigh, +divergence);

	/* Is the reciprocal lattice point close to either extreme of
	 * the sphere, maybe just outside the "Ewald volume"? */
	close = (fabs(rlow) < image->profile_radius)
	     || (fabs(rhigh) < image->profile_radius);

	/* Is the reciprocal lattice point somewhere between the
	 * extremes of the sphere, i.e. inside the "Ewald volume"? */
	inside = signbit(rlow) ^ signbit(rhigh);

	/* Can't be both inside and close */
	if ( inside ) close = 0;

	/* Neither?  Skip it. */
	if ( !(close || inside) ) return 0;

	/* If the "lower" Ewald sphere is a long way away, use the
	 * position at which the Ewald sphere would just touch the
	 * reflection. */
	if ( rlow < -image->profile_radius ) {
		rlow = -image->profile_radius;
		clamp_low = -1;
	}
	if ( rlow > +image->profile_radius ) {
		rlow = +image->profile_radius;
		clamp_low = +1;
	}
	/* Likewise the "higher" Ewald sphere */
	if ( rhigh < -image->profile_radius ) {
		rhigh = -image->profile_radius;
		clamp_high = -1;
	}
	if ( rhigh > +image->profile_radius ) {
		rhigh = +image->profile_radius;
		clamp_high = +1;
	}
	/* The six possible combinations of clamp_{low,high} (including
	 * zero) correspond to the six situations in Table 3 of Rossmann
	 * et al. (1979). */

	/* Calculate partiality and reject if too small */
	part = partiality(rlow, rhigh, image->profile_radius);
	if ( part < 0.1 ) return 0;

	/* Locate peak on detector. */
	p = locate_peak(xl, yl, zl, kcen, image->det, &xda, &yda);
	if ( p == -1 ) return 0;

	/* Add peak to list */
	refl = add_refl(reflections, h, k, l);
	set_detector_pos(refl, 0.0, xda, yda);
	set_partial(refl, rlow, rhigh, part, clamp_low, clamp_high);

	if ( output ) {
		printf("%3i %3i %3i %6f (at %5.2f,%5.2f) %5.2f\n",
		       h, k, l, 0.0, xda, yda, part);
	}

	return 1;
}


RefList *find_intersections(struct image *image, UnitCell *cell,
                            int output)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	RefList *reflections;
	int hmax, kmax, lmax;
	double mres;
	signed int h, k, l;

	reflections = reflist_new();

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	/* FIXME: Get this from image */
	mres = 1.0 / 8.0e-10;  /* 8 Angstroms */
	hmax = mres / modulus(asx, asy, asz);
	kmax = mres / modulus(bsx, bsy, bsz);
	lmax = mres / modulus(csx, csy, csz);

	for ( h=-hmax; h<hmax; h++ ) {
	for ( k=-kmax; k<kmax; k++ ) {
	for ( l=-lmax; l<lmax; l++ ) {
		check_reflection(image, mres, output, reflections, h, k, l,
		                 asx,asy,asz,bsx,bsy,bsz,csx,csy,csz);
	}
	}
	}

	return reflections;
}


double integrate_all(struct image *image, RefList *reflections)
{
	double itot = 0.0;
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		float x, y, intensity;
		double xp, yp;
		get_detector_pos(refl, &xp, &yp);

		if ( integrate_peak(image, xp, yp, &x, &y,
                                    &intensity, NULL, NULL, 0, 0) ) continue;

		itot += intensity;
	}

	return itot;
}
