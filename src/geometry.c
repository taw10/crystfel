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
#include <gsl/gsl_poly.h>
#include <assert.h>

#include "utils.h"
#include "cell.h"
#include "image.h"
#include "peaks.h"
#include "beam-parameters.h"


#define MAX_CPEAKS (256 * 256)


static signed int locate_peak(double x, double y, double z, double lambda,
                              struct detector *det, double *xdap, double *ydap)
{
	int p;
	signed int found = -1;
	const double den = 1.0/lambda + z;

	for ( p=0; p<det->n_panels; p++ ) {

		double xd, yd, cl;
		double xda, yda;

		/* Camera length for this panel */
		cl = det->panels[p].clen;

		/* Coordinates of peak relative to central beam, in m */
		xd = cl * x / den;
		yd = cl * y / den;

		/* Convert to pixels */
		xd *= det->panels[p].res;
		yd *= det->panels[p].res;

		/* Add the coordinates of the central beam */
		xda = xd + det->panels[p].cx;
		yda = yd + det->panels[p].cy;

		/* Now, is this on this panel? */
		if ( xda < det->panels[p].min_x ) continue;
		if ( xda > det->panels[p].max_x ) continue;
		if ( yda < det->panels[p].min_y ) continue;
		if ( yda > det->panels[p].max_y ) continue;

		/* If peak appears on multiple panels, reject it */
		if ( found != -1 ) return -1;

		/* Woohoo! */
		found = p;
		*xdap = xda;
		*ydap = yda;

	}

	return found;
}


static double excitation_error(double xl, double yl, double zl,
                               double ds_sq, double lambda)
{
	double tt;
	double a, b, c;
	double r1, r2;
	int n;

	tt = angle_between(0.0, 0.0, 1.0,  xl, yl, zl+1.0/lambda);

	a = 1.0;
	b = - cos(tt) * 2.0/lambda;
	c = pow(lambda, -2.0) - ds_sq;

	/* FIXME: I don't trust GSL's quadratic solver */
	n = gsl_poly_solve_quadratic(a, b, c, &r1, &r2);
	assert(n == 2);

	r1 -= 1.0/lambda;
	r2 -= 1.0/lambda;

	if ( r1 > r2 ) return r2;
	return r1;
}


struct cpeak *find_intersections(struct image *image, UnitCell *cell,
                                 int *n, int output)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct cpeak *cpeaks;
	int np = 0;
	int hmax, kmax, lmax;
	double mres;
	signed int h, k, l;
	double bandwidth = image->bw;
	double divergence = image->div;
	double lambda = image->lambda;
	const double profile_cutoff = 0.05e9;  /* 0.1 nm^-1 */
	double llow, lhigh;    /* Wavelength */

	cpeaks = malloc(sizeof(struct cpeak)*MAX_CPEAKS);
	if ( cpeaks == NULL ) {
		*n = 0;
		return NULL;
	}

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	/* FIXME: Account for left-handed indexing */
	asz = -asz;  bsz = -bsz;  csz = -csz;

	mres = 1.0 / 8.0e-10;  /* 8 Angstroms */
	hmax = mres / modulus(asx, asy, asz);
	kmax = mres / modulus(bsx, bsy, bsz);
	lmax = mres / modulus(csx, csy, csz);

	/* "low" gives the largest Ewald sphere,
	 * "high" gives the smallest Ewald sphere. */
	llow = lambda - lambda*bandwidth/2.0;
	lhigh = lambda + lambda*bandwidth/2.0;

	for ( h=-hmax; h<hmax; h++ ) {
	for ( k=-kmax; k<kmax; k++ ) {
	for ( l=-lmax; l<lmax; l++ ) {

		double xl, yl, zl;
		double ds, ds_sq;
		double rlow, rhigh;    /* "Excitation error" */
		signed int plow, phigh;
		double xdalow, xdahigh;
		double ydalow, ydahigh;
		double xda, yda;

		/* Ignore central beam */
		if ( (h==0) && (k==0) && (l==0) ) continue;

		/* Get the coordinates of the reciprocal lattice point */
		zl = h*asz + k*bsz + l*csz;
		if ( zl > 0.0 ) continue;  /* Throw out if it's "in front" */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;

		/* Calculate reciprocal lattice point modulus (and square) */
		ds_sq = modulus_squared(xl, yl, zl);  /* d*^2 */
		ds = sqrt(ds_sq);
		if ( ds > mres ) continue;  /* Outside resolution range */

		/* Calculate excitation errors */
		rlow = excitation_error(xl, yl, zl, ds_sq, llow);
		rhigh = excitation_error(xl, yl, zl, ds_sq, lhigh);

		if ( (fabs(rlow) > profile_cutoff)
		  && (fabs(rhigh) > profile_cutoff) ) {
			/* Reflection is not close to Bragg at either of
			 * the two extremes of wavelength, so skip it. */
			continue;
		}

		/* Locate peak on detector, and check it doesn't span panels */
		plow = locate_peak(xl, yl, zl, llow, image->det,
		                   &xdalow, &ydalow);
		if ( plow == -1 ) continue;
		phigh = locate_peak(xl, yl, zl, lhigh, image->det,
		                    &xdahigh, &ydahigh);
		if ( phigh == -1 ) continue;
		if ( plow != phigh ) continue;

		xda = xdahigh;
		yda = ydahigh;
		cpeaks[np].h = h;
		cpeaks[np].k = k;
		cpeaks[np].l = l;
		cpeaks[np].x = xdahigh;
		cpeaks[np].y = ydahigh;
		np++;

		if ( output ) {
			printf("%3i %3i %3i %6f (at %5.2f,%5.2f) %9e %9e\n",
			       h, k, l, 0.0, xda, yda, rlow, rhigh);
		}

		if ( np == MAX_CPEAKS ) goto out;

	}
	}
	}

out:
	*n = np;
	return cpeaks;
}


double integrate_all(struct image *image, struct cpeak *cpeaks, int n)
{
	double itot = 0.0;
	int i;

	for ( i=0; i<n; i++ ) {

		float x, y, intensity;

		if ( integrate_peak(image, cpeaks[i].x, cpeaks[i].y, &x, &y,
                                    &intensity, NULL, NULL, 0, 0) ) continue;

		itot += intensity;
	}

	return itot;
}
