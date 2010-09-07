/*
 * geometry.h
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

#include "utils.h"
#include "cell.h"
#include "image.h"
#include "peaks.h"


#define MAX_HITS (1024)


struct reflhit *find_intersections(struct image *image, UnitCell *cell,
                                   double divergence, double bandwidth,
                                   int *n, int output)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct reflhit *hits;
	int np = 0;
	int hmax, kmax, lmax;
	double mres;
	signed int h, k, l;

	hits = malloc(sizeof(struct reflhit)*MAX_HITS);
	if ( hits == NULL ) {
		*n = 0;
		return NULL;
	}

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	mres = 1.0 / 8.0e-10;  /* 8 Angstroms */
	hmax = mres / modulus(asx, asy, asz);
	kmax = mres / modulus(bsx, bsy, bsz);
	lmax = mres / modulus(csx, csy, csz);

	for ( h=-hmax; h<hmax; h++ ) {
	for ( k=-kmax; k<kmax; k++ ) {
	for ( l=-lmax; l<lmax; l++ ) {

		double xl, yl, zl;
		double ds_sq, dps_sq;
		double delta, divfact;
		double llow, lhigh;
		double xd, yd, cl;
		double xda, yda;
		int p;
		int found = 0;

		if ( (h==0) && (k==0) && (l==0) ) continue;

		llow = image->lambda - image->lambda*bandwidth/2.0;
		lhigh = image->lambda + image->lambda*bandwidth/2.0;

		/* Get the coordinates of the reciprocal lattice point */
		zl = h*asz + k*bsz + l*csz;
		if ( zl < 0.0 ) continue;  /* Do this check very early */
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;

		ds_sq = modulus_squared(xl, yl, zl);  /* d*^2 */
		delta = divergence/image->lambda;
		dps_sq = ds_sq + pow(delta, 2.0);  /* d'*^2 */

		/* In range? */
		divfact = 2.0 * delta * sqrt(xl*xl + yl*yl);
		if ( ds_sq - 2.0*zl/llow > 0.0 ) continue;
		if ( ds_sq - 2.0*zl/lhigh < 0.0 ) continue;

		/* Work out which panel this peak would fall on */
		for ( p=0; p<image->det->n_panels; p++ ) {

			/* Camera length for this panel */
			cl = image->det->panels[p].clen;

			/* Coordinates of peak relative to central beam, in m */
			xd = cl*xl / (ds_sq/(2.0*zl) - zl);
			yd = cl*yl / (ds_sq/(2.0*zl) - zl);

			/* Convert to pixels */
			xd *= image->det->panels[p].res;
			yd *= image->det->panels[p].res;

			/* Add the coordinates of the central beam */
			xda = xd + image->det->panels[p].cx;
			yda = yd + image->det->panels[p].cy;

			/* Now, is this on this panel? */
			if ( xda < image->det->panels[p].min_x ) continue;
			if ( xda > image->det->panels[p].max_x ) continue;
			if ( yda < image->det->panels[p].min_y ) continue;
			if ( yda > image->det->panels[p].max_y ) continue;

			/* Woohoo! */
			found = 1;
			break;

		}

		if ( !found ) continue;

		hits[np].h = h;
		hits[np].k = k;
		hits[np].l = l;
		hits[np].x = xda;
		hits[np].y = yda;
		np++;

		if ( output ) {
			printf("%i %i %i 0.0 (at %f,%f)\n", h, k, l, xda, yda);
		}

	}
	}
	}

	*n = np;
	return hits;
}


double integrate_all(struct image *image, struct reflhit *hits, int n)
{
	double itot = 0.0;
	int i;

	for ( i=0; i<n; i++ ) {

		float x, y, intensity;

		if ( integrate_peak(image, hits[i].x, hits[i].y, &x, &y,
                                    &intensity, 0, 0) ) continue;

		itot += intensity;
	}

	return itot;
}
