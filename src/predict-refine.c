/*
 * predict-refine.c
 *
 * Prediction refinement
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

#include "image.h"
#include "post-refinement.h"
#include "cell-utils.h"


/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (10)

/* Weighting of excitation error term (m^-1) compared to position term (m) */
#define EXC_WEIGHT (4e-20)

struct reflpeak {
	Reflection *refl;
	struct imagefeature *peak;
	double Ih;   /* normalised */
	struct panel *panel;  /* panel the reflection appears on
                               * (we assume this never changes) */
};

static int pair_peaks(ImageFeatureList *flist, UnitCell *cell, RefList *reflist,
                      struct reflpeak *rps, struct detector *det)
{
	int i;
	const double min_dist = 0.05;
	int n_acc = 0;
	int n_notintegrated = 0;

	for ( i=0; i<image_feature_count(flist); i++ ) {

		struct imagefeature *f;
		double h, k, l, hd, kd, ld;

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(flist, i);
		if ( f == NULL ) continue;

		double ax, ay, az;
		double bx, by, bz;
		double cx, cy, cz;

		cell_get_cartesian(cell,
		                   &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

		/* Decimal and fractional Miller indices of nearest
		 * reciprocal lattice point */
		hd = f->rx * ax + f->ry * ay + f->rz * az;
		kd = f->rx * bx + f->ry * by + f->rz * bz;
		ld = f->rx * cx + f->ry * cy + f->rz * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		/* Check distance */
		if ( (fabs(h - hd) < min_dist)
		  && (fabs(k - kd) < min_dist)
		  && (fabs(l - ld) < min_dist) )
		{
			Reflection *refl;

			/* Dig out the reflection */
			refl = find_refl(reflist, h, k, l);
			if ( refl == NULL ) {
				n_notintegrated++;
				continue;
			}

			rps[n_acc].refl = refl;
			rps[n_acc].peak = f;
			rps[n_acc].panel = find_panel(det, f->fs, f->ss);
			n_acc++;
		}

	}

	return n_acc;
}


static void twod_mapping(double fs, double ss, double *px, double *py,
                         struct panel *p)
{
	double xs, ys;

	xs = fs*p->fsx + ss*p->ssx;
	ys = fs*p->fsy + ss*p->ssy;

	*px = (xs + p->cnx) / p->res;
	*py = (ys + p->cny) / p->res;
}


static double r_gradient(UnitCell *cell, int k, Reflection *refl,
                         struct image *image)
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

		case REF_ASX :
		return - hs * sin(phi) * cos(azi);

		case REF_BSX :
		return - ks * sin(phi) * cos(azi);

		case REF_CSX :
		return - ls * sin(phi) * cos(azi);

		case REF_ASY :
		return - hs * sin(phi) * sin(azi);

		case REF_BSY :
		return - ks * sin(phi) * sin(azi);

		case REF_CSY :
		return - ls * sin(phi) * sin(azi);

		case REF_ASZ :
		return - hs * cos(phi);

		case REF_BSZ :
		return - ks * cos(phi);

		case REF_CSZ :
		return - ls * cos(phi);

	}

	ERROR("No gradient defined for parameter %i\n", k);
	abort();
}


/* Returns d(xh-xpk)/dP + d(yh-ypk)/dP, where P = any parameter */
static double x_gradient(int param, struct reflpeak *rp, struct detector *det,
                         double lambda, UnitCell *cell)
{
	signed int h, k, l;
	double xpk, ypk, xh, yh;
	double fsh, ssh;
	double tt, clen, azi, azf;

	twod_mapping(rp->peak->fs, rp->peak->ss, &xpk, &ypk, rp->panel);
	get_detector_pos(rp->refl, &fsh, &ssh);
	twod_mapping(fsh, ssh, &xh, &yh, rp->panel);
	get_indices(rp->refl, &h, &k, &l);

	tt = asin(lambda * resolution(cell, h, k, l));
	clen = rp->panel->clen;
	azi = atan2(yh, xh);
	azf = 2.0*cos(azi);  /* FIXME: Why factor of 2? */

	switch ( param ) {

		case REF_ASX :
		return h * lambda * clen / cos(tt);

		case REF_BSX :
		return k * lambda * clen / cos(tt);

		case REF_CSX :
		return l * lambda * clen / cos(tt);

		case REF_ASY :
		return 0.0;

		case REF_BSY :
		return 0.0;

		case REF_CSY :
		return 0.0;

		case REF_ASZ :
		return -h * lambda * clen * azf * sin(tt) / (cos(tt)*cos(tt));

		case REF_BSZ :
		return -k * lambda * clen * azf * sin(tt) / (cos(tt)*cos(tt));

		case REF_CSZ :
		return -l * lambda * clen * azf * sin(tt) / (cos(tt)*cos(tt));

	}

	ERROR("Positional gradient requested for parameter %i?\n", param);
	abort();
}


/* Returns d(yh-ypk)/dP, where P = any parameter */
static double y_gradient(int param, struct reflpeak *rp, struct detector *det,
                         double lambda, UnitCell *cell)
{
	signed int h, k, l;
	double xpk, ypk, xh, yh;
	double fsh, ssh;
	double tt, clen, azi, azf;

	twod_mapping(rp->peak->fs, rp->peak->ss, &xpk, &ypk, rp->panel);
	get_detector_pos(rp->refl, &fsh, &ssh);
	twod_mapping(fsh, ssh, &xh, &yh, rp->panel);
	get_indices(rp->refl, &h, &k, &l);

	tt = asin(lambda * resolution(cell, h, k, l));
	clen = rp->panel->clen;
	azi = atan2(yh, xh);
	azf = 2.0*sin(azi);  /* FIXME: Why factor of 2? */

	switch ( param ) {

		case REF_ASX :
		return 0.0;

		case REF_BSX :
		return 0.0;

		case REF_CSX :
		return 0.0;

		case REF_ASY :
		return h * lambda * clen / cos(tt);

		case REF_BSY :
		return k * lambda * clen / cos(tt);

		case REF_CSY :
		return l * lambda * clen / cos(tt);

		case REF_ASZ :
		return -h * lambda * clen * azf * sin(tt) / (cos(tt)*cos(tt));

		case REF_BSZ :
		return -k * lambda * clen * azf * sin(tt) / (cos(tt)*cos(tt));

		case REF_CSZ :
		return -l * lambda * clen * azf * sin(tt) / (cos(tt)*cos(tt));

	}

	ERROR("Positional gradient requested for parameter %i?\n", param);
	abort();
}


static double r_dev(struct reflpeak *rp)
{
	/* Excitation error term */
	double rlow, rhigh, p;
	get_partial(rp->refl, &rlow, &rhigh, &p);
	return (rlow+rhigh)/2.0;
}


static double x_dev(struct reflpeak *rp, struct detector *det)
{
	/* Peak position term */
	double xpk, ypk, xh, yh;
	double fsh, ssh;
	twod_mapping(rp->peak->fs, rp->peak->ss, &xpk, &ypk, rp->panel);
	get_detector_pos(rp->refl, &fsh, &ssh);
	twod_mapping(fsh, ssh, &xh, &yh, rp->panel);
	return xh-xpk;
}


static double y_dev(struct reflpeak *rp, struct detector *det)
{
	/* Peak position term */
	double xpk, ypk, xh, yh;
	double fsh, ssh;
	twod_mapping(rp->peak->fs, rp->peak->ss, &xpk, &ypk, rp->panel);
	get_detector_pos(rp->refl, &fsh, &ssh);
	twod_mapping(fsh, ssh, &xh, &yh, rp->panel);
	return yh-ypk;
}


static int iterate(struct reflpeak *rps, int n, UnitCell *cell,
                   struct image *image)
{
	int i;
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	/* Number of parameters to refine */
	M = gsl_matrix_calloc(9, 9);
	v = gsl_vector_calloc(9);

	for ( i=0; i<n; i++ ) {

		int k;
		double gradients[9];
		double w;

		/* Excitation error terms */
		w = EXC_WEIGHT * rps[i].Ih;

		for ( k=0; k<9; k++ ) {
			gradients[k] = r_gradient(cell, k, rps[i].refl, image);
		}

		for ( k=0; k<9; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<9; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = w * gradients[g] * gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = w * r_dev(&rps[i]);
			v_c *= -gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		/* Positional x terms */
		for ( k=0; k<9; k++ ) {
			gradients[k] = x_gradient(k, &rps[i], image->det,
			                          image->lambda, cell);
		}

		for ( k=0; k<9; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<9; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = gradients[g] * gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = x_dev(&rps[i], image->det);
			v_c *= -gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		/* Positional y terms */
		for ( k=0; k<9; k++ ) {
			gradients[k] = y_gradient(k, &rps[i], image->det,
			                          image->lambda, cell);
		}

		for ( k=0; k<9; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<9; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = gradients[g] * gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = y_dev(&rps[i], image->det);
			v_c *= -gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

	}

	show_matrix_eqn(M, v);
	shifts = solve_svd(v, M, NULL, 1);

	for ( i=0; i<9; i++ ) {
		STATUS("Shift %i = %e\n", i, gsl_vector_get(shifts, i));
	}

	/* Apply shifts */
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	asx += gsl_vector_get(shifts, 0);
	asy += gsl_vector_get(shifts, 1);
	asz += gsl_vector_get(shifts, 2);
	bsx += gsl_vector_get(shifts, 3);
	bsy += gsl_vector_get(shifts, 4);
	bsz += gsl_vector_get(shifts, 5);
	csx += gsl_vector_get(shifts, 6);
	csy += gsl_vector_get(shifts, 7);
	csz += gsl_vector_get(shifts, 8);
	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	gsl_vector_free(shifts);
	gsl_matrix_free(M);
	gsl_vector_free(v);

	return 0;
}


static double residual(struct reflpeak *rps, int n, struct detector *det)
{
	int i;
	double res = 0.0;
	double r;

	r = 0.0;
	for ( i=0; i<n; i++ ) {
		r += EXC_WEIGHT * rps[i].Ih * pow(r_dev(&rps[i]), 2.0);
	}
	printf("%e ", r);
	res += r;

	r = 0.0;
	for ( i=0; i<n; i++ ) {
		r += pow(x_dev(&rps[i], det), 2.0);
	}
	printf("%e ", r);
	res += r;

	r = 0.0;
	for ( i=0; i<n; i++ ) {
		r += pow(y_dev(&rps[i], det), 2.0);
	}
	printf("%e ", r);
	res += r;

	return res;
}


int refine_prediction(struct image *image, Crystal *cr)
{
	int n;
	int i;
	struct reflpeak *rps;
	double max_I;

	rps = malloc(image_feature_count(image->features)
	                       * sizeof(struct reflpeak));
	if ( rps == NULL ) return 1;

	n = pair_peaks(image->features, crystal_get_cell(cr),
	               crystal_get_reflections(cr), rps, image->det);
	STATUS("%i peaks\n", n);
	if ( n < 10 ) {
		ERROR("Too few paired peaks to refine orientation.\n");
		free(rps);
		return 1;
	}

	/* Normalise the intensities to max 1 */
	max_I = -INFINITY;
	for ( i=0; i<n; i++ ) {
		double cur_I = rps[i].peak->intensity;
		if ( cur_I > max_I ) max_I = cur_I;
	}
	if ( max_I <= 0.0 ) {
		ERROR("All peaks negative?\n");
		free(rps);
		return 1;
	}
	for ( i=0; i<n; i++ ) {
		rps[i].Ih = rps[i].peak->intensity / max_I;
	}

	/* Refine */
	STATUS("Initial residual = %e\n", residual(rps, n, image->det));
	for ( i=0; i<MAX_CYCLES; i++ ) {
		iterate(rps, n, crystal_get_cell(cr), image);
		update_partialities(cr, PMODEL_SCSPHERE);
		STATUS("Residual after iteration %i = %e\n",
		        i, residual(rps, n, image->det));
	}

	free(rps);
	return 0;
}
