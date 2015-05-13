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
#include "geometry.h"
#include "cell-utils.h"


/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (10)

/* Weighting of excitation error term (m^-1) compared to position term (m) */
#define EXC_WEIGHT (4e-20)

/* Parameters to refine */
static const enum gparam rv[] =
{
	GPARAM_ASX,
	GPARAM_ASY,
	GPARAM_ASZ,
	GPARAM_BSX,
	GPARAM_BSY,
	GPARAM_BSZ,
	GPARAM_CSX,
	GPARAM_CSY,
	GPARAM_CSZ,
	GPARAM_DETX,
	GPARAM_DETY,
};

static const int num_params = 11;

struct reflpeak {
	Reflection *refl;
	struct imagefeature *peak;
	double Ih;   /* normalised */
	struct panel *panel;  /* panel the reflection appears on
                               * (we assume this never changes) */
};


static void twod_mapping(double fs, double ss, double *px, double *py,
                         struct panel *p)
{
	double xs, ys;

	xs = fs*p->fsx + ss*p->ssx;
	ys = fs*p->fsy + ss*p->ssy;

	*px = (xs + p->cnx) / p->res;
	*py = (ys + p->cny) / p->res;
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


static void UNUSED write_pairs(const char *filename, struct reflpeak *rps,
                               int n, struct detector *det)
{
	int i;
	FILE *fh;

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", filename);
		return;
	}

	for ( i=0; i<n; i++ ) {

		double write_fs, write_ss;
		double fs, ss;
		struct panel *p;

		fs = rps[i].peak->fs;
		ss = rps[i].peak->ss;

		p = find_panel(det, fs, ss);
		write_fs = fs - p->min_fs + p->orig_min_fs;
		write_ss = ss - p->min_ss + p->orig_min_ss;

		fprintf(fh, "%7.2f %7.2f dev r,x,y: %9f %9f %9f %9f\n",
		        write_fs, write_ss,
		        r_dev(&rps[i])/1e9, fabs(r_dev(&rps[i])/1e9),
		        x_dev(&rps[i], det),
		        y_dev(&rps[i], det));

	}

	fclose(fh);

	STATUS("Wrote %i pairs to %s\n", n, filename);
}


static int cmpd2(const void *av, const void *bv)
{
	struct reflpeak *a, *b;

	a = (struct reflpeak *)av;
	b = (struct reflpeak *)bv;

	if ( fabs(r_dev(a)) < fabs(r_dev(b)) ) return -1;
	return 1;
}


static int check_outlier_transition(struct reflpeak *rps, int n,
                                    struct detector *det)
{
	int i;

	if ( n < 3 ) return n;

	qsort(rps, n, sizeof(struct reflpeak), cmpd2);
	//write_pairs("pairs-before-outlier.lst", rps, n, det);

	for ( i=1; i<n-1; i++ ) {

		int j;
		double grad = fabs(r_dev(&rps[i])) / i;

		for ( j=i+1; j<n; j++ ) {
			if ( fabs(r_dev(&rps[j])) < 0.001e9+grad*j ) {
				break;
			}
		}
		if ( j == n ) {
			//STATUS("Outlier transition found at position %i / %i\n",
			//       i, n);
			return i;
		}
	}

	//STATUS("No outlier transition found.\n");
	return n;
}


/* Associate a Reflection with each peak in "image" which is close to Bragg.
 * Reflections will be added to "reflist", which can be NULL if this is not
 * needed.  "rps" must be an array of sufficient size for all the peaks */
static int pair_peaks(struct image *image, Crystal *cr,
                      RefList *reflist, struct reflpeak *rps)
{
	int i;
	int n_acc = 0;
	int n = 0;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	RefList *all_reflist;

	all_reflist = reflist_new();
	cell_get_cartesian(crystal_get_cell(cr),
	                   &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* First, create a RefList containing the most likely indices for each
	 * peak, with no exclusion criteria */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double h, k, l, hd, kd, ld;
		Reflection *refl;
		struct panel *p;

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		/* Decimal and fractional Miller indices of nearest reciprocal
		 * lattice point */
		hd = f->rx * ax + f->ry * ay + f->rz * az;
		kd = f->rx * bx + f->ry * by + f->rz * bz;
		ld = f->rx * cx + f->ry * cy + f->rz * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		refl = reflection_new(h, k, l);
		if ( refl == NULL ) {
			ERROR("Failed to create reflection\n");
			return 0;
		}

		add_refl_to_list(refl, all_reflist);
		set_symmetric_indices(refl, h, k, l);

		/* It doesn't matter if the actual predicted location
		 * doesn't fall on this panel.  We're only interested
		 * in how far away it is from the peak location.
		 * The predicted position and excitation errors will be
		 * filled in by update_partialities(). */
		p = find_panel(image->det, f->fs, f->ss);
		set_panel(refl, p);

		rps[n].refl = refl;
		rps[n].peak = f;
		rps[n].panel = p;
		n++;

	}

	/* Get the excitation errors and detector positions for the candidate
	 * reflections */
	crystal_set_reflections(cr, all_reflist);
	update_partialities(cr, PMODEL_SCSPHERE);

	/* Pass over the peaks again, keeping only the ones which look like
	 * good pairings */
	for ( i=0; i<n; i++ ) {

		double fs, ss, pd;
		signed int h, k, l;
		Reflection *refl = rps[i].refl;

		get_indices(refl, &h, &k, &l);

		/* Is the supposed reflection anywhere near the peak? */
		get_detector_pos(refl, &fs, &ss);
		pd = pow(fs - rps[i].peak->fs, 2.0)
		   + pow(ss - rps[i].peak->ss, 2.0);
		if ( pd > 10.0 * 10.0 ) continue;

		rps[n_acc] = rps[i];
		rps[n_acc].refl = reflection_new(h, k, l);
		copy_data(rps[n_acc].refl, refl);
		if ( reflist != NULL ) {
			add_refl_to_list(rps[n_acc].refl, reflist);
		}
		n_acc++;

	}
	reflist_free(all_reflist);

	/* Sort the pairings by excitation error and look for a transition
	 * between good pairings and outliers */
	n_acc = check_outlier_transition(rps, n_acc, image->det);

	return n_acc;
}


void refine_radius(Crystal *cr, struct image *image)
{
	int n, n_acc;
	struct reflpeak *rps;
	RefList *reflist;

	/* Maximum possible size */
	rps = malloc(image_feature_count(image->features)
	                  * sizeof(struct reflpeak));
	if ( rps == NULL ) return;

	reflist = reflist_new();
	n_acc = pair_peaks(image, cr, reflist, rps);
	if ( n_acc < 3 ) {
		ERROR("Too few paired peaks (%i) to determine radius\n", n_acc);
		free(rps);
		return;
	}
	crystal_set_reflections(cr, reflist);
	update_partialities(cr, PMODEL_SCSPHERE);
	crystal_set_reflections(cr, NULL);

	qsort(rps, n_acc, sizeof(struct reflpeak), cmpd2);
	n = (n_acc-1) - n_acc/50;
	if ( n < 2 ) n = 2; /* n_acc is always >= 2 */
	crystal_set_profile_radius(cr, fabs(r_dev(&rps[n])));

	reflist_free(reflist);
	free(rps);
}


/* Returns d(xh-xpk)/dP, where P = any parameter */
static double x_gradient(int param, struct reflpeak *rp, struct detector *det,
                         double lambda, UnitCell *cell)
{
	signed int h, k, l;
	double xpk, ypk, xh, yh;
	double fsh, ssh;
	double x, z, wn;
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	twod_mapping(rp->peak->fs, rp->peak->ss, &xpk, &ypk, rp->panel);
	get_detector_pos(rp->refl, &fsh, &ssh);
	twod_mapping(fsh, ssh, &xh, &yh, rp->panel);
	get_indices(rp->refl, &h, &k, &l);

	wn = 1.0 / lambda;

	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	x = h*ax + k*bx + l*cx;
	z = h*az + k*bz + l*cz;

	switch ( param ) {

		case GPARAM_ASX :
		return h * rp->panel->clen / (wn+z);

		case GPARAM_BSX :
		return k * rp->panel->clen / (wn+z);

		case GPARAM_CSX :
		return l * rp->panel->clen / (wn+z);

		case GPARAM_ASY :
		return 0.0;

		case GPARAM_BSY :
		return 0.0;

		case GPARAM_CSY :
		return 0.0;

		case GPARAM_ASZ :
		return -h * x * rp->panel->clen / (wn*wn + 2*wn*z + z*z);

		case GPARAM_BSZ :
		return -k * x * rp->panel->clen / (wn*wn + 2*wn*z + z*z);

		case GPARAM_CSZ :
		return -l * x * rp->panel->clen / (wn*wn + 2*wn*z + z*z);

		case GPARAM_DETX :
		return -1;

		case GPARAM_DETY :
		return 0;

		case GPARAM_CLEN :
		return x / (wn+z);

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
	double y, z, wn;
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	twod_mapping(rp->peak->fs, rp->peak->ss, &xpk, &ypk, rp->panel);
	get_detector_pos(rp->refl, &fsh, &ssh);
	twod_mapping(fsh, ssh, &xh, &yh, rp->panel);
	get_indices(rp->refl, &h, &k, &l);

	wn = 1.0 / lambda;

	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	y = h*ay + k*by + l*cy;
	z = h*az + k*bz + l*cz;

	switch ( param ) {

		case GPARAM_ASX :
		return 0.0;

		case GPARAM_BSX :
		return 0.0;

		case GPARAM_CSX :
		return 0.0;

		case GPARAM_ASY :
		return h * rp->panel->clen / (wn+z);

		case GPARAM_BSY :
		return k * rp->panel->clen / (wn+z);

		case GPARAM_CSY :
		return l * rp->panel->clen / (wn+z);

		case GPARAM_ASZ :
		return -h * y * rp->panel->clen / (wn*wn + 2*wn*z + z*z);

		case GPARAM_BSZ :
		return -k * y * rp->panel->clen / (wn*wn + 2*wn*z + z*z);

		case GPARAM_CSZ :
		return -l * y * rp->panel->clen / (wn*wn + 2*wn*z + z*z);

		case GPARAM_DETX :
		return 0;

		case GPARAM_DETY :
		return -1;

		case GPARAM_CLEN :
		return y / (wn+z);

	}

	ERROR("Positional gradient requested for parameter %i?\n", param);
	abort();
}


static void update_detector(struct detector *det, double xoffs, double yoffs,
                            double coffs)
{
	int i;

	for ( i=0; i<det->n_panels; i++ ) {
		struct panel *p = &det->panels[i];
		p->cnx += xoffs * p->res;
		p->cny += yoffs * p->res;
		p->clen += coffs;
	}
}


static int iterate(struct reflpeak *rps, int n, UnitCell *cell,
                   struct image *image,
                   double *total_x, double *total_y, double *total_z)
{
	int i;
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	/* Number of parameters to refine */
	M = gsl_matrix_calloc(num_params, num_params);
	v = gsl_vector_calloc(num_params);

	for ( i=0; i<n; i++ ) {

		int k;
		double gradients[num_params];
		double w;

		/* Excitation error terms */
		w = EXC_WEIGHT * rps[i].Ih;

		for ( k=0; k<num_params; k++ ) {
			gradients[k] = r_gradient(cell, rv[k], rps[i].refl,
			                          image);
		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

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
		for ( k=0; k<num_params; k++ ) {
			gradients[k] = x_gradient(rv[k], &rps[i], image->det,
			                          image->lambda, cell);
		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

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
		for ( k=0; k<num_params; k++ ) {
			gradients[k] = y_gradient(rv[k], &rps[i], image->det,
			                          image->lambda, cell);
		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

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

	//show_matrix_eqn(M, v);
	shifts = solve_svd(v, M, NULL, 0);
	if ( shifts == NULL ) {
		ERROR("Failed to solve equations.\n");
		gsl_matrix_free(M);
		gsl_vector_free(v);
		return 1;
	}

	for ( i=0; i<num_params; i++ ) {
	//	STATUS("Shift %i = %e\n", i, gsl_vector_get(shifts, i));
		if ( isnan(gsl_vector_get(shifts, i)) ) {
			gsl_vector_set(shifts, i, 0.0);
		}
	}

	/* Apply shifts */
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	/* Ensure the order here matches the order in rv[] */
	asx += gsl_vector_get(shifts, 0);
	asy += gsl_vector_get(shifts, 1);
	asz += gsl_vector_get(shifts, 2);
	bsx += gsl_vector_get(shifts, 3);
	bsy += gsl_vector_get(shifts, 4);
	bsz += gsl_vector_get(shifts, 5);
	csx += gsl_vector_get(shifts, 6);
	csy += gsl_vector_get(shifts, 7);
	csz += gsl_vector_get(shifts, 8);
	update_detector(image->det, gsl_vector_get(shifts, 9),
	                            gsl_vector_get(shifts, 10),
	                            gsl_vector_get(shifts, 11));
	*total_x += gsl_vector_get(shifts, 9);
	*total_y += gsl_vector_get(shifts, 10);
	*total_z += gsl_vector_get(shifts, 11);

	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	gsl_vector_free(shifts);
	gsl_matrix_free(M);
	gsl_vector_free(v);

	return 0;
}


static double UNUSED residual(struct reflpeak *rps, int n, struct detector *det)
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
	printf("%e\n", r);
	res += r;

	return res;
}


int refine_prediction(struct image *image, Crystal *cr)
{
	int n;
	int i;
	struct reflpeak *rps;
	double max_I;
	RefList *reflist;
	double total_x = 0.0;
	double total_y = 0.0;
	double total_z = 0.0;
	char tmp[1024];

	rps = malloc(image_feature_count(image->features)
	                       * sizeof(struct reflpeak));
	if ( rps == NULL ) return 1;

	reflist = reflist_new();
	n = pair_peaks(image, cr, reflist, rps);
	if ( n < 10 ) {
		ERROR("Too few paired peaks (%i) to refine orientation.\n", n);
		free(rps);
		reflist_free(reflist);
		return 1;
	}
	crystal_set_reflections(cr, reflist);

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

	//STATUS("Initial residual = %e\n", residual(rps, n, image->det));

	/* Refine */
	for ( i=0; i<MAX_CYCLES; i++ ) {
		update_partialities(cr, PMODEL_SCSPHERE);
		if ( iterate(rps, n, crystal_get_cell(cr), image,
		             &total_x, &total_y, &total_z) ) return 1;
		//STATUS("Residual after %i = %e\n", i,
		//       residual(rps, n, image->det));
	}
	//STATUS("Final residual = %e\n", residual(rps, n, image->det));

	snprintf(tmp, 1024, "predict_refine/det_shift x = %.3f y = %.3f mm",
	                    total_x*1e3, total_y*1e3);
	crystal_add_notes(cr, tmp);

	crystal_set_reflections(cr, NULL);
	reflist_free(reflist);

	n = pair_peaks(image, cr, NULL, rps);
	free(rps);
	if ( n < 10 ) {
		ERROR("Too few paired peaks (%i) after refinement.\n", n);
		return 1;
	}

	return 0;
}
