/*
 * predict-refine.c
 *
 * Prediction refinement
 *
 * Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2023 Thomas White <taw@physics.org>
 *   2016      Valerio Mariani
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "image.h"
#include "geometry.h"
#include "cell-utils.h"
#include "predict-refine.h"
#include "profile.h"
#include "crystfel-mille.h"


/** \file predict-refine.h */


/* Weighting of excitation error term (m^-1) compared to position term (m) */
#define EXC_WEIGHT (4e-20)


double r_dev(struct reflpeak *rp)
{
	/* Excitation error term */
	return get_exerr(rp->refl);
}


double fs_dev(struct reflpeak *rp, struct detgeom *det)
{
	double fsh, ssh;
	get_detector_pos(rp->refl, &fsh, &ssh);
	return fsh - rp->peak->fs;
}


double ss_dev(struct reflpeak *rp, struct detgeom *det)
{
	double fsh, ssh;
	get_detector_pos(rp->refl, &fsh, &ssh);
	return ssh - rp->peak->ss;
}


double r_gradient(int param, Reflection *refl, UnitCell *cell, double wavelength)
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
	phi = angle_between_2d(tl, zl+1.0/wavelength, 0.0, 1.0); /* 2theta */
	azi = atan2(yl, xl); /* azimuth */

	switch ( param ) {

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

		/* Detector movements don't affect excitation error */
		case GPARAM_DET_TX :
		case GPARAM_DET_TY :
		case GPARAM_DET_TZ :
		case GPARAM_DET_RX :
		case GPARAM_DET_RY :
		case GPARAM_DET_RZ :
		return 0.0;

	}

	ERROR("No r gradient defined for parameter %i\n", param);
	abort();
}


/* Spot position gradients for diffraction physics (anything that changes the
 * diffracted ray direction) */
int fs_ss_gradient_physics(int param, Reflection *refl, UnitCell *cell,
                         struct detgeom_panel *p, gsl_matrix *Minv,
                         double fs, double ss, double mu,
                         float *fsg, float *ssg)
{
	signed int h, k, l;
	double xl, yl, zl, kpred;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	gsl_vector *dRdp;
	gsl_vector *v;

	get_indices(refl, &h, &k, &l);
	kpred = get_kpred(refl);
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	xl = h*asx + k*bsx + l*csx;
	yl = h*asy + k*bsy + l*csy;
	zl = h*asz + k*bsz + l*csz;

	dRdp = gsl_vector_calloc(3);

	switch ( param ) {

		case GPARAM_ASX :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_BSX :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_CSX :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_ASY :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_BSY :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_CSY :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_ASZ :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_BSZ :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		case GPARAM_CSZ :
		gsl_vector_set(dRdp, 0, 0.0);
		gsl_vector_set(dRdp, 1, 0.0);
		gsl_vector_set(dRdp, 2, 0.0);
		break;

		default :
		ERROR("Invalid physics gradient %i\n", param);
		return 1;
	}

	v = gsl_vector_calloc(3);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Minv, dRdp, 0.0, v);

	*fsg = mu*(gsl_vector_get(v, 1) - fs*gsl_vector_get(v, 0));
	*ssg = mu*(gsl_vector_get(v, 2) - ss*gsl_vector_get(v, 0));

	return 0;
}


/* Spot position gradients for panel motions (translation or rotation) */
int fs_ss_gradient_panel(int param, Reflection *refl, UnitCell *cell,
                         struct detgeom_panel *p, gsl_matrix *Minv,
                         double fs, double ss, double mu,
                         gsl_vector *t,
                         float *fsg, float *ssg)
{
	gsl_vector *v;
	gsl_matrix *gM;  /* M^-1 * dM/dx * M^-1 */
	gsl_matrix *dMdp = gsl_matrix_calloc(3, 3);

	switch ( param ) {

		case GPARAM_DET_TX :
		gsl_matrix_set(dMdp, 0, 0, 1.0);
		break;

		case GPARAM_DET_TY :
		gsl_matrix_set(dMdp, 1, 0, 1.0);
		break;

		case GPARAM_DET_TZ :
		gsl_matrix_set(dMdp, 2, 0, 1.0);
		break;

		case GPARAM_DET_RX :
		*fsg = 0.0;
		*ssg = 0.0;
		return 0;

		case GPARAM_DET_RY :
		*fsg = 0.0;
		*ssg = 0.0;
		return 0;

		case GPARAM_DET_RZ :
		*fsg = 0.0;
		*ssg = 0.0;
		return 0;

		default:
		ERROR("Invalid panel gradient %i\n", param);
		return 1;

	}

	gM = matrix_mult3(Minv, dMdp, Minv);
	gsl_matrix_free(dMdp);

	v = gsl_vector_calloc(3);
	gsl_blas_dgemv(CblasNoTrans, -1.0, gM, t, 0.0, v);
	gsl_vector_free(t);
	gsl_matrix_free(gM);

	*fsg = mu*(gsl_vector_get(v, 1) - fs*gsl_vector_get(v, 0));
	*ssg = mu*(gsl_vector_get(v, 2) - ss*gsl_vector_get(v, 0));

	gsl_vector_free(v);

	return 0;
}


/* Returns dfs_refl/dP and fss_refl/dP, where P = any parameter
 * fs_refl is the fast scan position of refl in panel 'p',
 * in metres (not pixels) */
int fs_ss_gradient(int param, Reflection *refl, UnitCell *cell,
                   struct detgeom_panel *p, gsl_matrix *Minv,
                   float *fsg, float *ssg)
{
	signed int h, k, l;
	double xl, yl, zl, kpred;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	double tta, ctt, phi;
	gsl_vector *t;
	gsl_vector *v;
	gsl_matrix *M;
	double mu;
	double fs, ss;

	get_indices(refl, &h, &k, &l);
	kpred = get_kpred(refl);
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	xl = h*asx + k*bsx + l*csx;
	yl = h*asy + k*bsy + l*csy;
	zl = h*asz + k*bsz + l*csz;

	/* Calculate 2theta (scattering angle) and azimuth (phi) */
	tta = atan2(sqrt(xl*xl+yl*yl), kpred+zl);
	ctt = cos(tta);
	phi = atan2(yl, xl);

	/* Set up matrix equation */
	M = gsl_matrix_alloc(3, 3);
	v = gsl_vector_alloc(3);
	t = gsl_vector_alloc(3);
	if ( (M==NULL) || (v==NULL) || (t==NULL) ) {
		ERROR("Failed to allocate vectors for gradient calculation\n");
		return 1;
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
	gsl_matrix_set(M, 2, 0, p->cnz);
	gsl_matrix_set(M, 2, 1, p->fsz);
	gsl_matrix_set(M, 2, 2, p->ssz);

	if ( gsl_linalg_HH_solve(M, t, v) ) {
		ERROR("Failed to solve gradient equation\n");
		return 1;
	}
	gsl_matrix_free(M);

	mu = 1.0 / gsl_vector_get(v, 0);
	fs = mu*gsl_vector_get(v, 1);
	ss = mu*gsl_vector_get(v, 2);
	gsl_vector_free(v);

	if ( param <= GPARAM_CSZ ) {
		gsl_vector_free(t);
		return fs_ss_gradient_physics(param, refl, cell, p,
		                              Minv, fs, ss, mu,
		                              fsg, ssg);
	} else {
		return fs_ss_gradient_panel(param, refl, cell, p,
		                            Minv, fs, ss, mu, t,
		                            fsg, ssg);
	}
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
                                    struct detgeom *det)
{
	int i;

	if ( n < 3 ) return n;

	qsort(rps, n, sizeof(struct reflpeak), cmpd2);

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
	int n_final;
	int n = 0;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double dx, dy;
	RefList *all_reflist;
	double lowest_one_over_d;

	all_reflist = reflist_new();
	cell_get_cartesian(crystal_get_cell(cr),
	                   &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	lowest_one_over_d = lowest_reflection(crystal_get_cell(cr));

	crystal_get_det_shift(cr, &dx, &dy);

	/* First, create a RefList containing the most likely indices for each
	 * peak, with no exclusion criteria */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double h, k, l, hd, kd, ld;
		Reflection *refl;
		double r[3];

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		detgeom_transform_coords(&image->detgeom->panels[f->pn],
		                         f->fs, f->ss, image->lambda,
		                         dx, dy, r);

		/* Decimal and fractional Miller indices of nearest reciprocal
		 * lattice point */
		hd = r[0] * ax + r[1] * ay + r[2] * az;
		kd = r[0] * bx + r[1] * by + r[2] * bz;
		ld = r[0] * cx + r[1] * cy + r[2] * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		/* Don't pair with 000, because that can cause trouble later */
		if ( (h==0) && (k==0) && (l==0) ) continue;

		if ( (fabs(h)>=512) || (fabs(k)>=512) || (fabs(l)>=512) ) {
			ERROR("Peak %i (on panel %s at %.2f,%.2f) has indices too "
			      "large for pairing (%.0f %.0f %.0f)\n",
			      i, image->detgeom->panels[f->pn].name,
			      f->fs, f->ss, h, k, l);
			continue;
		}

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
		 * filled in by update_predictions(). */
		set_panel_number(refl, f->pn);

		rps[n].refl = refl;
		rps[n].peak = f;
		n++;

	}

	/* Get the excitation errors and detector positions for the candidate
	 * reflections */
	crystal_set_reflections(cr, all_reflist);
	update_predictions(cr);

	/* Pass over the peaks again, keeping only the ones which look like
	 * good pairings */
	for ( i=0; i<n; i++ ) {

		double fs, ss;
		signed int h, k, l;
		int pnl;
		double refl_r[3];
		double pk_r[3];
		Reflection *refl = rps[i].refl;

		get_indices(refl, &h, &k, &l);

		/* Is the supposed reflection anywhere near the peak? */
		get_detector_pos(refl, &fs, &ss);

		pnl = get_panel_number(refl);
		detgeom_transform_coords(&image->detgeom->panels[pnl],
		                         fs, ss,
		                         image->lambda, dx, dy, refl_r);
		detgeom_transform_coords(&image->detgeom->panels[pnl],
		                         rps[i].peak->fs, rps[i].peak->ss,
		                         image->lambda, dx, dy, pk_r);

		if ( modulus(refl_r[0] - pk_r[0],
		             refl_r[1] - pk_r[1],
		             refl_r[2] - pk_r[2]) > lowest_one_over_d / 3.0 )
		{
			continue;
		}

		rps[n_acc] = rps[i];
		rps[n_acc].refl = reflection_new(h, k, l);
		copy_data(rps[n_acc].refl, refl);
		n_acc++;

	}
	reflist_free(all_reflist);
	crystal_set_reflections(cr, NULL);

	/* Sort the pairings by excitation error and look for a transition
	 * between good pairings and outliers */
	n_final = check_outlier_transition(rps, n_acc, image->detgeom);

	/* Add the final accepted reflections to the caller's list */
	if ( reflist != NULL ) {
		for ( i=0; i<n_final; i++ ) {
			add_refl_to_list(rps[i].refl, reflist);
		}
	}

	/* Free the reflections beyond the outlier cutoff */
	for ( i=n_final; i<n_acc; i++ ) {
		reflection_free(rps[i].refl);
	}

	return n_final;
}


int refine_radius(Crystal *cr, struct image *image)
{
	int n, n_acc;
	struct reflpeak *rps;
	RefList *reflist;

	/* Maximum possible size */
	rps = malloc(image_feature_count(image->features)
	                  * sizeof(struct reflpeak));
	if ( rps == NULL ) return 1;

	reflist = reflist_new();
	n_acc = pair_peaks(image, cr, reflist, rps);
	if ( n_acc < 3 ) {
		free(rps);
		reflist_free(reflist);
		return 1;
	}
	crystal_set_reflections(cr, reflist);
	update_predictions(cr);
	crystal_set_reflections(cr, NULL);

	qsort(rps, n_acc, sizeof(struct reflpeak), cmpd2);
	n = (n_acc-1) - n_acc/50;
	if ( n < 2 ) n = 2; /* n_acc is always >= 2 */
	crystal_set_profile_radius(cr, fabs(r_dev(&rps[n])));

	reflist_free(reflist);
	free(rps);

	return 0;
}


static int iterate(struct reflpeak *rps, int n, UnitCell *cell,
                   struct image *image, gsl_matrix **Minvs)
{
	int i;
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	const enum gparam rv[] =
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
	};
	int num_params = 9;

	/* Number of parameters to refine */
	M = gsl_matrix_calloc(num_params, num_params);
	v = gsl_vector_calloc(num_params);

	for ( i=0; i<n; i++ ) {

		int k;
		double r_gradients[num_params];
		float fs_gradients[num_params];
		float ss_gradients[num_params];
		double w;

		w = rps[i].Ih;

		/* Calculate the gradients w.r.t. all parameters, for the
		 * three measurements (r, fs, ss) of this reflpeak */
		for ( k=0; k<num_params; k++ ) {
			r_gradients[k] = r_gradient(rv[k], rps[i].refl, cell,
			                            image->lambda);
			fs_ss_gradient(rv[k], rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn],
			               &fs_gradients[k], &ss_gradients[k]);

		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = w * r_gradients[g] * r_gradients[k];
				M_c += w * fs_gradients[g] * fs_gradients[k];
				M_c += w * ss_gradients[g] * ss_gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = -w * r_dev(&rps[i]) * r_gradients[k];
			v_c += -fs_dev(&rps[i], image->detgeom) * fs_gradients[k];
			v_c += -ss_dev(&rps[i], image->detgeom) * ss_gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

	}

	int k;
	for ( k=0; k<num_params; k++ ) {
		double M_curr;
		M_curr = gsl_matrix_get(M, k, k);
		gsl_matrix_set(M, k, k, M_curr);
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

	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	gsl_vector_free(shifts);
	gsl_matrix_free(M);
	gsl_vector_free(v);

	return 0;
}


static double pred_residual(struct reflpeak *rps, int n, struct detgeom *det)
{
	int i;
	double res = 0.0;
	double r;

	r = 0.0;
	for ( i=0; i<n; i++ ) {
		r += EXC_WEIGHT * rps[i].Ih * pow(r_dev(&rps[i]), 2.0);
	}
	res += r;

	r = 0.0;
	for ( i=0; i<n; i++ ) {
		r += pow(det->panels[rps[i].peak->pn].pixel_pitch*fs_dev(&rps[i], det), 2.0);
	}
	res += r;

	r = 0.0;
	for ( i=0; i<n; i++ ) {
		r += pow(det->panels[rps[i].peak->pn].pixel_pitch*ss_dev(&rps[i], det), 2.0);
	}
	res += r;

	return res;
}


/* NB Only for use when the list of reflpeaks was created without a RefList.
 * If a RefList was used, then reflist_free the list then just free() the rps */
static void free_rps_noreflist(struct reflpeak *rps, int n)
{
	int i;

	for ( i=0; i<n; i++ ) {
		reflection_free(rps[i].refl);
	}
	free(rps);
}


int refine_prediction(struct image *image, Crystal *cr, Mille *mille)
{
	int n;
	int i;
	struct reflpeak *rps;
	double max_I;
	RefList *reflist;
	char tmp[256];
	gsl_matrix **Minvs;

	rps = malloc(image_feature_count(image->features)
	                       * sizeof(struct reflpeak));
	if ( rps == NULL ) return 1;

	reflist = reflist_new();
	n = pair_peaks(image, cr, reflist, rps);
	if ( n < 10 ) {
		free(rps);
		reflist_free(reflist);
		return 1;
	}
	crystal_set_reflections(cr, reflist);

	Minvs = make_panel_minvs(image->detgeom);

	/* Normalise the intensities to max 1 */
	max_I = -INFINITY;
	for ( i=0; i<n; i++ ) {
		double cur_I = rps[i].peak->intensity;
		if ( cur_I > max_I ) max_I = cur_I;
	}
	if ( max_I <= 0.0 ) {
		ERROR("All peaks negative?\n");
		free(rps);
		crystal_set_reflections(cr, NULL);
		return 1;
	}
	for ( i=0; i<n; i++ ) {
		if ( rps[i].peak->intensity > 0.0 ) {
			rps[i].Ih = rps[i].peak->intensity / max_I;
		} else {
			rps[i].Ih = 0.0;
		}
	}

	//STATUS("Initial residual = %e\n",
	//       pred_residual(rps, n, image->detgeom));

	/* Refine (max 10 cycles) */
	for ( i=0; i<10; i++ ) {
		update_predictions(cr);
		if ( iterate(rps, n, crystal_get_cell(cr), image, Minvs) )
		{
			crystal_set_reflections(cr, NULL);
			return 1;
		}
		//STATUS("Residual after %i = %e\n", i,
		//       pred_residual(rps, n, image->detgeom));
	}
	//STATUS("Final residual = %e\n",
	//       pred_residual(rps, n, image->detgeom));

	snprintf(tmp, 255, "predict_refine/final_residual = %e",
	         pred_residual(rps, n, image->detgeom));
	crystal_add_notes(cr, tmp);

	#ifdef HAVE_MILLEPEDE
	if ( mille != NULL ) {
		profile_start("mille-calc");
		write_mille(mille, n, crystal_get_cell(cr), rps, image, Minvs);
		profile_end("mille-calc");
	}
	#endif

	for ( i=0; i<image->detgeom->n_panels; i++ ) {
		gsl_matrix_free(Minvs[i]);
	}
	free(Minvs);

	crystal_set_reflections(cr, NULL);
	reflist_free(reflist);

	n = pair_peaks(image, cr, NULL, rps);
	free_rps_noreflist(rps, n);
	if ( n < 10 ) {
		#ifdef HAVE_MILLEPEDE
		if ( mille != NULL ) {
			crystfel_mille_delete_last_record(mille);
		}
		#endif
		return 1;
	}

	#ifdef HAVE_MILLEPEDE
	if ( mille != NULL ) {
		profile_start("mille-write");
		crystfel_mille_write_record(mille);
		profile_end("mille-write");
	}
	#endif

	return 0;
}
