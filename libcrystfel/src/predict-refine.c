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

double r_dev(struct reflpeak *rp)
{
	/* Excitation error term */
	return EXC_WEIGHT * get_exerr(rp->refl);
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


void crossp_norm(double c1[3], double c2[3], double u[3])
{
	double nrm;

	u[0] = c1[1] * c2[2] - c1[2] * c2[1];
	u[1] = - c1[0] * c2[2] + c1[2] * c2[0];
	u[2] = c1[0] * c2[1] - c1[1] * c2[0];
	nrm = modulus(u[0], u[1], u[2]);
	u[0] /= nrm;
	u[1] /= nrm;
	u[2] /= nrm;
}


void rotate3d(double vec[3], double axis[3], double ang)
{
	double t = 1.0 - cos(ang);
	double c = cos(ang);
	double s = sin(ang);
	double nx, ny, nz;
	double ux, uy, uz;
	double x, y, z;

	x = vec[0];  y = vec[1];  z = vec[2];
	ux = axis[0]; uy = axis[1]; uz = axis[2];

	nx = (t*ux*ux+c)*x    + (t*ux*uy-s*uz)*y + (t*ux*uz+s*uy)*z;
	ny = (t*ux*uy+s*uz)*x + (t*uy*uy+c)*y    + (t*uy*uz-s*ux)*z;
	nz = (t*ux*uz-s*uy)*x + (t*uz*uy+s*ux)*y + (t*uz*uz+c)*z;

	vec[0] = nx; vec[1] = ny; vec[2] = nz;
}


static gsl_vector *vec3(double x, double y, double z)
{
	gsl_vector *v = gsl_vector_alloc(3);
	gsl_vector_set(v, 0, x);
	gsl_vector_set(v, 1, y);
	gsl_vector_set(v, 2, z);
	return v;
}


/* Calculate dR/d(param), where R is the vector from the center of the Ewald
 * sphere through the reciprocal lattice point */
static gsl_vector *ray_vector_gradient(int param, UnitCell *cell, Reflection *refl,
                                       gsl_vector **q)
{
	double as[3], bs[3], cs[3];
	double u[3];
	double xl, yl, zl;
	signed int h, k, l;
	double m;

	get_symmetric_indices(refl, &h, &k, &l);

	cell_get_reciprocal(cell, &as[0], &as[1], &as[2],
	                          &bs[0], &bs[1], &bs[2],
	                          &cs[0], &cs[1], &cs[2]);

	xl = h*as[0] + k*bs[0] + l*cs[0];
	yl = h*as[1] + k*bs[1] + l*cs[1];
	zl = h*as[2] + k*bs[2] + l*cs[2];

	if ( q != NULL ) {
		*q = vec3(xl, yl, zl);
	}

	switch ( param ) {

		case GPARAM_A_STAR :
		m = modulus(as[0], as[1], as[2]);
		return vec3(h*as[0]/m, h*as[1]/m, h*as[2]/m);

		case GPARAM_B_STAR :
		m = modulus(bs[0], bs[1], bs[2]);
		return vec3(k*bs[0]/m, k*bs[1]/m, k*bs[2]/m);

		case GPARAM_C_STAR :
		m = modulus(cs[0], cs[1], cs[2]);
		return vec3(l*cs[0]/m, l*cs[1]/m, l*cs[2]/m);

		case GPARAM_AL_STAR :
		crossp_norm(cs, bs, u);
		return vec3(u[1]*bs[2]*k - u[2]*bs[1]*k,
		            u[2]*bs[0]*k - u[0]*bs[2]*k,
		            u[0]*bs[1]*k - u[1]*bs[0]*k);

		case GPARAM_BE_STAR :
		crossp_norm(as, cs, u);
		return vec3(u[1]*cs[2]*l - u[2]*cs[1]*l,
		            u[2]*cs[0]*l - u[0]*cs[2]*l,
		            u[0]*cs[1]*l - u[1]*cs[0]*l);

		case GPARAM_GA_STAR :
		crossp_norm(bs, as, u);
		return vec3(u[1]*as[2]*h - u[2]*as[1]*h,
		            u[2]*as[0]*h - u[0]*as[2]*h,
		            u[0]*as[1]*h - u[1]*as[0]*h);

		case GPARAM_CELL_RX :
		return vec3(0.0, -zl, yl);

		case GPARAM_CELL_RY :
		return vec3(zl, 0.0, -xl);

		case GPARAM_CELL_RZ :
		return vec3(-yl, xl, 0.0);

		default :
		ERROR("Invalid physics gradient parameter %i\n", param);
		abort();

	}
}


static void add_and_free(gsl_vector *total, gsl_vector *n)
{
	gsl_vector_add(total, n);
	gsl_vector_free(n);
}


static gsl_vector *ray_vector_gradient_bravais(int param, UnitCell *cell,
                                               Reflection *refl,
                                               gsl_vector **q)
{
	if ( cell_get_lattice_type(cell) == L_RHOMBOHEDRAL ) {
		if ( param == GPARAM_AL_STAR ) {
			gsl_vector *tot = gsl_vector_calloc(3);
			add_and_free(tot, ray_vector_gradient(param, cell, refl, q));
			add_and_free(tot, ray_vector_gradient(GPARAM_BE_STAR, cell, refl, q));
			add_and_free(tot, ray_vector_gradient(GPARAM_GA_STAR, cell, refl, q));
			return tot;
		}
	}

	switch ( cell_get_lattice_type(cell) ) {

		case L_TRICLINIC :
		case L_MONOCLINIC :
		case L_ORTHORHOMBIC :
		break;

		case L_TETRAGONAL :
		case L_HEXAGONAL :
		switch ( cell_get_unique_axis(cell) ) {

			case 'a' :
			if ( param == GPARAM_B_STAR ) {
				gsl_vector *tot = gsl_vector_calloc(3);
				add_and_free(tot, ray_vector_gradient(param, cell, refl, q));
				add_and_free(tot, ray_vector_gradient(GPARAM_C_STAR, cell, refl, q));
				return tot;
			}
			break;

			case 'b' :
			if ( param == GPARAM_A_STAR ) {
				gsl_vector *tot = gsl_vector_calloc(3);
				add_and_free(tot, ray_vector_gradient(param, cell, refl, q));
				add_and_free(tot, ray_vector_gradient(GPARAM_C_STAR, cell, refl, q));
				return tot;
			}
			break;

			case 'c' :
			if ( param == GPARAM_A_STAR ) {
				gsl_vector *tot = gsl_vector_calloc(3);
				add_and_free(tot, ray_vector_gradient(param, cell, refl, q));
				add_and_free(tot, ray_vector_gradient(GPARAM_B_STAR, cell, refl, q));
				return tot;
			}
			break;

			default :
			ERROR("Invalid unique axis\n");
			return NULL;
		}
		break;

		case L_CUBIC :
		case L_RHOMBOHEDRAL :
		if ( param == GPARAM_A_STAR ) {
			gsl_vector *tot = gsl_vector_calloc(3);
			add_and_free(tot, ray_vector_gradient(param, cell, refl, q));
			add_and_free(tot, ray_vector_gradient(GPARAM_B_STAR, cell, refl, q));
			add_and_free(tot, ray_vector_gradient(GPARAM_C_STAR, cell, refl, q));
			return tot;
		}
		break;

	}

	return ray_vector_gradient(param, cell, refl, q);
}


double r_gradient(int param, Reflection *refl, UnitCell *cell, double wavelength)
{
	gsl_vector *dRdp;
	gsl_vector *q;
	double qdotd, modq;

	if ( param >= GPARAM_CELL_RZ ) {
		/* This is a panel parameter, not a physics parameter, so has
		 * no effect on excitation error */
		return 0;
	}
	dRdp = ray_vector_gradient_bravais(param, cell, refl, &q);

	gsl_vector_set(q, 2, gsl_vector_get(q,2)+1.0/wavelength);
	gsl_blas_ddot(q, dRdp, &qdotd);
	modq = gsl_blas_dnrm2(q);

	gsl_vector_free(q);
	gsl_vector_free(dRdp);

	return -EXC_WEIGHT * qdotd / modq;
}


/* Spot position gradients for diffraction physics (anything that changes the
 * diffracted ray direction) */
int fs_ss_gradient_physics(int param, Reflection *refl, UnitCell *cell,
                           struct detgeom_panel *p, gsl_matrix *Minv,
                           double fs, double ss, double mu,
                           float *fsg, float *ssg)
{
	gsl_vector *dRdp;
	gsl_vector *v;

	dRdp = ray_vector_gradient_bravais(param, cell, refl, NULL);

	v = gsl_vector_calloc(3);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Minv, dRdp, 0.0, v);

	*fsg = mu*(gsl_vector_get(v, 1) - fs*gsl_vector_get(v, 0));
	*ssg = mu*(gsl_vector_get(v, 2) - ss*gsl_vector_get(v, 0));

	gsl_vector_free(v);
	gsl_vector_free(dRdp);

	return 0;
}


/* Spot position gradients for panel motions (translation or rotation) */
int fs_ss_gradient_panel(int param, Reflection *refl, UnitCell *cell,
                         struct detgeom_panel *p, gsl_matrix *Minv,
                         double fs, double ss, double mu,
                         gsl_vector *t, double cx, double cy, double cz,
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
		gsl_matrix_set(dMdp, 1, 0, cz-p->pixel_pitch*p->cnz);
		gsl_matrix_set(dMdp, 2, 0, p->pixel_pitch*p->cny-cy);
		gsl_matrix_set(dMdp, 1, 1, -p->pixel_pitch*p->fsz);
		gsl_matrix_set(dMdp, 2, 1, p->pixel_pitch*p->fsy);
		gsl_matrix_set(dMdp, 1, 2, -p->pixel_pitch*p->ssz);
		gsl_matrix_set(dMdp, 2, 2, p->pixel_pitch*p->ssy);
		break;

		case GPARAM_DET_RY :
		gsl_matrix_set(dMdp, 0, 0, p->pixel_pitch*p->cnz-cz);
		gsl_matrix_set(dMdp, 2, 0, cx-p->pixel_pitch*p->cnx);
		gsl_matrix_set(dMdp, 0, 1, p->pixel_pitch*p->fsz);
		gsl_matrix_set(dMdp, 2, 1, -p->pixel_pitch*p->fsx);
		gsl_matrix_set(dMdp, 0, 2, p->pixel_pitch*p->ssz);
		gsl_matrix_set(dMdp, 2, 2, -p->pixel_pitch*p->ssx);
		break;

		case GPARAM_DET_RZ :
		gsl_matrix_set(dMdp, 0, 0, cy-p->pixel_pitch*p->cny);
		gsl_matrix_set(dMdp, 1, 0, p->pixel_pitch*p->cnx-cx);
		gsl_matrix_set(dMdp, 0, 1, -p->pixel_pitch*p->fsy);
		gsl_matrix_set(dMdp, 1, 1, p->pixel_pitch*p->fsx);
		gsl_matrix_set(dMdp, 0, 2, -p->pixel_pitch*p->ssy);
		gsl_matrix_set(dMdp, 1, 2, p->pixel_pitch*p->ssx);
		break;

		default:
		ERROR("Invalid panel gradient parameter %i\n", param);
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


/* Returns the gradient of fs_dev and ss_dev w.r.t. any parameter.
 * cx,cy,cz are the rotation axis coordinates (only 2 in use at any time)
 * in metres (not pixels) */
int fs_ss_gradient(int param, Reflection *refl, UnitCell *cell,
                   struct detgeom_panel *p, gsl_matrix *Minv,
                   double cx, double cy, double cz,
                   float *fsg, float *ssg)
{
	signed int h, k, l;
	double xl, yl, zl, kpred;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	gsl_vector *t;
	gsl_vector *v;
	gsl_matrix *M;
	double mu;
	double fs, ss;

	get_symmetric_indices(refl, &h, &k, &l);
	kpred = get_kpred(refl);
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	xl = h*asx + k*bsx + l*csx;
	yl = h*asy + k*bsy + l*csy;
	zl = h*asz + k*bsz + l*csz;

	/* Set up matrix equation */
	M = gsl_matrix_alloc(3, 3);
	v = gsl_vector_alloc(3);
	t = gsl_vector_alloc(3);
	if ( (M==NULL) || (v==NULL) || (t==NULL) ) {
		ERROR("Failed to allocate vectors for gradient calculation\n");
		return 1;
	}

	gsl_vector_set(t, 0, xl);
	gsl_vector_set(t, 1, yl);
	gsl_vector_set(t, 2, kpred+zl);

	gsl_matrix_set(M, 0, 0, (p->cnx)*p->pixel_pitch);
	gsl_matrix_set(M, 0, 1, (p->fsx)*p->pixel_pitch);
	gsl_matrix_set(M, 0, 2, (p->ssx)*p->pixel_pitch);
	gsl_matrix_set(M, 1, 0, (p->cny)*p->pixel_pitch);
	gsl_matrix_set(M, 1, 1, (p->fsy)*p->pixel_pitch);
	gsl_matrix_set(M, 1, 2, (p->ssy)*p->pixel_pitch);
	gsl_matrix_set(M, 2, 0, (p->cnz)*p->pixel_pitch);
	gsl_matrix_set(M, 2, 1, (p->fsz)*p->pixel_pitch);
	gsl_matrix_set(M, 2, 2, (p->ssz)*p->pixel_pitch);

	if ( gsl_linalg_HH_solve(M, t, v) ) {
		ERROR("Failed to solve gradient equation\n");
		return 1;
	}
	gsl_matrix_free(M);

	mu = 1.0 / gsl_vector_get(v, 0);
	fs = mu*gsl_vector_get(v, 1);
	ss = mu*gsl_vector_get(v, 2);
	gsl_vector_free(v);

	if ( param <= GPARAM_CELL_RZ ) {
		gsl_vector_free(t);
		return fs_ss_gradient_physics(param, refl, cell, p,
		                              Minv, fs, ss, mu,
		                              fsg, ssg);
	} else {
		return fs_ss_gradient_panel(param, refl, cell, p,
		                            Minv, fs, ss, mu, t,
		                            cx, cy, cz,
		                            fsg, ssg);
	}
}


static int cmpd2(const void *av, const void *bv)
{
	struct reflpeak *a, *b;

	a = (struct reflpeak *)av;
	b = (struct reflpeak *)bv;

	if ( fabs(get_exerr(a->refl)) < fabs(get_exerr(b->refl)) ) return -1;
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
		double grad = fabs(get_exerr(rps[i].refl)) / i;

		for ( j=i+1; j<n; j++ ) {
			if ( fabs(get_exerr(rps[j].refl)) < 0.001e9+grad*j ) {
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
	update_predictions(all_reflist, cr, image);

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
	rps = cfmalloc(image_feature_count(image->features)
	                    * sizeof(struct reflpeak));
	if ( rps == NULL ) return 1;

	reflist = reflist_new();
	n_acc = pair_peaks(image, cr, reflist, rps);
	if ( n_acc < 3 ) {
		cffree(rps);
		reflist_free(reflist);
		return 1;
	}
	update_predictions(reflist, cr, image);

	qsort(rps, n_acc, sizeof(struct reflpeak), cmpd2);
	n = (n_acc-1) - n_acc/50;
	if ( n < 2 ) n = 2; /* n_acc is always >= 2 */
	crystal_set_profile_radius(cr, fabs(get_exerr(rps[n].refl)));

	reflist_free(reflist);
	cffree(rps);

	return 0;
}


void adjust_vector_length(double vec[3], double adj)
{
	double m = modulus(vec[0], vec[1], vec[2]);
	vec[0] += adj*(vec[0])/m;
	vec[1] += adj*(vec[1])/m;
	vec[2] += adj*(vec[2])/m;
}


static void adjust_astar(double bs[3], double cs[3], double shift)
{
	double u[3];
	crossp_norm(cs, bs, u);
	rotate3d(bs, u, shift);
}


static void adjust_bstar(double as[3], double cs[3], double shift)
{
	double u[3];
	crossp_norm(as, cs, u);
	rotate3d(cs, u, shift);
}


static void adjust_cstar(double as[3], double bs[3], double shift)
{
	double u[3];
	crossp_norm(bs, as, u);
	rotate3d(as, u, shift);
}


static int iterate(struct reflpeak *rps, int n,
                   enum gparam *rv, int num_params, UnitCell *cell,
                   struct image *image, gsl_matrix **Minvs)
{
	int i;
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	double as[3], bs[3], cs[3];

	M = gsl_matrix_calloc(num_params, num_params);
	v = gsl_vector_calloc(num_params);

	for ( i=0; i<n; i++ ) {

		int k;
		float fs_gradients[num_params];
		float ss_gradients[num_params];
		float r_gradients[num_params];

		/* Calculate all gradients for this parameter */
		for ( k=0; k<num_params; k++ ) {
			int pn = rps[i].peak->pn;
			r_gradients[k] = r_gradient(rv[k], rps[i].refl, cell, image->lambda);
			fs_ss_gradient(rv[k], rps[i].refl, cell, &image->detgeom->panels[pn],
			               Minvs[pn], 0, 0, 0, &fs_gradients[k], &ss_gradients[k]);
		}

		/* Excitation error terms */
		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = r_gradients[g] * r_gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = r_dev(&rps[i]);
			v_c *= -r_gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		/* Positional fs terms */
		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = fs_gradients[g] * fs_gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = fs_dev(&rps[i], image->detgeom);
			v_c *= -fs_gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		/* Positional ss terms */
		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = ss_gradients[g] * ss_gradients[k];
				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = ss_dev(&rps[i], image->detgeom);
			v_c *= -ss_gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

	}

	int k;
	for ( k=0; k<num_params; k++ ) {
		double M_curr = gsl_matrix_get(M, k, k);
		gsl_matrix_set(M, k, k, M_curr+1e-18);
	}

	//show_matrix_eqn(M, v);
	shifts = solve_svd(v, M, NULL, 0);
	if ( shifts == NULL ) {
		ERROR("Failed to solve equations.\n");
		gsl_matrix_free(M);
		gsl_vector_free(v);
		return 1;
	}

	/* Apply shifts */
	cell_get_reciprocal(cell, &as[0], &as[1], &as[2],
	                          &bs[0], &bs[1], &bs[2],
	                          &cs[0], &cs[1], &cs[2]);
	for ( i=0; i<num_params; i++ ) {

		double shift = gsl_vector_get(shifts, i);
		if ( isnan(shift) ) shift = 0.0;

		LatticeType lt = cell_get_lattice_type(cell);
		char ua = cell_get_unique_axis(cell);

		switch ( rv[i] ) {

			case GPARAM_A_STAR :
			adjust_vector_length(as, shift);
			if ( (lt == L_TETRAGONAL) || (lt == L_HEXAGONAL) ) {
				if ( ua == 'b' ) {
					adjust_vector_length(cs, shift);
				} else if ( ua == 'c' ) {
					adjust_vector_length(bs, shift);
				}
			}
			if ( (lt == L_CUBIC) || (lt == L_RHOMBOHEDRAL) ) {
				adjust_vector_length(bs, shift);
				adjust_vector_length(cs, shift);
			}
			break;

			case GPARAM_B_STAR :
			adjust_vector_length(bs, shift);
			if ( (lt == L_TETRAGONAL) || (lt == L_HEXAGONAL) ) {
				if ( ua == 'a' ) {
					adjust_vector_length(cs, shift);
				}
			}
			break;

			case GPARAM_C_STAR :
			adjust_vector_length(cs, shift);
			break;

			case GPARAM_AL_STAR :
			adjust_astar(bs, cs, shift);
			if ( lt == L_RHOMBOHEDRAL ) {
				adjust_bstar(as, cs, shift);
				adjust_cstar(as, bs, shift);
			}
			break;

			case GPARAM_BE_STAR :
			adjust_bstar(as, cs, shift);
			break;

			case GPARAM_GA_STAR :
			adjust_cstar(as, bs, shift);
			break;

			case GPARAM_CELL_RX :
			rotate2d(&as[1], &as[2], 0.0, 0.0, shift);
			rotate2d(&bs[1], &bs[2], 0.0, 0.0, shift);
			rotate2d(&cs[1], &cs[2], 0.0, 0.0, shift);
			break;

			case GPARAM_CELL_RY :
			rotate2d(&as[2], &as[0], 0.0, 0.0, shift);
			rotate2d(&bs[2], &bs[0], 0.0, 0.0, shift);
			rotate2d(&cs[2], &cs[0], 0.0, 0.0, shift);
			break;

			case GPARAM_CELL_RZ :
			rotate2d(&as[0], &as[1], 0.0, 0.0, shift);
			rotate2d(&bs[0], &bs[1], 0.0, 0.0, shift);
			rotate2d(&cs[0], &cs[1], 0.0, 0.0, shift);
			break;

			default :
			ERROR("Unrecognised parameter %i\n", rv[i]);
			break;

		}

	}
	cell_set_reciprocal(cell, as[0], as[1], as[2],
	                          bs[0], bs[1], bs[2],
	                          cs[0], cs[1], cs[2]);
	gsl_vector_free(shifts);
	gsl_matrix_free(M);
	gsl_vector_free(v);

	return 0;
}


static double pred_residual(struct reflpeak *rps, int n, struct detgeom *det,
                            double *pres_r, double *pres_fs, double *pres_ss)
{
	int i;
	double res_r, res_fs, res_ss;

	res_r = 0.0;
	res_fs = 0.0;
	res_ss = 0.0;
	for ( i=0; i<n; i++ ) {
		res_r += sq(r_dev(&rps[i]));
		res_fs += sq(fs_dev(&rps[i], det));
		res_ss += sq(ss_dev(&rps[i], det));
	}

	if ( pres_r != NULL ) *pres_r = res_r;
	if ( pres_fs != NULL ) *pres_fs = res_fs;
	if ( pres_ss != NULL ) *pres_ss = res_ss;

	return res_r + res_fs + res_ss;
}


/* NB Only for use when the list of reflpeaks was created without a RefList.
 * If a RefList was used, then reflist_free the list then just cffree() the rps */
static void free_rps_noreflist(struct reflpeak *rps, int n)
{
	int i;

	for ( i=0; i<n; i++ ) {
		reflection_free(rps[i].refl);
	}
	cffree(rps);
}


static int parameters_to_refine(UnitCell *cell, enum gparam *rv)
{
	int num_params = 0;

	switch ( cell_get_lattice_type(cell) ) {

		case L_TRICLINIC :
		case L_MONOCLINIC :
		case L_ORTHORHOMBIC :
		rv[num_params++] = GPARAM_A_STAR;
		rv[num_params++] = GPARAM_B_STAR;
		rv[num_params++] = GPARAM_C_STAR;
		break;

		case L_TETRAGONAL :
		case L_HEXAGONAL :
		if ( cell_get_unique_axis(cell) == 'a' ) {
			rv[num_params++] = GPARAM_A_STAR;
			rv[num_params++] = GPARAM_B_STAR; /* == GPARAM_C_STAR */
		} else if ( cell_get_unique_axis(cell) == 'b' ) {
			rv[num_params++] = GPARAM_A_STAR; /* == GPARAM_C_STAR */
			rv[num_params++] = GPARAM_B_STAR;
		} else if ( cell_get_unique_axis(cell) == 'c' ) {
			rv[num_params++] = GPARAM_A_STAR; /* == GPARAM_B_STAR */
			rv[num_params++] = GPARAM_C_STAR;
		} else {
			ERROR("Unrecognised unique axis:\n");
			cell_print(cell);
			return 0;
		}
		break;

		case L_CUBIC :
		case L_RHOMBOHEDRAL :
		rv[num_params++] = GPARAM_A_STAR;
		break;

		default :
		ERROR("Unrecognised lattice type:\n");
		cell_print(cell);
		return 0;

	};

	switch ( cell_get_lattice_type(cell) ) {

		case L_TRICLINIC :
		rv[num_params++] = GPARAM_AL_STAR;
		rv[num_params++] = GPARAM_BE_STAR;
		rv[num_params++] = GPARAM_GA_STAR;
		break;

		case L_MONOCLINIC :
		if ( cell_get_unique_axis(cell) == 'a' ) {
			rv[num_params++] = GPARAM_AL_STAR;
		} else if ( cell_get_unique_axis(cell) == 'b' ) {
			rv[num_params++] = GPARAM_BE_STAR;
		} else if ( cell_get_unique_axis(cell) == 'c' ) {
			rv[num_params++] = GPARAM_GA_STAR;
		} else {
			ERROR("Unrecognised unique axis:\n");
			cell_print(cell);
			return 0;
		}

		case L_RHOMBOHEDRAL :
		rv[num_params++] = GPARAM_AL_STAR; /* == beta and gamma */
		break;

		case L_ORTHORHOMBIC :
		case L_TETRAGONAL :
		case L_HEXAGONAL :
		case L_CUBIC :
		break;

		default :
		ERROR("Unrecognised lattice type:\n");
		cell_print(cell);
		return 0;
	};

	/* Always refine orientation */
	rv[num_params++] = GPARAM_CELL_RX;
	rv[num_params++] = GPARAM_CELL_RY;
	rv[num_params++] = GPARAM_CELL_RZ;
	return num_params;
}


int refine_prediction(struct image *image, Crystal *cr,
                      Mille *mille, int max_mille_level,
                      UnitCell *target)
{
	int n;
	int i;
	struct reflpeak *rps;
	double max_I;
	RefList *reflist;
	char tmp[256];
	gsl_matrix **Minvs;
	double res_r, res_fs, res_ss, res_overall;

	enum gparam rv[] = {
		GPARAM_A_STAR,
		GPARAM_B_STAR,
		GPARAM_C_STAR,
		GPARAM_AL_STAR,
		GPARAM_BE_STAR,
		GPARAM_GA_STAR,
		GPARAM_CELL_RX,
		GPARAM_CELL_RY,
		GPARAM_CELL_RZ,
	};
	int num_params = 9;

	rps = cfmalloc(image_feature_count(image->features)
	                         * sizeof(struct reflpeak));
	if ( rps == NULL ) return 1;

	reflist = reflist_new();
	n = pair_peaks(image, cr, reflist, rps);
	if ( n < 3 ) {
		cffree(rps);
		reflist_free(reflist);
		return 1;
	}

	Minvs = make_panel_minvs(image->detgeom);

	/* Normalise the intensities to max 1 */
	max_I = -INFINITY;
	for ( i=0; i<n; i++ ) {
		double cur_I = rps[i].peak->intensity;
		if ( cur_I > max_I ) max_I = cur_I;
	}
	if ( max_I <= 0.0 ) {
		ERROR("All peaks negative?\n");
		cffree(rps);
		return 1;
	}
	for ( i=0; i<n; i++ ) {
		if ( rps[i].peak->intensity > 0.0 ) {
			rps[i].Ih = rps[i].peak->intensity / max_I;
		} else {
			rps[i].Ih = 0.0;
		}
	}

	res_overall = pred_residual(rps, n, image->detgeom, &res_r, &res_fs, &res_ss);
	snprintf(tmp, 255, "predict_refine/initial_residual = %f (%f %f %f)",
	         res_overall, res_r, res_fs, res_ss);
	crystal_add_notes(cr, tmp);
	//STATUS("Initial residual = %f (%f %f %f)\n",
	//       res_overall, res_r, res_fs, res_ss);

	/* Pretend it's triclinic for now */
	cell_set_lattice_type(crystal_get_cell(cr), L_TRICLINIC);

	/* Refine (max 5 cycles) */
	for ( i=0; i<5; i++ ) {
		update_predictions(reflist, cr, image);
		if ( iterate(rps, n, rv, num_params, crystal_get_cell(cr), image, Minvs) )
		{
			return 1;
		}

		res_overall = pred_residual(rps, n, image->detgeom, &res_r, &res_fs, &res_ss);
		//STATUS("Residual after %i = %f (%f %f %f)\n",
		//       i, res_overall, res_r, res_fs, res_ss);
	}

	if ( target == NULL ) goto done;

	UnitCell *nc = impose_bravais(crystal_get_cell(cr),
	                              cell_get_lattice_type(target),
	                              cell_get_unique_axis(target));
	if ( nc == NULL ) {
		ERROR("Failed to impose Bravais conditions\n");
		return 1;
	}
	crystal_set_cell(cr, nc);
	num_params = parameters_to_refine(nc, rv);
	if ( num_params == 0 ) {
		ERROR("Couldn't determine which parameters to refine\n");
		return 1;
	}

	update_predictions(reflist, cr, image);
	res_overall = pred_residual(rps, n, image->detgeom, &res_r, &res_fs, &res_ss);
	//STATUS("After applying Bravais constraints = %f (%f %f %f)\n",
	//       res_overall, res_r, res_fs, res_ss);

	/* Refine again, with Bravais constraints (max 5 cycles) */
	for ( i=0; i<5; i++ ) {
		update_predictions(reflist, cr, image);
		if ( iterate(rps, n, rv, num_params, crystal_get_cell(cr), image, Minvs) )
		{
			return 1;
		}

		res_overall = pred_residual(rps, n, image->detgeom, &res_r, &res_fs, &res_ss);
		//STATUS("Residual after %i = %f (%f %f %f)\n",
		//       i, res_overall, res_r, res_fs, res_ss);
	}

done:
	res_overall = pred_residual(rps, n, image->detgeom, &res_r, &res_fs, &res_ss);
	snprintf(tmp, 255, "predict_refine/final_residual = %f (%f %f %f)",
	         res_overall, res_r, res_fs, res_ss);
	crystal_add_notes(cr, tmp);
	//STATUS("Final residual = %f (%f %f %f)\n",
	//       res_overall, res_r, res_fs, res_ss);

	if ( (mille != NULL) && (n>4) ) {
		crystfel_mille_delete_last_record(mille);
		profile_start("mille-calc");
		write_mille(mille, n, crystal_get_cell(cr), rv, num_params,
		            rps, image, max_mille_level, Minvs);
		profile_end("mille-calc");
	}

	for ( i=0; i<image->detgeom->n_panels; i++ ) {
		gsl_matrix_free(Minvs[i]);
	}
	cffree(Minvs);

	reflist_free(reflist);

	n = pair_peaks(image, cr, NULL, rps);
	free_rps_noreflist(rps, n);

	if ( n < 3 ) {
		return 1;
	} else {
		return 0;
	}
}
