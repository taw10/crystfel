/*
 * integration.c
 *
 * Integration of intensities
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <ncurses.h>

#include "reflist.h"
#include "cell.h"
#include "crystal.h"
#include "cell-utils.h"
#include "geometry.h"
#include "image.h"
#include "peaks.h"
#include "integration.h"


#define VERBOSITY (0)
// ((h==-6) && (k==0) && (l==-8))


static void check_eigen(gsl_vector *e_val)
{
	int i;
	double vmax, vmin;
	const int n = e_val->size;
	const double max_condition = 1e6;
	const int verbose = 0;

	if ( verbose ) STATUS("Eigenvalues:\n");
	vmin = +INFINITY;
	vmax = 0.0;
	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( verbose ) STATUS("%i: %e\n", i, val);
		if ( val > vmax ) vmax = val;
		if ( val < vmin ) vmin = val;
	}

	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( val < vmax/max_condition ) {
			gsl_vector_set(e_val, i, 0.0);
		}
	}

	vmin = +INFINITY;
	vmax = 0.0;
	for ( i=0; i<n; i++ ) {
		double val = gsl_vector_get(e_val, i);
		if ( val == 0.0 ) continue;
		if ( val > vmax ) vmax = val;
		if ( val < vmin ) vmin = val;
	}
	if ( verbose ) {
		STATUS("Condition number: %e / %e = %5.2f\n",
		       vmax, vmin, vmax/vmin);
	}
}


static gsl_vector *solve_svd(gsl_vector *v, gsl_matrix *M)
{
	gsl_matrix *s_vec;
	gsl_vector *s_val;
	int err, n;
	gsl_vector *shifts;

	n = v->size;
	if ( v->size != M->size1 ) return NULL;
	if ( v->size != M->size2 ) return NULL;

	s_val = gsl_vector_calloc(n);
	s_vec = gsl_matrix_calloc(n, n);

	err = gsl_linalg_SV_decomp_jacobi(M, s_vec, s_val);
	if ( err ) {
		ERROR("SVD failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		return NULL;
	}
	/* "M" is now "U" */

	check_eigen(s_val);

	shifts = gsl_vector_calloc(n);
	err = gsl_linalg_SV_solve(M, s_vec, s_val, v, shifts);
	if ( err ) {
		ERROR("Matrix solution failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		gsl_vector_free(shifts);
		return NULL;
	}

	gsl_matrix_free(s_vec);
	gsl_vector_free(s_val);

	return shifts;
}


enum boxmask_val
{
	BM_IG,  /* "Soft" ignore */
	BM_BH,  /* "Hard" ignore (black hole) */
	BM_BG,  /* Background */
	BM_PK   /* Peak */
};


struct intcontext
{
	int halfw;
	int w;
	enum boxmask_val *bm;  /* Box mask */
	struct image *image;
	gsl_matrix *bgm;  /* Background estimation matrix */

	struct peak_box *boxes;
	int n_boxes;
	int max_boxes;

	UnitCell *cell;
	double k;

	/* Peak region sums */
	double pks_p2;
	double pks_q2;
	double pks_pq;
	double pks_p;
	double pks_q;
	int m;

	int n_reference_profiles;
	double **reference_profiles;
	double **reference_den;
	int *n_profiles_in_reference;

	/* For summation integration only */
	int ir_inn;
	int ir_mid;
	int ir_out;
};


struct peak_box
{
	int cfs;   /* Coordinates of corner */
	int css;

	enum boxmask_val *bm;  /* Box mask */

	int pn;           /* Panel number */
	struct panel *p;  /* The panel itself */

	/* Fitted background parameters */
	double a;
	double b;
	double c;

	/* Measured intensity (tentative, profile fitted or otherwise) */
	double intensity;
	double sigma;
	double J;  /* Profile scaling factor */

	/* Offsets to final observed position */
	double offs_fs;
	double offs_ss;

	int rp;   /* Reference profile number */

	Reflection *refl;

	int verbose;
};


static void addm(gsl_matrix *m, int i, int j, double val)
{
	double v = gsl_matrix_get(m, i, j);
	gsl_matrix_set(m, i, j, v+val);
}


static void addv(gsl_vector *v, int i, double val)
{
	double k = gsl_vector_get(v, i);
	gsl_vector_set(v, i, k+val);
}


static float boxi(struct intcontext *ic, struct peak_box *bx, int p, int q)
{
	int fs, ss;

	fs = bx->cfs + p;
	ss = bx->css + q;

	assert(fs >= 0);
	assert(fs < bx->p->w);
	assert(ss >= 0);
	assert(ss < bx->p->h);

	assert(p >= 0);
	assert(p < ic->w);
	assert(q >= 0);
	assert(q < ic->w);

	return ic->image->dp[bx->pn][fs + bx->p->w*ss];
}


static void colour_on(enum boxmask_val b)
{
	switch ( b ) {

		case BM_BG :
		attron(COLOR_PAIR(1));
		break;

		case BM_PK :
		attron(COLOR_PAIR(2));
		break;

		case BM_BH :
		attron(COLOR_PAIR(3));
		break;

		default:
		break;

	}
}


static void colour_off(enum boxmask_val b)
{
	switch ( b ) {

		case BM_BG :
		attroff(COLOR_PAIR(1));
		break;

		case BM_PK :
		attroff(COLOR_PAIR(2));
		break;

		case BM_BH :
		attroff(COLOR_PAIR(3));
		break;

		default:
		break;

	}
}


static void show_peak_box(struct intcontext *ic, struct peak_box *bx)
{
	int q;

	initscr();
	clear();
	start_color();
	init_pair(1, COLOR_WHITE, COLOR_BLUE) ;  /* Background */
	init_pair(2, COLOR_WHITE, COLOR_RED);    /* Peak */
	init_pair(3, COLOR_BLACK, COLOR_CYAN);   /* Blackhole */

	printw("Pixel values:\n");
	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(bx->bm[p+q*ic->w]);
			printw("%5.0f ", boxi(ic, bx, p, q));
			colour_off(bx->bm[p+q*ic->w]);

		}

		printw("\n");
	}

	printw("\nFitted background:\n");
	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(bx->bm[p+q*ic->w]);
			printw("%5.0f ", bx->a*p + bx->b*q + bx->c);
			colour_off(bx->bm[p+q*ic->w]);

		}

		printw("\n");
	}

	printw("Reference profile number %i, ", bx->rp);
	printw("Background parameters: a=%.2f, b=%.2f, c=%.2f\n",
	       bx->a, bx->b, bx->c);
	getch();
	refresh();
	endwin();
}


static void show_reference_profile(struct intcontext *ic, int i)
{
	int q;

	initscr();
	clear();
	start_color();
	init_pair(1, COLOR_WHITE, COLOR_BLUE) ;  /* Background */
	init_pair(2, COLOR_WHITE, COLOR_RED);    /* Peak */

	printw("Reference profile number %i (%i contributions):\n", i,
	       ic->n_profiles_in_reference[i]);

	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(ic->bm[p+q*ic->w]);
			printw("%3.0f ", ic->reference_profiles[i][p+ic->w*q]);
			colour_off(ic->bm[p+q*ic->w]);

		}

		printw("\n");
	}

	getch();
	refresh();
	endwin();
}


static void fit_bg(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	gsl_vector *v;
	gsl_vector *ans;
	gsl_matrix *M;

	v = gsl_vector_calloc(3);

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		if ( bx->bm[p + ic->w*q] == BM_BG ) {

			bi = boxi(ic, bx, p, q);

			addv(v, 0, bi*p);
			addv(v, 1, bi*q);
			addv(v, 2, bi);

		}

	}
	}

	if ( bx->verbose ) {
		show_matrix_eqn(ic->bgm, v, 3);
	}

	M = gsl_matrix_alloc(3, 3);
	gsl_matrix_memcpy(M, ic->bgm);
	ans = solve_svd(v, M);
	gsl_vector_free(v);
	gsl_matrix_free(M);

	bx->a = gsl_vector_get(ans, 0);
	bx->b = gsl_vector_get(ans, 1);
	bx->c = gsl_vector_get(ans, 2);

	gsl_vector_free(ans);

	if ( bx->verbose ) {
		show_peak_box(ic, bx);
	}
}


static void zero_profiles(struct intcontext *ic)
{
	int i;

	for ( i=0; i<ic->n_reference_profiles; i++ ) {

		int p, q;

		for ( p=0; p<ic->w; p++ ) {
		for ( q=0; q<ic->w; q++ ) {
			ic->reference_profiles[i][p+ic->w*q] = 0.0;
			ic->reference_den[i][p+ic->w*q] = 0.0;
		}
		}

		ic->n_profiles_in_reference[i] = 0;
	}
}



static int alloc_boxes(struct intcontext *ic, int new_max_boxes)
{
	struct peak_box *boxes_new;

	boxes_new = realloc(ic->boxes, sizeof(struct peak_box)*new_max_boxes);
	if ( boxes_new == NULL ) return 1;

	ic->boxes = boxes_new;
	ic->max_boxes = new_max_boxes;
	return 0;
}


static int init_intcontext(struct intcontext *ic)
{
	int p, q;
	int i;

	ic->w = 2*ic->halfw + 1;

	ic->bm = malloc(ic->w * ic->w * sizeof(enum boxmask_val));
	if ( ic->bm == NULL ) {
		ERROR("Failed to allocate box mask.\n");
		return 1;
	}

	ic->bgm = gsl_matrix_calloc(3, 3);
	if ( ic->bgm == NULL ) {
		ERROR("Failed to initialise matrix.\n");
		return 1;
	}

	ic->pks_p2 = 0.0;
	ic->pks_q2 = 0.0;
	ic->pks_pq = 0.0;
	ic->pks_p = 0.0;
	ic->pks_q = 0.0;
	ic->m = 0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		if ( (p==0) || (q==0) || (p==ic->w-1) || (q==ic->w-1) ) {
			ic->bm[p + ic->w*q] = BM_BG;
		} else {
			ic->bm[p + ic->w*q] = BM_PK;
		}

		switch ( ic->bm[p + ic->w*q] ) {

			case BM_IG :
			case BM_BH :
			break;

			case BM_BG :
			addm(ic->bgm, 0, 0, p*p);
			addm(ic->bgm, 0, 1, p*q);
			addm(ic->bgm, 0, 2, p);
			addm(ic->bgm, 1, 0, p*q);
			addm(ic->bgm, 1, 1, q*q);
			addm(ic->bgm, 1, 2, q);
			addm(ic->bgm, 2, 0, p);
			addm(ic->bgm, 2, 1, q);
			addm(ic->bgm, 2, 2, 1);
			break;

			case BM_PK :
			ic->pks_p2 += p*p;
			ic->pks_q2 += q*q;
			ic->pks_pq += p*q;
			ic->pks_p += p;
			ic->pks_q += q;
			ic->m++;
			break;

		}

	}
	}

	/* How many reference profiles? */
	ic->n_reference_profiles = ic->image->det->n_panels;
	ic->reference_profiles = calloc(ic->n_reference_profiles,
	                                sizeof(double *));
	if ( ic->reference_profiles == NULL ) return 1;
	ic->reference_den = calloc(ic->n_reference_profiles, sizeof(double *));
	if ( ic->reference_den == NULL ) return 1;
	ic->n_profiles_in_reference = calloc(ic->n_reference_profiles,
	                                     sizeof(int));
	if ( ic->n_profiles_in_reference == NULL ) return 1;
	for ( i=0; i<ic->n_reference_profiles; i++ ) {
		ic->reference_profiles[i] = malloc(ic->w*ic->w*sizeof(double));
		if ( ic->reference_profiles[i] == NULL ) return 1;
		ic->reference_den[i] = malloc(ic->w*ic->w*sizeof(double));
		if ( ic->reference_den[i] == NULL ) return 1;
	}
	zero_profiles(ic);

	ic->boxes = NULL;
	ic->n_boxes = 0;
	ic->max_boxes = 0;
	if ( alloc_boxes(ic, 32) ) {
		return 1;
	}

	return 0;
}


static void free_intcontext(struct intcontext *ic)
{
	int i;

	free(ic->boxes);
	for ( i=0; i<ic->n_reference_profiles; i++ ) {
		free(ic->reference_profiles[i]);
		free(ic->reference_den[i]);
	}
	free(ic->reference_profiles);
	free(ic->reference_den);
	free(ic->n_profiles_in_reference);
	free(ic->bm);
	gsl_matrix_free(ic->bgm);
}


static void setup_ring_masks(struct intcontext *ic,
                             double ir_inn, double ir_mid, double ir_out)
{
	double lim_sq, out_lim_sq, mid_lim_sq;
	int p, q;

	lim_sq = pow(ir_inn, 2.0);
	mid_lim_sq = pow(ir_mid, 2.0);
	out_lim_sq = pow(ir_out, 2.0);

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		int rsq;

		rsq = (p-ic->halfw)*(p-ic->halfw) + (q-ic->halfw)*(q-ic->halfw);

		if ( rsq > out_lim_sq ) {
			/* Outside outer radius */
			ic->bm[p + ic->w*q] = BM_IG;
		} else {

			if ( rsq >= mid_lim_sq ) {
				/* Inside outer radius, outside middle radius */
				ic->bm[p + ic->w*q] = BM_BG;
			} else if ( rsq <= lim_sq ) {
				/* Inside inner radius */
				ic->bm[p + ic->w*q] = BM_PK;
			} else {
				/* Outside inner radius, inside middle radius */
				ic->bm[p + ic->w*q] = BM_IG;
			}

		}

	}
	}

}


static struct peak_box *add_box(struct intcontext *ic)
{
	int idx;

	if ( ic->n_boxes == ic->max_boxes ) {
		if ( alloc_boxes(ic, ic->max_boxes+32) ) {
			return NULL;
		}
	}

	idx = ic->n_boxes++;

	ic->boxes[idx].cfs = 0;
	ic->boxes[idx].css = 0;
	ic->boxes[idx].bm = NULL;
	ic->boxes[idx].pn = -1;
	ic->boxes[idx].p = NULL;
	ic->boxes[idx].a = 0.0;
	ic->boxes[idx].b = 0.0;
	ic->boxes[idx].c = 0.0;
	ic->boxes[idx].intensity = 0.0;
	ic->boxes[idx].sigma = 0.0;
	ic->boxes[idx].J = 0.0;
	ic->boxes[idx].rp = -1;
	ic->boxes[idx].refl = NULL;
	ic->boxes[idx].verbose = 0;

	return &ic->boxes[idx];
}


static void delete_box(struct intcontext *ic, struct peak_box *bx)
{
	int i;
	int found = 0;

	for ( i=0; i<ic->n_boxes; i++ ) {
		if ( &ic->boxes[i] == bx ) {
			found = 1;
			break;
		}
	}

	if ( !found ) {
		ERROR("Couldn't find box %p in context %p\n", bx, ic);
		return;
	}

	memmove(&ic->boxes[i], &ic->boxes[i+1],
	        (ic->n_boxes-i-1)*sizeof(struct peak_box));
	ic->n_boxes--;
}


static double tentative_intensity(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double intensity = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;
		intensity += boxi(ic, bx, p, q);

	}
	}

	intensity -= bx->a * ic->pks_p;
	intensity -= bx->b * ic->pks_q;
	intensity -= bx->c * ic->m;

	return intensity;
}


static void observed_position(struct intcontext *ic, struct peak_box *bx,
                              double *pos_p, double *pos_q)
{
	int p, q;
	double num_p = 0.0;
	double num_q = 0.0;
	double den = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		int bi;
		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;
		bi = boxi(ic, bx, p, q);

		num_p += bi*(p - ic->halfw);
		num_q += bi*(q - ic->halfw);
		den += bi;

	}
	}

	num_p += -bx->a*ic->pks_p2 - bx->b*ic->pks_pq - bx->c*ic->pks_p;
	num_q += -bx->a*ic->pks_q2 - bx->b*ic->pks_pq - bx->c*ic->pks_q;
	den += -bx->a*ic->pks_p - bx->b*ic->pks_q - bx->c;

	*pos_p = num_p / den;
	*pos_q = num_q / den;
}


static void add_to_reference_profile(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double val;
		float bi;

		if ( bx->bm[p + ic->w*q] == BM_IG ) continue;
		if ( bx->bm[p + ic->w*q] == BM_BH ) continue;
		bi = boxi(ic, bx, p, q);

		val = bi*bx->intensity;
		val -= p*bx->intensity*bx->a + q*bx->intensity*bx->b
		       + bx->intensity*bx->c;

		ic->reference_profiles[bx->rp][p+ic->w*q] += val;
		ic->reference_den[bx->rp][p+ic->w*q] += pow(bx->intensity, 2.0);
	}
	}

	ic->n_profiles_in_reference[bx->rp]++;
}


static void calculate_reference_profiles(struct intcontext *ic)
{
	int i;

	for ( i=0; i<ic->n_reference_profiles; i++ ) {

		int p, q;
		double max = 0.0;

		for ( p=0; p<ic->w; p++ ) {
		for ( q=0; q<ic->w; q++ ) {

			double den;

			den = ic->reference_den[i][p+ic->w*q];
			ic->reference_profiles[i][p+ic->w*q] /= den;
			if ( ic->reference_profiles[i][p+ic->w*q] > max ) {
				max = ic->reference_profiles[i][p+ic->w*q];
			}

		}
		}

		max /= 100.0;

		for ( p=0; p<ic->w; p++ ) {
		for ( q=0; q<ic->w; q++ ) {

			ic->reference_profiles[i][p+ic->w*q] /= max;

		}
		}

	}

	//for ( i=0; i<ic->n_reference_profiles; i++ ) {
	//	show_reference_profile(ic, i);
	//}
}


static int check_box(struct intcontext *ic, struct peak_box *bx, int *sat)
{
	int p, q;
	int n_pk = 0;
	int n_bg = 0;
	double adx, ady, adz;
	double bdx, bdy, bdz;
	double cdx, cdy, cdz;
	signed int hr, kr, lr;

	if ( sat != NULL ) *sat = 0;

	bx->bm = malloc(ic->w*ic->w*sizeof(int));
	if ( bx->bm == NULL ) return 1;

	cell_get_cartesian(ic->cell,
	                   &adx, &ady, &adz,
	                   &bdx, &bdy, &bdz,
	                   &cdx, &cdy, &cdz);
	get_indices(bx->refl, &hr, &kr, &lr);

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		int fs, ss;
		double hd, kd, ld;
		signed int h, k, l;
		struct rvec dv;

		fs = bx->cfs + p;
		ss = bx->css + q;

		if ( (fs < 0) || (fs >= bx->p->w)
		  || (ss < 0) || (ss >= bx->p->h) ) return 1;

		if ( (p < 0) || (p >= ic->w) || (q < 0) || (q >= ic->w) ) {
			return 1;
		}

		bx->bm[p+ic->w*q] = ic->bm[p+ic->w*q];

		if ( ic->image->bad[bx->pn][fs + bx->p->w*ss] ) {
			bx->bm[p+ic->w*q] = BM_BH;
		}

		if ( (bx->bm[p+ic->w*q] != BM_IG)
		  && (bx->bm[p+ic->w*q] != BM_BH)
		  && (boxi(ic, bx, p, q) > bx->p->max_adu) ) {
			if ( sat != NULL ) *sat = 1;
		}

		/* Ignore if this pixel is closer to the next reciprocal lattice
		 * point */
		dv = get_q_for_panel(bx->p, fs, ss, NULL, ic->k);
		hd = dv.u * adx + dv.v * ady + dv.w * adz;
		kd = dv.u * bdx + dv.v * bdy + dv.w * bdz;
		ld = dv.u * cdx + dv.v * cdy + dv.w * cdz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);
		if ( (h != hr) || (k != kr) || (l != lr) ) {
			bx->bm[p+ic->w*q] = BM_BH;
		}

		if ( bx->bm[p+ic->w*q] == BM_PK ) n_pk++;
		if ( bx->bm[p+ic->w*q] == BM_BG ) n_bg++;

	}
	}

	if ( n_pk < 4 ) return 1;
	if ( n_bg < 4 ) return 1;

	return 0;
}


static double fit_J(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double sum = 0.0;
	double den = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi, P;

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;

		bi = boxi(ic, bx, p, q);
		P = ic->reference_profiles[bx->rp][p+ic->w*q];

		sum += bi*P;
		sum += - bx->a*p*P - bx->b*q*P - bx->c*P;
		den += pow(P, 2.0);

	}
	}

	return sum / den;
}


static int center_and_check_box(struct intcontext *ic, struct peak_box *bx,
                                int *sat)
{
	int i;

	bx->offs_fs = 0.0;
	bx->offs_ss = 0.0;

	if ( check_box(ic, bx, sat) ) return 1;

	for ( i=0; i<10; i++ ) {

		int p, q;
		double sum_fs = 0.0;
		double sum_ss = 0.0;
		double den = 0.0;
		int t_offs_fs = 0;
		int t_offs_ss = 0;
		double offs_fs, offs_ss;
		int ifs, iss;

		for ( p=0; p<ic->w; p++ ) {
		for ( q=0; q<ic->w; q++ ) {

			double bi = boxi(ic, bx, p, q);

			if ( bx->bm[p + ic->w*q] == BM_BH ) continue;

			sum_fs += bi * (p-ic->halfw);
			sum_ss += bi * (q-ic->halfw);
			den += bi;

		}
		}

		offs_fs = sum_fs / den;
		offs_ss = sum_ss / den;

		ifs = rint(offs_fs);
		iss = rint(offs_ss);
		bx->offs_fs += ifs;
		bx->offs_ss += iss;
		bx->cfs += ifs;
		bx->css += iss;

		t_offs_fs += rint(offs_fs);
		t_offs_ss += rint(offs_fs);

		if ( check_box(ic, bx, sat) ) return 1;

		if ( t_offs_fs*t_offs_fs + t_offs_ss*t_offs_ss > ic->w*ic->w ) {
			return 1;
		}
	}

	return 0;
}


static double fit_intensity(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double J = fit_J(ic, bx);
	double sum = 0.0;

	bx->J = J;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double P;

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;

		P = ic->reference_profiles[bx->rp][p+ic->w*q];
		sum += P;

	}
	}

	if ( bx->verbose ) {
		STATUS("J = %f\n", J);
	}

	return J * sum;
}



static double calc_sigma(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double sum = 0.0;
	double mb = 0.0;
	int nb = 0;
	double sigb2 = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		bi = boxi(ic, bx, p, q);

		if ( bx->bm[p + ic->w*q] == BM_PK ) {

			double p1, p2;

			p1 = bx->J * ic->reference_profiles[bx->rp][p+ic->w*q];

			p2 = bi - bx->a*p - bx->b*q - bx->c;
			sum += pow(p1-p2, 2.0);

		} else if ( bx->bm[p + ic->w*q] == BM_BG ) {

			mb += bi;
			nb++;

		}

	}
	}

	mb /= nb;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		if ( bx->bm[p + ic->w*q] != BM_BG ) continue;
		bi = boxi(ic, bx, p, q);
		sigb2 += pow(bi - mb, 2.0);

	}
	}

	return sqrt(sum + sigb2);
}


static double average_background(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double sum = 0.0;
	int n = 0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;

		sum += bx->a*p + bx->b*q + bx->c;
		n++;

	}
	}

	return sum/n;
}


static int bg_ok(struct peak_box *bx)
{
	if ( (fabs(bx->a) > 10.0) || (fabs(bx->b) > 10.0) ) {
		return 0;
	} else {
		return 1;
	}
}


static int suitable_reference(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double max = 0.0;
	int height_ok;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;

		bi = boxi(ic, bx, p, q);
		if ( bi > max ) max = bi;

	}
	}

	height_ok = max > 10.0 * average_background(ic, bx);

	return bg_ok(bx) && height_ok;
}


static void measure_all_intensities(IntegrationMethod meth, RefList *list,
                                    struct image *image, UnitCell *cell)
{
	Reflection *refl;
	RefListIterator *iter;
	struct intcontext ic;
	int i;
	int n_saturated = 0;

	ic.halfw = 4;
	ic.image = image;
	ic.k = 1.0/image->lambda;
	ic.cell = cell;
	if ( init_intcontext(&ic) ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pfs, pss;
		signed int h, k, l;
		struct peak_box *bx;
		int pn;
		struct panel *p;
		int fid_fs, fid_ss;  /* Center coordinates, rounded,
		                      * in overall data block */
		int cfs, css;  /* Corner coordinates */
		int saturated;
		int r;

		set_redundancy(refl, 0);

		get_detector_pos(refl, &pfs, &pss);
		fid_fs = lrint(pfs);
		fid_ss = lrint(pss);
		pn = find_panel_number(image->det, fid_fs, fid_ss);
		p = &image->det->panels[pn];

		cfs = (fid_fs-p->min_fs) - ic.halfw;
		css = (fid_ss-p->min_ss) - ic.halfw;

		if ( (cfs + ic.w >= p->w) || (css + ic.w >= p->h ) ) {
			continue;
		}
		if ( (cfs < 0) || (css < 0 ) ) continue;

		bx = add_box(&ic);
		bx->refl = refl;
		bx->cfs = cfs;
		bx->css = css;
		bx->p = p;
		bx->pn = pn;

		/* Which reference profile? */
		bx->rp = 0;//bx->pn;

		if ( meth & INTEGRATION_CENTER ) {
			r = center_and_check_box(&ic, bx, &saturated);
		} else {
			r = check_box(&ic, bx, &saturated);
			bx->offs_fs = 0.0;
			bx->offs_ss = 0.0;
		}
		if ( r ) {
			delete_box(&ic, bx);
			continue;
		}

		if ( saturated ) {
			n_saturated++;
			if ( !(meth & INTEGRATION_SATURATED) ) {
				delete_box(&ic, bx);
				continue;
			}
		}

		get_indices(refl, &h, &k, &l);
		if ( VERBOSITY ) {
			bx->verbose = 1;
		}

		fit_bg(&ic, bx);

		observed_position(&ic, bx, &bx->offs_fs, &bx->offs_ss);

		bx->intensity = tentative_intensity(&ic, bx);
		set_intensity(refl, bx->intensity);

		if ( suitable_reference(&ic, bx) ) {
			add_to_reference_profile(&ic, bx);
		}
	}

	calculate_reference_profiles(&ic);

	for ( i=0; i<ic.n_boxes; i++ ) {

		struct peak_box *bx;

		bx = &ic.boxes[i];
		if ( bx->verbose ) {
			show_reference_profile(&ic, bx->rp);
			STATUS("%f -> ", bx->intensity);
		}
		bx->intensity = fit_intensity(&ic, bx);
		bx->sigma = calc_sigma(&ic, bx);

#if 0
		if ( isnan(bx->intensity) ) {
			signed int h, k, l;
			get_indices(bx->refl, &h, &k, &l);
			STATUS("NaN intensity for %i %i %i !\n", h, k, l);
			STATUS("panel %s\n", image->det->panels[bx->pn].name);
			show_peak_box(&ic, bx);
			show_reference_profile(&ic, bx->rp);
		}
		if ( bx->intensity < 0.0 ) {
			signed int h, k, l;
			get_indices(bx->refl, &h, &k, &l);
			STATUS("Negative intensity (%f) for %i %i %i !\n",
			       bx->intensity, h, k, l);
			STATUS("panel %s\n", image->det->panels[bx->pn].name);
			show_peak_box(&ic, bx);
			show_reference_profile(&ic, bx->rp);
		}
#endif

		if ( bg_ok(bx) ) {

			double pfs, pss;

			set_intensity(bx->refl, bx->intensity);
			set_esd_intensity(bx->refl, bx->sigma);
			set_redundancy(bx->refl, 1);

			/* Update position */
			get_detector_pos(bx->refl, &pfs, &pss);
			pfs += bx->offs_fs;
			pss += bx->offs_ss;
			set_detector_pos(bx->refl, 0.0, pfs, pss);

		}
	}

	free_intcontext(&ic);

	image->num_saturated_peaks = n_saturated;
}


static void estimate_mosaicity(IntegrationMethod meth, Crystal *cr,
                               struct image *image)
{
	int msteps = 50;
	int i;
	const double mest = crystal_get_mosaicity(cr);
	const double mmax = 2.0 * mest;
	RefList *list;

	STATUS("Initial estimate: m = %f\n", mest);

	crystal_set_mosaicity(cr, mmax);
	list = find_intersections(image, cr);
	crystal_set_reflections(cr, list);
	measure_all_intensities(meth, list, image, crystal_get_cell(cr));

	for ( i=1; i<=msteps; i++ ) {

		/* "m" varies from just over zero up to 2x the given estimate */
		Reflection *refl;
		RefListIterator *iter;
		const double m = mmax*((double)i/msteps);
		int n_gained = 0;
		int n_lost = 0;
		double i_gained = 0.0;
		double i_lost = 0.0;

		crystal_set_mosaicity(cr, m);
		update_partialities(cr, PMODEL_SPHERE);

		for ( refl = first_refl(list, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			if ( get_redundancy(refl) == 0 ) {
				if ( get_temp1(refl) > 0.0 ) {
					i_lost += get_intensity(refl);
					n_lost++;
				}
				set_temp1(refl, -1.0);
			} else if ( get_temp1(refl) < 0.0 ) {
				i_gained += get_intensity(refl);
				n_gained++;
				set_temp1(refl, 1.0);
			}
		}

		if ( i > 1 ) {
			STATUS("%.2e %10.2f %4i %10.2f %4i %10.2f\n", m,
			       i_gained, n_gained, i_lost, n_lost,
		               i_gained - i_lost);
		}

	}
}


struct integr_ind
{
	double res;
	Reflection *refl;
};


static int compare_resolution(const void *av, const void *bv)
{
	const struct integr_ind *a = av;
	const struct integr_ind *b = bv;

	return a->res > b->res;
}


static struct integr_ind *sort_reflections(RefList *list, UnitCell *cell,
                                           int *np)
{
	struct integr_ind *il;
	Reflection *refl;
	RefListIterator *iter;
	int i, n;

	*np = 0;  /* For now */

	n = num_reflections(list);
	if ( n == 0 ) return NULL;

	il = calloc(n, sizeof(struct integr_ind));
	if ( il == NULL ) return NULL;

	i = 0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double res;

		if ( get_redundancy(refl) == 0 ) continue;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);

		il[i].res = res;
		il[i].refl = refl;

		i++;
	}

	qsort(il, i, sizeof(struct integr_ind), compare_resolution);

	*np = i;
	return il;
}


static void estimate_resolution(RefList *reflections, Crystal *cr,
                                struct image *image)
{
	struct integr_ind *il;
	int n, i;
	int score = 1000;  /* FIXME */
	int cutoff = 0;
	double limit = 0.0;

	if ( num_reflections(reflections) == 0 ) return;

	il = sort_reflections(reflections, crystal_get_cell(cr), &n);
	if ( il == NULL ) {
		ERROR("Couldn't sort reflections\n");
		return;
	}

	for ( i=0; i<n; i++ ) {

		double intensity, sigma, snr;
		Reflection *refl;

		refl = il[i].refl;
		intensity = get_intensity(refl);
		sigma = get_esd_intensity(refl);

		/* I/sigma(I) cutoff
		 * Rejects reflections below --min-integration-snr, or if the
		 * SNR is clearly silly.  Silly indicates that the intensity
		 * was zero. */
		snr = fabs(intensity)/sigma;

		/* Record intensity and set redundancy to 1 on success */
		if ( cutoff ) {
			set_redundancy(refl, 0);
		}

		if ( snr > 3.0 ) {
			score++;
		} else {
			score--;
		}

		//STATUS("%5.2f A, %5.2f, %i\n", 1e10/il[i].res, snr, score);
		if ( score == 0 ) {
			limit = il[i].res;
			cutoff = 1;
		}

	}

	crystal_set_resolution_limit(cr, limit);

	free(il);
}


static void integrate_prof2d(IntegrationMethod meth, Crystal *cr,
                             struct image *image,
                             double ir_inn, double ir_mid, double ir_out)
{
	RefList *reflections;
	UnitCell *cell;

	cell = crystal_get_cell(cr);

	/* Create initial list of reflections with nominal parameters */
	reflections = find_intersections(image, cr);
	measure_all_intensities(meth, reflections, image, cell);

	/* Find resolution limit of pattern using this list */
	//estimate_resolution(reflections, cr, image);

	//reflist_free(reflections);

	STATUS("Initial resolution estimate = %.2f nm^-1 or %.2f A\n",
	       crystal_get_resolution_limit(cr)/1e9,
	       1e9 / crystal_get_resolution_limit(cr));

	/* Estimate the mosaicity of the crystal using this resolution limit */
	//estimate_mosaicity(cr, image);

	/* Create new list of reflections with refined mosaicity */
	//reflections = find_intersections(image, cr);
	//measure_all_intensities(reflections, image);
	crystal_set_reflections(cr, reflections);

	//estimate_resolution(reflections, cr, image);
}


static void integrate_box(struct intcontext *ic, struct peak_box *bx,
                          double *intensity, double *sigma)
{
	double pk_total;
	int pk_counts;
	double fsct, ssct;
	double bg_tot;
	double bg_tot_sq;
	int bg_counts;
	double bg_mean, bg_var;
	double var;
	double aduph;
	int p, q;

	aduph = bx->p->adu_per_eV * ph_lambda_to_eV(ic->image->lambda);

	/* Measure the background */
	bg_tot = 0.0;
	bg_tot_sq = 0.0;
	bg_counts = 0;
	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		if ( bx->bm[p + ic->w*q] != BM_BG ) continue;

		bi = boxi(ic, bx, p, q);
		bg_tot += bi;
		bg_tot_sq += pow(bi, 2.0);
		bg_counts++;

	}
	}

	bg_mean = bg_tot / bg_counts;
	bg_var = (bg_tot_sq/bg_counts) - pow(bg_mean, 2.0);

	/* Measure the peak */
	pk_total = 0.0;
	pk_counts = 0;
	fsct = 0.0;  ssct = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;

		bi = boxi(ic, bx, p, q);
		pk_counts++;
		pk_total += (bi - bg_mean);
		fsct += (bi-bg_mean)*(p - ic->halfw);
		ssct += (bi-bg_mean)*(q - ic->halfw);

	}
	}

	/* This offset is in addition to the offset from center_and_check_box */
	bx->offs_fs += (double)fsct / pk_total;
	bx->offs_ss += (double)ssct / pk_total;

	var = pk_counts * bg_var;
	var += aduph * pk_total;

	if ( intensity != NULL ) *intensity = pk_total;
	if ( sigma != NULL ) *sigma = sqrt(var);
}


static void integrate_rings(IntegrationMethod meth, Crystal *cr,
                            struct image *image,
                            double ir_inn, double ir_mid, double ir_out)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	UnitCell *cell;
	struct intcontext ic;
	int n_saturated = 0;
	double limit = 0.0;

	list = find_intersections(image, cr);
	if ( list == NULL ) return;

	if ( num_reflections(list) == 0 ) return;

	cell = crystal_get_cell(cr);

	ic.halfw = ir_out;
	ic.image = image;
	ic.k = 1.0/image->lambda;
	ic.cell = cell;
	ic.ir_inn = ir_inn;
	ic.ir_mid = ir_mid;
	ic.ir_out = ir_out;
	if ( init_intcontext(&ic) ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}
	setup_ring_masks(&ic, ir_inn, ir_mid, ir_out);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pfs, pss;
		signed int h, k, l;
		struct peak_box *bx;
		int pn;
		struct panel *p;
		int fid_fs, fid_ss;  /* Center coordinates, rounded,
		                      * in overall data block */
		int cfs, css;  /* Corner coordinates */
		double intensity;
		double sigma;
		int saturated;
		double one_over_d;
		int r;

		set_redundancy(refl, 0);

		get_detector_pos(refl, &pfs, &pss);
		fid_fs = lrint(pfs);
		fid_ss = lrint(pss);
		pn = find_panel_number(image->det, fid_fs, fid_ss);
		p = &image->det->panels[pn];

		cfs = (fid_fs-p->min_fs) - ic.halfw;
		css = (fid_ss-p->min_ss) - ic.halfw;

		if ( (cfs + ic.w >= p->w) || (css + ic.w >= p->h ) ) {
			continue;
		}
		if ( (cfs < 0) || (css < 0 ) ) continue;

		bx = add_box(&ic);
		bx->refl = refl;
		bx->cfs = cfs;
		bx->css = css;
		bx->p = p;
		bx->pn = pn;

		if ( meth & INTEGRATION_CENTER ) {
			r = center_and_check_box(&ic, bx, &saturated);
		} else {
			r = check_box(&ic, bx, &saturated);
			bx->offs_fs = 0.0;
			bx->offs_ss = 0.0;
		}
		if ( r ) {
			delete_box(&ic, bx);
			continue;
		}

		if ( saturated ) {
			n_saturated++;
			if ( !(meth & INTEGRATION_SATURATED) ) {
				delete_box(&ic, bx);
				continue;
			}
		}

		get_indices(refl, &h, &k, &l);
		if ( VERBOSITY ) {
			bx->verbose = 1;
		}

		if ( bx->verbose ) show_peak_box(&ic, bx);

		integrate_box(&ic, bx, &intensity, &sigma);

		/* Record intensity and set redundancy to 1 */
		set_intensity(refl, intensity);
		set_esd_intensity(refl, sigma);
		set_redundancy(refl, 1);

		one_over_d = resolution(cell, h, k, l);
		if ( one_over_d > limit ) limit = one_over_d;

		/* Update position */
		pfs += bx->offs_fs;
		pss += bx->offs_ss;
		set_detector_pos(refl, 0.0, pfs, pss);

	}

	free_intcontext(&ic);

	crystal_set_num_saturated_reflections(cr, n_saturated);
	crystal_set_resolution_limit(cr, limit);
	crystal_set_reflections(cr, list);
}


void integrate_all(struct image *image, IntegrationMethod meth,
                   double ir_inn, double ir_mid, double ir_out)
{
	int i;

	for ( i=0; i<image->n_crystals; i++ ) {

		switch ( meth & INTEGRATION_METHOD_MASK ) {

			case INTEGRATION_NONE :
			return;

			case INTEGRATION_RINGS :
			integrate_rings(meth, image->crystals[i], image,
			                ir_inn, ir_mid, ir_out);
			return;

			case INTEGRATION_PROF2D :
			integrate_prof2d(meth, image->crystals[i], image,
			                 ir_inn, ir_mid, ir_out);
			return;

			default :
			ERROR("Unrecognised integration method %i\n", meth);
			return;

		}

	}
}


IntegrationMethod integration_method(const char *str, int *err)
{
	int n, i;
	char **methods;
	IntegrationMethod meth = INTEGRATION_NONE;

	if ( err != NULL ) *err = 0;
	n = assplode(str, ",-", &methods, ASSPLODE_NONE);

	for ( i=0; i<n; i++ ) {

		if ( strcmp(methods[i], "rings") == 0) {
			meth = INTEGRATION_DEFAULTS_RINGS;

		} else if ( strcmp(methods[i], "prof2d") == 0) {
			meth = INTEGRATION_DEFAULTS_PROF2D;

		} else if ( strcmp(methods[i], "none") == 0) {
			return INTEGRATION_NONE;

		} else if ( strcmp(methods[i], "sat") == 0) {
			meth |= INTEGRATION_SATURATED;

		} else if ( strcmp(methods[i], "nosat") == 0) {
			meth &= ~INTEGRATION_SATURATED;

		} else if ( strcmp(methods[i], "cen") == 0) {
			meth |= INTEGRATION_CENTER;

		} else if ( strcmp(methods[i], "nocen") == 0) {
			meth &= ~INTEGRATION_CENTER;

		} else {
			ERROR("Bad integration method: '%s'\n", str);
			if ( err != NULL ) *err = 1;
			return INTEGRATION_NONE;
		}

		free(methods[i]);

	}
	free(methods);

	return meth;

}
