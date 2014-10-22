/*
 * integration.c
 *
 * Integration of intensities
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
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
#include <unistd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#ifdef HAVE_CURSES_COLOR
#include <ncurses.h>
#endif

#include "reflist.h"
#include "reflist-utils.h"
#include "cell.h"
#include "crystal.h"
#include "cell-utils.h"
#include "geometry.h"
#include "image.h"
#include "peaks.h"
#include "integration.h"


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


static gsl_vector *solve_svd(gsl_vector *v, gsl_matrix *Mp)
{
	gsl_matrix *s_vec;
	gsl_vector *s_val;
	int err, n;
	gsl_vector *shifts;
	gsl_matrix *M;

	n = v->size;
	if ( v->size != Mp->size1 ) return NULL;
	if ( v->size != Mp->size2 ) return NULL;

	M = gsl_matrix_alloc(n, n);
	if ( M == NULL ) return NULL;
	gsl_matrix_memcpy(M, Mp);

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
	gsl_matrix_free(M);

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
	IntegrationMethod meth;

	int halfw;
	int w;
	enum boxmask_val *bm;  /* Box mask */
	struct image *image;
	int **masks;  /* Peak location mask from make_BgMask() */

	struct peak_box *boxes;
	int n_boxes;
	int max_boxes;

	UnitCell *cell;
	double k;

	int n_reference_profiles;
	double **reference_profiles;
	double **reference_den;
	int *n_profiles_in_reference;

	int ir_inn;
	int ir_mid;
	int ir_out;

	int n_saturated;
	int n_implausible;

	IntDiag int_diag;
	signed int int_diag_h;
	signed int int_diag_k;
	signed int int_diag_l;
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

	/* Peak region sums */
	double pks_p2;
	double pks_q2;
	double pks_pq;
	double pks_p;
	double pks_q;
	int m;

	gsl_matrix *bgm;  /* Background estimation matrix */

	/* Measured intensity (tentative, profile fitted or otherwise) */
	double intensity;
	double sigma;
	double J;  /* Profile scaling factor */

	/* Highest (non-ignored) value in peak */
	double peak;

	/* Offsets to final observed position */
	double offs_fs;
	double offs_ss;

	int rp;   /* Reference profile number */

	Reflection *refl;
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


#ifdef HAVE_CURSES_COLOR
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
#endif


#ifdef HAVE_CURSES_COLOR
static void show_reference_profile(struct intcontext *ic, int i)
{
	int q;

	printw("Reference profile number %i (%i contributions):\n", i,
	       ic->n_profiles_in_reference[i]);

	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(ic->bm[p+q*ic->w]);
			printw("%4.0f ", ic->reference_profiles[i][p+ic->w*q]);
			colour_off(ic->bm[p+q*ic->w]);

		}

		printw("\n");
	}
}
#endif


static void show_peak_box(struct intcontext *ic, struct peak_box *bx,
                          int results_pipe)
{
#ifdef HAVE_CURSES_COLOR
	int q;
	signed int h, k, l;
	double fs, ss;

	if ( results_pipe != 0 ) write(results_pipe, "SUSPEND\n", 8);

	initscr();
	clear();
	start_color();
	init_pair(1, COLOR_WHITE, COLOR_BLUE) ;  /* Background */
	init_pair(2, COLOR_WHITE, COLOR_RED);    /* Peak */
	init_pair(3, COLOR_BLACK, COLOR_CYAN);   /* Blackhole */

	get_indices(bx->refl, &h, &k, &l);
	get_detector_pos(bx->refl, &fs, &ss);
	printw("Indices %i %i %i\nPosition fs = %.1f, ss = %.1f\n\n",
	       h, k, l, fs, ss);

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

	printw("\nFitted background (parameters a=%.2f, b=%.2f, c=%.2f)\n",
	       bx->a, bx->b, bx->c);
	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(bx->bm[p+q*ic->w]);
			printw("%5.0f ", bx->a*p + bx->b*q + bx->c);
			colour_off(bx->bm[p+q*ic->w]);

		}

		printw("\n");
	}

	if ( ic->meth & INTEGRATION_PROF2D ) {
		printw("\n");
		show_reference_profile(ic, bx->rp);
	}

	printw("\nIntensity = %.2f +/- %.2f\n", get_intensity(bx->refl),
	                                        get_esd_intensity(bx->refl));

	printw("\n\nPress any key to continue processing...\n\n");

	refresh();
	getch();
	endwin();

	if ( results_pipe != 0 ) write(results_pipe, "RELEASE\n", 8);
#endif
}


static void fit_bg(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	gsl_vector *v;
	gsl_vector *ans;

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

	/* SVD is massive overkill here */
	ans = solve_svd(v, bx->bgm);
	gsl_vector_free(v);

	bx->a = gsl_vector_get(ans, 0);
	bx->b = gsl_vector_get(ans, 1);
	bx->c = gsl_vector_get(ans, 2);

	gsl_vector_free(ans);
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
	int i;

	ic->w = 2*ic->halfw + 1;

	ic->bm = malloc(ic->w * ic->w * sizeof(enum boxmask_val));
	if ( ic->bm == NULL ) {
		ERROR("Failed to allocate box mask.\n");
		return 1;
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

	for ( i=0; i<ic->n_boxes; i++ ) {
		free(ic->boxes[i].bm);
		gsl_matrix_free(ic->boxes[i].bgm);
	}
	free(ic->boxes);

	for ( i=0; i<ic->n_reference_profiles; i++ ) {
		free(ic->reference_profiles[i]);
		free(ic->reference_den[i]);
	}
	free(ic->reference_profiles);
	free(ic->reference_den);
	free(ic->n_profiles_in_reference);
	free(ic->bm);
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

	ic->boxes[idx].bgm = gsl_matrix_calloc(3, 3);
	if ( ic->boxes[idx].bgm == NULL ) {
		ERROR("Failed to initialise matrix.\n");
		return NULL;
	}

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

	free(bx->bm);
	gsl_matrix_free(bx->bgm);

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

	intensity -= bx->a * bx->pks_p;
	intensity -= bx->b * bx->pks_q;
	intensity -= bx->c * bx->m;

	return intensity;
}


static void UNUSED observed_position(struct intcontext *ic, struct peak_box *bx,
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

	num_p += -bx->a*bx->pks_p2 - bx->b*bx->pks_pq - bx->c*bx->pks_p;
	num_q += -bx->a*bx->pks_q2 - bx->b*bx->pks_pq - bx->c*bx->pks_q;
	den += -bx->a*bx->pks_p - bx->b*bx->pks_q - bx->c;

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


static void setup_peak_integrals(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;

	bx->pks_p2 = 0.0;
	bx->pks_q2 = 0.0;
	bx->pks_pq = 0.0;
	bx->pks_p = 0.0;
	bx->pks_q = 0.0;
	bx->m = 0;

	gsl_matrix_set(bx->bgm, 0, 0, 0.0);
	gsl_matrix_set(bx->bgm, 0, 1, 0.0);
	gsl_matrix_set(bx->bgm, 0, 2, 0.0);
	gsl_matrix_set(bx->bgm, 1, 0, 0.0);
	gsl_matrix_set(bx->bgm, 1, 1, 0.0);
	gsl_matrix_set(bx->bgm, 1, 2, 0.0);
	gsl_matrix_set(bx->bgm, 2, 0, 0.0);
	gsl_matrix_set(bx->bgm, 2, 1, 0.0);
	gsl_matrix_set(bx->bgm, 2, 2, 0.0);

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		switch ( bx->bm[p + ic->w*q] ) {

			case BM_IG :
			case BM_BH :
			break;

			case BM_BG :
			addm(bx->bgm, 0, 0, p*p);
			addm(bx->bgm, 0, 1, p*q);
			addm(bx->bgm, 0, 2, p);
			addm(bx->bgm, 1, 0, p*q);
			addm(bx->bgm, 1, 1, q*q);
			addm(bx->bgm, 1, 2, q);
			addm(bx->bgm, 2, 0, p);
			addm(bx->bgm, 2, 1, q);
			addm(bx->bgm, 2, 2, 1);
			break;

			case BM_PK :
			bx->pks_p2 += p*p;
			bx->pks_q2 += q*q;
			bx->pks_pq += p*q;
			bx->pks_p += p;
			bx->pks_q += q;
			bx->m++;
			break;

		}

	}
	}
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
	if ( bx->bm == NULL ) {
		ERROR("Failed to allocate box mask\n");
		return 1;
	}

	cell_get_cartesian(ic->cell,
	                   &adx, &ady, &adz,
	                   &bdx, &bdy, &bdz,
	                   &cdx, &cdy, &cdz);
	get_indices(bx->refl, &hr, &kr, &lr);

	bx->peak = -INFINITY;
	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		int fs, ss;
		double hd, kd, ld;
		signed int h, k, l;
		struct rvec dv;

		fs = bx->cfs + p;
		ss = bx->css + q;

		if ( (fs < 0) || (fs >= bx->p->w)
		  || (ss < 0) || (ss >= bx->p->h) ) {
			return 1;
		}

		bx->bm[p+ic->w*q] = ic->bm[p+ic->w*q];

		if ( ic->image->bad[bx->pn][fs + bx->p->w*ss] ) {
			bx->bm[p+ic->w*q] = BM_BH;
		}

		/* If this is a background pixel, it shouldn't contain any
		 * pixels which are in the peak region of ANY reflection */
		if ( ic->masks != NULL ) {

			switch ( bx->bm[p+ic->w*q] ) {

				case BM_BG:
				case BM_IG:
				if ( ic->masks[bx->pn][fs + bx->p->w*ss] > 0 ) {
					bx->bm[p+ic->w*q] = BM_BH;
				}
				break;

				case BM_PK:
				if ( ic->masks[bx->pn][fs + bx->p->w*ss] > 1 ) {
					bx->bm[p+ic->w*q] = BM_BH;
				}
				break;

				case BM_BH:
				break;

			}
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

		if ( (bx->bm[p+ic->w*q] != BM_IG)
		  && (bx->bm[p+ic->w*q] != BM_BH)
		  && (boxi(ic, bx, p, q) > bx->peak) ) {
			bx->peak = boxi(ic, bx, p, q);
		}

		if ( bx->bm[p+ic->w*q] == BM_PK ) n_pk++;
		if ( bx->bm[p+ic->w*q] == BM_BG ) n_bg++;

	}
	}

	setup_peak_integrals(ic, bx);

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
	int i, nstrong;

	bx->offs_fs = 0.0;
	bx->offs_ss = 0.0;

	if ( check_box(ic, bx, sat) ) return 1;
	fit_bg(ic, bx);

	nstrong = 0;
	for ( i=0; i<10; i++ ) {

		int p, q;
		double max = -INFINITY;
		int t_offs_fs = 0;
		int t_offs_ss = 0;
		int ifs = 0;
		int iss = 0;

		for ( q=0; q<ic->w; q++ ) {
		for ( p=0; p<ic->w; p++ ) {

			double bi, bg;

			if ( bx->bm[p + ic->w*q] == BM_BH ) continue;

			bi = boxi(ic, bx, p, q);
			bg = bx->a*p + bx->b*q + bx->c;

			if ( bi <= 3.0*bg ) continue;
			nstrong++;

			if ( bi > max ) {
				max = bi;
				ifs = p - ic->halfw;
				iss = q - ic->halfw;
			}

		}
		}

		/* We require at least two bright pixels in the peak region,
		 * otherwise we might just be centering on a pixel which is
		 * just at the tail of the distribution */
		if ( nstrong < 2 ) return 0;

		bx->offs_fs += ifs;
		bx->offs_ss += iss;
		bx->cfs += ifs;
		bx->css += iss;
		t_offs_fs += ifs;
		t_offs_ss += iss;

		free(bx->bm);
		if ( check_box(ic, bx, sat) ) {
			return 1;
		}

		if ( t_offs_fs*t_offs_fs + t_offs_ss*t_offs_ss > ic->w*ic->w ) {
			return 1;
		}

		fit_bg(ic, bx);

		if ( (ifs==0) && (iss==0) ) break;

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


static void mean_var_area(struct intcontext *ic, struct peak_box *bx,
                          enum boxmask_val v, double *pmean, double *pvar)
{
	int p, q;
	double sum = 0.0;
	double var = 0.0;
	int n = 0;
	double mean;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {
		if ( bx->bm[p + ic->w*q] != v ) continue;
		sum += boxi(ic, bx, p, q);
		n++;
	}
	}
	mean = sum/n;

	n = 0;
	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {
		if ( bx->bm[p + ic->w*q] != v ) continue;
		var += pow(boxi(ic, bx, p, q) - mean, 2.0);
		n++;
	}
	}
	var = var/n;

	*pmean = mean;
	*pvar = var;
}


static double bg_under_peak(struct intcontext *ic, struct peak_box *bx)
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

	height_ok = max > 10.0 * bg_under_peak(ic, bx);

	return bg_ok(bx) && height_ok;
}


static void add_to_rg_matrix(struct intcontext *ic, struct panel *p,
                             gsl_matrix *M, gsl_vector *vx, gsl_vector *vy,
                             int *pn)
{
	int i;

	for ( i=0; i<ic->n_boxes; i++ ) {

		double x, y, w;
		double fs, ss;
		double offs_x, offs_y;

		if ( ic->boxes[i].p != p ) continue;

		fs = ic->boxes[i].cfs + ic->halfw;
		ss = ic->boxes[i].css + ic->halfw;

		twod_mapping(fs, ss, &x, &y, ic->boxes[i].p);

		w = ic->boxes[i].intensity;

		addm(M, 0, 0, w*fs*fs);
		addm(M, 0, 1, w*fs*ss);
		addm(M, 0, 2, w*fs);
		addm(M, 1, 0, w*fs*ss);
		addm(M, 1, 1, w*ss*ss);
		addm(M, 1, 2, w*ss);
		addm(M, 2, 0, w*fs);
		addm(M, 2, 1, w*ss);
		addm(M, 2, 2, w);

		/* Offsets in lab coordinate system */
		offs_x = p->fsx * ic->boxes[i].offs_fs
		       + p->ssx * ic->boxes[i].offs_ss;

		offs_y = p->fsy * ic->boxes[i].offs_fs
		       + p->ssy * ic->boxes[i].offs_ss;

		addv(vx, 0, w*offs_x*fs);
		addv(vx, 1, w*offs_x*ss);
		addv(vx, 2, w*offs_x);

		addv(vy, 0, w*offs_y*fs);
		addv(vy, 1, w*offs_y*ss);
		addv(vy, 2, w*offs_y);

		(*pn)++;
	}
}


static void UNUSED refine_rigid_groups(struct intcontext *ic)
{
	int i;

	for ( i=0; i<ic->image->det->n_rigid_groups; i++ ) {

		struct rigid_group *rg;
		gsl_matrix *M;
		gsl_vector *vx;
		gsl_vector *vy;
		gsl_vector *dq1;
		gsl_vector *dq2;
		int j;
		int n;

		M = gsl_matrix_calloc(3, 3);
		if ( M == NULL ) {
			ERROR("Failed to allocate matrix\n");
			return;
		}

		vx = gsl_vector_calloc(3);
		if ( vx == NULL ) {
			ERROR("Failed to allocate vector\n");
			return;
		}

		vy = gsl_vector_calloc(3);
		if ( vy == NULL ) {
			ERROR("Failed to allocate vector\n");
			return;
		}

		rg = ic->image->det->rigid_groups[i];

		n = 0;
		for ( j=0; j<rg->n_panels; j++ ) {
			add_to_rg_matrix(ic, rg->panels[j], M, vx, vy, &n);
		}

		if ( n > 10 ) {

			dq1 = solve_svd(vx, M);
			dq2 = solve_svd(vy, M);

			rg->d_fsx = gsl_vector_get(dq1, 0);
			rg->d_ssx = gsl_vector_get(dq1, 1);
			rg->d_cnx = gsl_vector_get(dq1, 2);
			rg->d_fsy = gsl_vector_get(dq2, 0);
			rg->d_ssy = gsl_vector_get(dq2, 1);
			rg->d_cny = gsl_vector_get(dq2, 2);
			rg->have_deltas = 1;

			gsl_vector_free(dq1);
			gsl_vector_free(dq2);

		} else {

			rg->d_fsx = 0.0;
			rg->d_ssx = 0.0;
			rg->d_cnx = 0.0;
			rg->d_fsy = 0.0;
			rg->d_ssy = 0.0;
			rg->d_cny = 0.0;
			rg->have_deltas = 0;

		}

		gsl_vector_free(vx);
		gsl_vector_free(vy);
		gsl_matrix_free(M);

	}
}


static int get_int_diag(struct intcontext *ic, Reflection *refl)
{
	if ( ic->int_diag == INTDIAG_NONE ) return 0;

	if ( ic->int_diag == INTDIAG_ALL ) return 1;

	if ( ic->int_diag == INTDIAG_RANDOM ) {
		return random() < RAND_MAX/100;
	}

	if ( ic->int_diag == INTDIAG_NEGATIVE ) {
		double i, sigi;
		i = get_intensity(refl);
		sigi = get_esd_intensity(refl);
		return i < -3.0*sigi;
	}

	if ( ic->int_diag == INTDIAG_IMPLAUSIBLE ) {
		double i, sigi;
		i = get_intensity(refl);
		sigi = get_esd_intensity(refl);
		return i < -5.0*sigi;
	}

	if ( ic->int_diag == INTDIAG_STRONG ) {
		double i, sigi;
		i = get_intensity(refl);
		sigi = get_esd_intensity(refl);
		return i > 3.0*sigi;
	}

	if ( ic->int_diag == INTDIAG_INDICES ) {
		signed int h, k, l;
		get_indices(refl, &h, &k, &l);
		if ( ic->int_diag_h != h ) return 0;
		if ( ic->int_diag_k != k ) return 0;
		if ( ic->int_diag_l != l ) return 0;
		return 1;
	}

	return 0;
}



static void integrate_prof2d_once(struct intcontext *ic, struct peak_box *bx,
                                  int results_pipe)
{
	bx->intensity = fit_intensity(ic, bx);
	bx->sigma = calc_sigma(ic, bx);

	if ( bg_ok(bx) ) {

		double pfs, pss;
		double bgmean;
		double sig2_bg;  /* unused */

		mean_var_area(ic, bx, BM_BG, &bgmean, &sig2_bg);

		set_intensity(bx->refl, bx->intensity);
		set_esd_intensity(bx->refl, bx->sigma);
		set_peak(bx->refl, bx->peak);
		set_mean_bg(bx->refl, bgmean);
		set_redundancy(bx->refl, 1);

		/* Update position */
		get_detector_pos(bx->refl, &pfs, &pss);
		pfs += bx->offs_fs;
		pss += bx->offs_ss;
		set_detector_pos(bx->refl, 0.0, pfs, pss);

		if ( bx->intensity < -5.0*bx->sigma ) {
			ic->n_implausible++;
			set_redundancy(bx->refl, 0);
		}

		if ( get_int_diag(ic, bx->refl) ) show_peak_box(ic, bx,
		                                                results_pipe);

	} else {

		set_redundancy(bx->refl, 0);

	}
}


static void setup_profile_boxes(struct intcontext *ic, RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pfs, pss;
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

		/* Explicit truncation of digits after the decimal point.
		 * This is actually the correct thing to do here, not
		 * e.g. lrint().  pfs/pss is the position of the spot, measured
		 * in numbers of pixels, from the panel corner (not the center
		 * of the first pixel).  So any coordinate from 2.0 to 2.9999
		 * belongs to pixel index 2. */
		fid_fs = pfs;
		fid_ss = pss;
		pn = find_panel_number(ic->image->det, fid_fs, fid_ss);
		p = &ic->image->det->panels[pn];

		cfs = (fid_fs-p->min_fs) - ic->halfw;
		css = (fid_ss-p->min_ss) - ic->halfw;

		bx = add_box(ic);
		bx->refl = refl;
		bx->cfs = cfs;
		bx->css = css;
		bx->p = p;
		bx->pn = pn;

		/* Which reference profile? */
		bx->rp = 0;//bx->pn;

		if ( ic->meth & INTEGRATION_CENTER ) {
			r = center_and_check_box(ic, bx, &saturated);
		} else {
			r = check_box(ic, bx, &saturated);
			bx->offs_fs = 0.0;
			bx->offs_ss = 0.0;
		}
		if ( r ) {
			delete_box(ic, bx);
			continue;
		}

		if ( saturated ) {
			ic->n_saturated++;
			if ( !(ic->meth & INTEGRATION_SATURATED) ) {
				delete_box(ic, bx);
				continue;
			}
		}

		fit_bg(ic, bx);

		bx->intensity = tentative_intensity(ic, bx);
		set_intensity(refl, bx->intensity);

		if ( suitable_reference(ic, bx) ) {
			add_to_reference_profile(ic, bx);
		}
	}
}


static void integrate_prof2d(IntegrationMethod meth,
                             Crystal *cr, struct image *image, IntDiag int_diag,
                             signed int idh, signed int idk, signed int idl,
                             double ir_inn, double ir_mid, double ir_out,
                             int results_pipe, int **masks)
{
	RefList *list;
	UnitCell *cell;
	struct intcontext ic;
	int i;
	int n_saturated = 0;

	list = crystal_get_reflections(cr);
	cell = crystal_get_cell(cr);

	ic.halfw = ir_out;
	ic.image = image;
	ic.k = 1.0/image->lambda;
	ic.meth = meth;
	ic.n_saturated = 0;
	ic.n_implausible = 0;
	ic.cell = cell;
	ic.int_diag = int_diag;
	ic.int_diag_h = idh;
	ic.int_diag_k = idk;
	ic.int_diag_l = idl;
	ic.masks = masks;
	if ( init_intcontext(&ic) ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}

	setup_ring_masks(&ic, ir_inn, ir_mid, ir_out);
	setup_profile_boxes(&ic, list);
	calculate_reference_profiles(&ic);

	for ( i=0; i<ic.n_reference_profiles; i++ ) {
		if ( ic.n_profiles_in_reference[i] == 0 ) {
			ERROR("Reference profile %i has no contributions.\n",
			      i);
			free_intcontext(&ic);
			return;
		}
	}

	for ( i=0; i<ic.n_boxes; i++ ) {
		struct peak_box *bx;
		bx = &ic.boxes[i];
		integrate_prof2d_once(&ic, bx, results_pipe);
	}

	//refine_rigid_groups(&ic);

	free_intcontext(&ic);

	image->num_saturated_peaks = n_saturated;
}


static void integrate_rings_once(Reflection *refl, struct image *image,
                                 struct intcontext *ic, UnitCell *cell,
                                 int results_pipe)
{
	double pfs, pss;
	struct peak_box *bx;
	int pn;
	struct panel *p;
	int fid_fs, fid_ss;  /* Center coordinates, rounded,
	                      * in overall data block */
	int cfs, css;  /* Corner coordinates */
	double intensity;
	double sigma;
	int saturated;
	int r;
	double bgmean, sig2_bg, sig2_poisson, aduph;

	set_redundancy(refl, 0);

	get_detector_pos(refl, &pfs, &pss);

	/* Explicit truncation of digits after the decimal point.
	 * This is actually the correct thing to do here, not
	 * e.g. lrint().  pfs/pss is the position of the spot, measured
	 * in numbers of pixels, from the panel corner (not the center
	 * of the first pixel).  So any coordinate from 2.0 to 2.9999
	 * belongs to pixel index 2. */
	fid_fs = pfs;
	fid_ss = pss;
	pn = find_panel_number(image->det, fid_fs, fid_ss);
	p = &image->det->panels[pn];

	cfs = (fid_fs-p->min_fs) - ic->halfw;
	css = (fid_ss-p->min_ss) - ic->halfw;

	bx = add_box(ic);
	bx->refl = refl;
	bx->cfs = cfs;
	bx->css = css;
	bx->p = p;
	bx->pn = pn;

	if ( ic->meth & INTEGRATION_CENTER ) {
		r = center_and_check_box(ic, bx, &saturated);
	} else {
		r = check_box(ic, bx, &saturated);
		if ( !r ) {
			fit_bg(ic, bx);
		}
		bx->offs_fs = 0.0;
		bx->offs_ss = 0.0;
	}
	if ( r ) {
		delete_box(ic, bx);
		return;
	}

	if ( saturated ) {
		ic->n_saturated++;
		if ( !(ic->meth & INTEGRATION_SATURATED) ) {
			delete_box(ic, bx);
			return;
		}
	}

	intensity = tentative_intensity(ic, bx);
	mean_var_area(ic, bx, BM_BG, &bgmean, &sig2_bg);

	aduph = bx->p->adu_per_eV * ph_lambda_to_eV(ic->image->lambda);
	sig2_poisson = aduph * intensity;

	/* If intensity is within one photon of nothing, set the Poisson
	 * error to be one photon */
	if ( fabs(intensity / aduph) < 1.0 ) {
		sig2_poisson = aduph;
	}

	/* If intensity is negative by more than one photon, assume that
	 * the peak is in the background and add a Poisson error of
	 * appropriate size */
	if ( intensity < -aduph ) {
		sig2_poisson = -aduph*intensity;
	}

	sigma = sqrt(sig2_poisson + bx->m*sig2_bg);

	/* Record intensity and set redundancy to 1 */
	bx->intensity = intensity;
	set_intensity(refl, intensity);
	set_esd_intensity(refl, sigma);
	set_redundancy(refl, 1);
	set_mean_bg(refl, bgmean);
	set_peak(refl, bx->peak);

	/* Update position */
	pfs += bx->offs_fs;
	pss += bx->offs_ss;
	set_detector_pos(refl, 0.0, pfs, pss);

	if ( get_int_diag(ic, refl) ) show_peak_box(ic, bx, results_pipe);

	if ( intensity < -5.0*sigma ) {
		ic->n_implausible++;
		set_redundancy(refl, 0);
	}
}


static int compare_double(const void *av, const void *bv)
{
	double a = *(double *)av;
	double b = *(double *)bv;
	if ( a > b ) return 1;
	if ( a < b ) return -1;
	return 0;
}


static double estimate_resolution(UnitCell *cell, ImageFeatureList *flist)
{
	int i;
	const double min_dist = 0.25;
	double max_res = 0.0;
	double *acc;
	int n_acc = 0;
	int max_acc = 1024;
	int n;

	acc = malloc(max_acc*sizeof(double));
	if ( acc == NULL ) {
		ERROR("Allocation failed during estimate_resolution!\n");
		return INFINITY;
	}

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
			double res = 2.0*resolution(cell, h, k, l); /* 1/d */
			acc[n_acc++] = res;
			if ( n_acc == max_acc ) {
				max_acc += 1024;
				acc = realloc(acc, max_acc*sizeof(double));
				if ( acc == NULL ) {
					ERROR("Allocation failed during"
					      " estimate_resolution!\n");
					return INFINITY;
				}
			}
		}

	}

	if ( n_acc < 3 ) {
		STATUS("WARNING: Too few peaks to estimate resolution.\n");
		return 0.0;
	}

	/* Slightly horrible outlier removal */
	qsort(acc, n_acc, sizeof(double), compare_double);
	n = n_acc/50;
	if ( n < 2 ) n = 2;
	max_res = acc[(n_acc-1)-n];

	free(acc);
	return max_res;
}


static void integrate_rings(IntegrationMethod meth,
                            Crystal *cr, struct image *image, IntDiag int_diag,
                            signed int idh, signed int idk, signed int idl,
                            double ir_inn, double ir_mid, double ir_out,
                            int results_pipe, int **masks)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	UnitCell *cell;
	struct intcontext ic;

	list = crystal_get_reflections(cr);
	cell = crystal_get_cell(cr);

	ic.halfw = ir_out;
	ic.image = image;
	ic.k = 1.0/image->lambda;
	ic.n_saturated = 0;
	ic.n_implausible = 0;
	ic.cell = cell;
	ic.ir_inn = ir_inn;
	ic.ir_mid = ir_mid;
	ic.ir_out = ir_out;
	ic.int_diag = int_diag;
	ic.int_diag_h = idh;
	ic.int_diag_k = idk;
	ic.int_diag_l = idl;
	ic.meth = meth;
	ic.masks = masks;
	if ( init_intcontext(&ic) ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}
	setup_ring_masks(&ic, ir_inn, ir_mid, ir_out);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		integrate_rings_once(refl, image, &ic, cell, results_pipe);
	}

	//refine_rigid_groups(&ic);

	free_intcontext(&ic);

	crystal_set_num_saturated_reflections(cr, ic.n_saturated);

	if ( ic.n_implausible ) {
		STATUS("Warning: %i implausibly negative reflection%s.\n",
		       ic.n_implausible, ic.n_implausible>1?"s":"");
	}
}


static void apply_resolution_cutoff(Crystal *cr, double res)
{
	Reflection *refl;
	RefListIterator *iter;
	UnitCell *cell;

	cell = crystal_get_cell(cr);

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		get_indices(refl, &h, &k, &l);
		if ( 2.0*resolution(cell, h, k, l) > res ) {
			set_redundancy(refl, 0);
		}
	}
}


void integrate_all_4(struct image *image, IntegrationMethod meth,
                     PartialityModel pmodel, double push_res,
                     double ir_inn, double ir_mid, double ir_out,
                     IntDiag int_diag,
                     signed int idh, signed int idk, signed int idl,
                     int results_pipe)
{
	int i;
	int *masks[image->det->n_panels];

	/* Predict all reflections */
	for ( i=0; i<image->n_crystals; i++ ) {
		RefList *list;
		list = find_intersections(image, image->crystals[i], pmodel);
		crystal_set_reflections(image->crystals[i], list);
	}

	for ( i=0; i<image->det->n_panels; i++ ) {
		masks[i] = make_BgMask(image, &image->det->panels[i], ir_inn);
	}

	for ( i=0; i<image->n_crystals; i++ ) {

		double res = INFINITY;
		Crystal *cr = image->crystals[i];

		switch ( meth & INTEGRATION_METHOD_MASK ) {

			case INTEGRATION_NONE :
			break;

			case INTEGRATION_RINGS :
			integrate_rings(meth, cr, image,
			                int_diag, idh, idk, idl,
			                ir_inn, ir_mid, ir_out,
			                results_pipe, masks);
			res = estimate_resolution(crystal_get_cell(cr),
			                          image->features);
			break;

			case INTEGRATION_PROF2D :
			integrate_prof2d(meth, cr, image,
			                 int_diag, idh, idk, idl,
			                 ir_inn, ir_mid, ir_out,
			                 results_pipe, masks);
			res = estimate_resolution(crystal_get_cell(cr),
			                          image->features);
			break;

			default :
			ERROR("Unrecognised integration method %i\n", meth);
			break;

		}

		crystal_set_resolution_limit(cr, res);
		if ( meth & INTEGRATION_RESCUT ) {
			apply_resolution_cutoff(cr, res+push_res);
		}

	}

	for ( i=0; i<image->det->n_panels; i++ ) {
		free(masks[i]);
	}
}


void integrate_all_3(struct image *image, IntegrationMethod meth,
                     PartialityModel pmodel, double push_res,
                     double ir_inn, double ir_mid, double ir_out,
                     IntDiag int_diag,
                     signed int idh, signed int idk, signed int idl)
{
	integrate_all_4(image, meth, pmodel, 0.0, ir_inn, ir_mid, ir_out,
	                int_diag, idh, idk, idl, 0);
}


void integrate_all_2(struct image *image, IntegrationMethod meth,
                     double push_res,
                     double ir_inn, double ir_mid, double ir_out,
                     IntDiag int_diag,
                     signed int idh, signed int idk, signed int idl)
{
	integrate_all_3(image, meth, PMODEL_SCSPHERE, 0.0,
	                ir_inn, ir_mid, ir_out,
	                int_diag, idh, idk, idl);
}


void integrate_all(struct image *image, IntegrationMethod meth,
                   double ir_inn, double ir_mid, double ir_out,
                   IntDiag int_diag,
                   signed int idh, signed int idk, signed int idl)
{
	integrate_all_2(image, meth, 0.0, ir_inn, ir_mid, ir_out,
	                int_diag, idh, idk, idl);
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

		} else if ( strcmp(methods[i], "rescut") == 0) {
			meth |= INTEGRATION_RESCUT;

		} else if ( strcmp(methods[i], "norescut") == 0) {
			meth &= ~INTEGRATION_RESCUT;

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

