/*
 * integration.c
 *
 * Integration of intensities
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
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
#include <unistd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "reflist.h"
#include "reflist-utils.h"
#include "cell.h"
#include "crystal.h"
#include "cell-utils.h"
#include "geometry.h"
#include "image.h"
#include "peaks.h"
#include "integration.h"
#include "detgeom.h"


/** \file integration.h */


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
	struct detgeom_panel *p;  /* The panel itself */

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


static void colour_on(enum boxmask_val b)
{
	switch ( b ) {

		case BM_BG :
		printf("\e[44m\e[37m");
		break;

		case BM_PK :
		printf("\e[41m\e[37m");
		break;

		case BM_BH :
		printf("\e[46m\e[30m");
		break;

		default:
		break;

	}
}


static void colour_off(enum boxmask_val b)
{
	switch ( b ) {

		case BM_BG :
		printf("\e[49m\e[39m");
		break;

		case BM_PK :
		printf("\e[49m\e[39m");
		break;

		case BM_BH :
		printf("\e[49m\e[39m");
		break;

		default:
		break;

	}
}


static void show_reference_profile(struct intcontext *ic, int i)
{
	int q;

	printf("Reference profile number %i (%i contributions):\n", i,
	       ic->n_profiles_in_reference[i]);

	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(ic->bm[p+q*ic->w]);
			printf("%4.0f ", ic->reference_profiles[i][p+ic->w*q]);
			colour_off(ic->bm[p+q*ic->w]);

		}

		printf("\n");
	}
}


static void show_peak_box(struct intcontext *ic, struct peak_box *bx,
                          pthread_mutex_t *term_lock)
{
	int q;
	signed int h, k, l;
	double fs, ss;

	if ( term_lock != NULL ) pthread_mutex_lock(term_lock);

	get_indices(bx->refl, &h, &k, &l);
	get_detector_pos(bx->refl, &fs, &ss);
	printf("-------- Start of integration diagnostics\n");
	printf("Indices %i %i %i\nPanel %s\nPosition fs = %.1f, ss = %.1f\n\n",
	       h, k, l, bx->p->name, fs, ss);

	printf("Pixel values:\n");
	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(bx->bm[p+q*ic->w]);
			printf("%5.0f ", boxi(ic, bx, p, q));
			colour_off(bx->bm[p+q*ic->w]);

		}

		printf("\n");
	}

	printf("\nFitted background (parameters a=%.2f, b=%.2f, c=%.2f)\n",
	       bx->a, bx->b, bx->c);
	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {

			colour_on(bx->bm[p+q*ic->w]);
			printf("%5.0f ", bx->a*p + bx->b*q + bx->c);
			colour_off(bx->bm[p+q*ic->w]);

		}

		printf("\n");
	}

	if ( ic->meth & INTEGRATION_PROF2D ) {
		printf("\n");
		show_reference_profile(ic, bx->rp);
	}

	printf("\nIntensity = %.2f +/- %.2f\n", get_intensity(bx->refl),
	                                        get_esd_intensity(bx->refl));
	printf("-------- End of integration diagnostics\n");

	if ( term_lock != NULL ) pthread_mutex_unlock(term_lock);
}


static void fit_gradient_bg(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	gsl_vector *v;
	gsl_vector *ans;

	v = gsl_vector_calloc(3);

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		if ( bx->bm[p + ic->w*q] == BM_BG ) {

			double bi;
			bi = boxi(ic, bx, p, q);

			addv(v, 0, bi*p);
			addv(v, 1, bi*q);
			addv(v, 2, bi);

		}

	}
	}

	/* SVD is massive overkill here, but the routine is right there. */
	ans = solve_svd(v, bx->bgm, NULL, 0);
	gsl_vector_free(v);

	bx->a = gsl_vector_get(ans, 0);
	bx->b = gsl_vector_get(ans, 1);
	bx->c = gsl_vector_get(ans, 2);

	gsl_vector_free(ans);
}


static void fit_bg(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double tbg = 0.0;
	int n = 0;

	if ( ic->meth & INTEGRATION_GRADIENTBG ) {
		fit_gradient_bg(ic, bx);
		return;
	}

	/* else do a flat background */
	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {
		if ( bx->bm[p + ic->w*q] == BM_BG ) {
			tbg += boxi(ic, bx, p, q);
			n++;
		}
	}
	}

	bx->a = 0.0;
	bx->b = 0.0;
	bx->c = tbg / n;
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


static void setup_ring_masks(struct intcontext *ic,
                             double ir_inn,
                             double ir_mid,
                             double ir_out)
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


void intcontext_set_diag(struct intcontext *ic,
                         IntDiag int_diag,
                         signed int idh,
                         signed int idk,
                         signed int idl)
{
	ic->int_diag = int_diag;
	ic->int_diag_h = idh;
	ic->int_diag_k = idk;
	ic->int_diag_l = idl;
}


struct intcontext *intcontext_new(struct image *image,
                                  UnitCell *cell,
                                  IntegrationMethod meth,
                                  int ir_inn, int ir_mid, int ir_out,
                                  int **masks)
{
	int i;
	struct intcontext *ic;

	ic = malloc(sizeof(struct intcontext));
	if ( ic == NULL ) return NULL;

	ic->halfw = ir_out;
	ic->image = image;
	ic->k = 1.0/image->lambda;
	ic->meth = meth;
	ic->n_saturated = 0;
	ic->n_implausible = 0;
	ic->cell = cell;
	ic->masks = masks;
	ic->int_diag = INTDIAG_NONE;
	ic->w = 2*ic->halfw + 1;

	ic->bm = malloc(ic->w * ic->w * sizeof(enum boxmask_val));
	if ( ic->bm == NULL ) {
		ERROR("Failed to allocate box mask.\n");
		free(ic);
		return NULL;
	}

	/* How many reference profiles? */
	ic->n_reference_profiles = 1;
	ic->reference_profiles = calloc(ic->n_reference_profiles,
	                                sizeof(double *));
	if ( ic->reference_profiles == NULL ) {
		free(ic);
		return NULL;
	}
	ic->reference_den = calloc(ic->n_reference_profiles, sizeof(double *));
	if ( ic->reference_den == NULL ) {
		free(ic);
		return NULL;
	}
	ic->n_profiles_in_reference = calloc(ic->n_reference_profiles,
	                                     sizeof(int));
	if ( ic->n_profiles_in_reference == NULL ) {
		free(ic);
		return NULL;
	}
	for ( i=0; i<ic->n_reference_profiles; i++ ) {
		ic->reference_profiles[i] = malloc(ic->w*ic->w*sizeof(double));
		if ( ic->reference_profiles[i] == NULL ) return NULL;
		ic->reference_den[i] = malloc(ic->w*ic->w*sizeof(double));
		if ( ic->reference_den[i] == NULL ) return NULL;
	}
	zero_profiles(ic);

	ic->boxes = NULL;
	ic->n_boxes = 0;
	ic->max_boxes = 0;
	if ( alloc_boxes(ic, 32) ) {
		free(ic);
		return NULL;
	}

	setup_ring_masks(ic, ir_inn, ir_mid, ir_out);

	return ic;
}


void intcontext_free(struct intcontext *ic)
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

	bx->bm = malloc(ic->w*ic->w*sizeof(enum boxmask_val));
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
		float lsat;

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

		/* Per-pixel saturation value */
		if ( ic->image->sat != NULL ) {
			lsat = ic->image->sat[bx->pn][fs + bx->p->w*ss];
		} else {
			lsat = INFINITY;
		}
		if ( (bx->bm[p+ic->w*q] != BM_IG)
		  && (bx->bm[p+ic->w*q] != BM_BH)
		  && ((boxi(ic, bx, p, q) > bx->p->max_adu)
		  || (boxi(ic, bx, p, q) > lsat)) )
		{
			if ( sat != NULL ) *sat = 1;
		}

		/* Find brightest pixel */
		if ( (bx->bm[p+ic->w*q] != BM_IG)
		  && (bx->bm[p+ic->w*q] != BM_BH)
		  && (boxi(ic, bx, p, q) > bx->peak) )
		{
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


static double peak_height(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double max = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		double bi;

		if ( bx->bm[p + ic->w*q] != BM_PK ) continue;

		bi = boxi(ic, bx, p, q);
		if ( bi > max ) max = bi;

	}
	}

	return max;
}


static int bg_ok(struct intcontext *ic, struct peak_box *bx)
{
	double max_grad;

	max_grad = fabs((peak_height(ic, bx) - bg_under_peak(ic, bx))) / 10.0;

	if ( (fabs(bx->a) > max_grad) || (fabs(bx->b) > max_grad) ) {
		return 0;
	} else {
		return 1;
	}
}


static int suitable_reference(struct intcontext *ic, struct peak_box *bx)
{
	int height_ok;
	double max = peak_height(ic, bx);
	height_ok = max > 10.0 * bg_under_peak(ic, bx);
	return bg_ok(ic, bx) && height_ok;
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
                                  pthread_mutex_t *term_lock)
{
	bx->intensity = fit_intensity(ic, bx);
	bx->sigma = calc_sigma(ic, bx);

	if ( bg_ok(ic, bx) ) {

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
		set_detector_pos(bx->refl, pfs, pss);

		if ( bx->intensity < -5.0*bx->sigma ) {
			ic->n_implausible++;
			set_redundancy(bx->refl, 0);
		}

		if ( get_int_diag(ic, bx->refl) ) {
			show_peak_box(ic, bx, term_lock);
		}

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
		int fid_fs, fid_ss;  /* Center coordinates, rounded,
		                      * in overall data block */
		int cfs, css;  /* Corner coordinates */
		int saturated;
		int r;

		set_redundancy(refl, 0);

		get_detector_pos(refl, &pfs, &pss);
		pn = get_panel_number(refl);

		/* Explicit truncation of digits after the decimal point.
		 * This is actually the correct thing to do here, not
		 * e.g. lrint().  pfs/pss is the position of the spot, measured
		 * in numbers of pixels, from the panel corner (not the center
		 * of the first pixel).  So any coordinate from 2.0 to 2.9999
		 * belongs to pixel index 2. */
		fid_fs = pfs;
		fid_ss = pss;

		cfs = fid_fs - ic->halfw;
		css = fid_ss - ic->halfw;

		/* Add the box */
		bx = add_box(ic);
		bx->refl = refl;
		bx->cfs = cfs;
		bx->css = css;
		bx->p = &ic->image->detgeom->panels[pn];
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


void integrate_prof2d(IntegrationMethod meth,
                      Crystal *cr, struct image *image, IntDiag int_diag,
                      signed int idh, signed int idk, signed int idl,
                      double ir_inn, double ir_mid, double ir_out,
                      pthread_mutex_t *term_lock, int **masks)
{
	RefList *list;
	UnitCell *cell;
	struct intcontext *ic;
	int i;

	list = crystal_get_reflections(cr);
	cell = crystal_get_cell(cr);

	ic = intcontext_new(image, cell, meth,
	                    ir_inn, ir_mid, ir_out,
	                    masks);
	if ( ic == NULL ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}

	intcontext_set_diag(ic, int_diag, idh, idk, idl);
	setup_profile_boxes(ic, list);
	calculate_reference_profiles(ic);

	for ( i=0; i<ic->n_reference_profiles; i++ ) {
		if ( ic->n_profiles_in_reference[i] == 0 ) {
			ERROR("Reference profile %i has no contributions.\n",
			      i);
			intcontext_free(ic);
			return;
		}
	}

	for ( i=0; i<ic->n_boxes; i++ ) {
		struct peak_box *bx;
		bx = &ic->boxes[i];
		integrate_prof2d_once(ic, bx, term_lock);
	}

	intcontext_free(ic);
}


int integrate_rings_once(Reflection *refl,
                         struct intcontext *ic,
                         pthread_mutex_t *term_lock)
{
	double pfs, pss;
	struct peak_box *bx;
	int pn;
	int fid_fs, fid_ss;  /* Center coordinates, rounded,
	                      * in overall data block */
	int cfs, css;  /* Corner coordinates */
	double intensity;
	double sigma;
	int saturated;
	int r;
	double bgmean, sig2_bg, sig2_poisson, aduph, bias;

	set_redundancy(refl, 0);

	get_detector_pos(refl, &pfs, &pss);
	pn = get_panel_number(refl);

	/* Explicit truncation of digits after the decimal point.
	 * This is actually the correct thing to do here, not
	 * e.g. lrint().  pfs/pss is the position of the spot, measured
	 * in numbers of pixels, from the panel corner (not the center
	 * of the first pixel).  So any coordinate from 2.0 to 2.9999
	 * belongs to pixel index 2. */
	fid_fs = pfs;
	fid_ss = pss;

	cfs = fid_fs - ic->halfw;
	css = fid_ss - ic->halfw;

	bx = add_box(ic);
	bx->refl = refl;
	bx->cfs = cfs;
	bx->css = css;
	bx->p = &ic->image->detgeom->panels[pn];
	bx->pn = pn;

	if ( ic->meth & INTEGRATION_CENTER ) {
		r = center_and_check_box(ic, bx, &saturated);
	} else {
		r = check_box(ic, bx, &saturated);
		if ( !r ) {
			fit_bg(ic, bx);
			if ( !bg_ok(ic, bx) ) r = 1;
		}
		bx->offs_fs = 0.0;
		bx->offs_ss = 0.0;
	}
	if ( r ) {
		delete_box(ic, bx);
		return 1;
	}

	if ( saturated ) {
		ic->n_saturated++;
		if ( !(ic->meth & INTEGRATION_SATURATED) ) {
			delete_box(ic, bx);
			return 1;
		}
	}

	intensity = tentative_intensity(ic, bx);
	mean_var_area(ic, bx, BM_BG, &bgmean, &sig2_bg);

	bias = bx->p->adu_bias;
	intensity -= bias;

	aduph = bx->p->adu_per_photon;
	sig2_poisson = aduph * intensity;

	/* If intensity is within one photon of nothing, set the Poisson
	 * error to be one photon.  Pretend I = 1 photon = 1*aduph,
	 * then sig2_posson = aduph*intensity = aduph^2. */
	if ( fabs(intensity / aduph) < 1.0 ) {
		sig2_poisson = aduph * aduph;
	}

	/* If intensity is negative by more than one photon, assume that
	 * the peak is in the background and add a Poisson error of
	 * appropriate size */
	if ( intensity < -aduph ) {
		sig2_poisson = -aduph*intensity;
	} else if ( intensity < 0.0 ) {
		/* If the intensity is negative (by less than one
		 * photon), assume the reflection is very weak and
		 * therefore has a Poisson error of one photon. */
		sig2_poisson = aduph;
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
	set_detector_pos(refl, pfs, pss);

	if ( get_int_diag(ic, refl) ) show_peak_box(ic, bx, term_lock);

	if ( intensity < -5.0*sigma ) {
		ic->n_implausible++;
		set_redundancy(refl, 0);
		return 1;
	}

	return 0;
}


static double estimate_resolution(Crystal *cr, struct image *image)
{
	int i;
	const double min_dist = 0.25;
	double max_res = 0.0;
	double *acc;
	int n_acc = 0;
	int max_acc = 1024;
	int n;
	double dx, dy;
	UnitCell *cell;


	acc = malloc(max_acc*sizeof(double));
	if ( acc == NULL ) {
		ERROR("Allocation failed during estimate_resolution!\n");
		return INFINITY;
	}

	cell = crystal_get_cell(cr);
	crystal_get_det_shift(cr, &dx, &dy);

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double h, k, l, hd, kd, ld;
		double ax, ay, az;
		double bx, by, bz;
		double cx, cy, cz;
		double r[3];

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		cell_get_cartesian(cell,
		                   &ax, &ay, &az,
		                   &bx, &by, &bz,
		                   &cx, &cy, &cz);

		detgeom_transform_coords(&image->detgeom->panels[f->pn],
		                         f->fs, f->ss, image->lambda,
		                         dx, dy, r);

		/* Decimal and fractional Miller indices of nearest
		 * reciprocal lattice point */
		hd = r[0] * ax + r[1] * ay + r[2] * az;
		kd = r[0] * bx + r[1] * by + r[2] * bz;
		ld = r[0] * cx + r[1] * cy + r[2] * cz;
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
				acc = srealloc(acc, max_acc*sizeof(double));
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
		free(acc);
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
                            pthread_mutex_t *term_lock, int **masks)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	UnitCell *cell;
	struct intcontext *ic;
	int n_rej = 0;
	int n_refl = 0;

	list = crystal_get_reflections(cr);
	cell = crystal_get_cell(cr);

	ic = intcontext_new(image, cell, meth,
	                    ir_inn, ir_mid, ir_out,
	                    masks);
	if ( ic == NULL ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}

	intcontext_set_diag(ic, int_diag, idh, idk, idl);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		n_refl++;
		n_rej += integrate_rings_once(refl, ic,
		                              term_lock);
	}

	intcontext_free(ic);

	if ( n_rej*4 > n_refl ) {
		ERROR("WARNING: %i reflections could not be integrated\n",
		      n_rej);
	}

	crystal_set_num_saturated_reflections(cr, ic->n_saturated);
	crystal_set_num_implausible_reflections(cr, ic->n_implausible);
}


void integrate_all_5(struct image *image, IntegrationMethod meth,
                     PartialityModel pmodel, double push_res,
                     double ir_inn, double ir_mid, double ir_out,
                     IntDiag int_diag,
                     signed int idh, signed int idk, signed int idl,
                     pthread_mutex_t *term_lock, int overpredict)
{
	int i;
	int *masks[image->detgeom->n_panels];

	/* Predict all reflections */
	for ( i=0; i<image->n_crystals; i++ ) {

		RefList *list;
		double res;
		double saved_R = crystal_get_profile_radius(image->crystals[i]);

		if ( overpredict ) {
			crystal_set_profile_radius(image->crystals[i],
			                           saved_R * 5);
		}

		res = estimate_resolution(image->crystals[i], image);
		crystal_set_resolution_limit(image->crystals[i], res);

		list = predict_to_res(image->crystals[i], res+push_res);
		crystal_set_reflections(image->crystals[i], list);

		if ( overpredict ) {
			crystal_set_profile_radius(image->crystals[i], saved_R);
		}

	}

	for ( i=0; i<image->detgeom->n_panels; i++ ) {
		masks[i] = make_BgMask(image, &image->detgeom->panels[i],
		                       i, ir_inn);
	}

	for ( i=0; i<image->n_crystals; i++ ) {

		Crystal *cr = image->crystals[i];

		switch ( meth & INTEGRATION_METHOD_MASK ) {

			case INTEGRATION_NONE :
			break;

			case INTEGRATION_RINGS :
			integrate_rings(meth, cr, image,
			                int_diag, idh, idk, idl,
			                ir_inn, ir_mid, ir_out,
			                term_lock, masks);
			break;

			case INTEGRATION_PROF2D :
			integrate_prof2d(meth, cr, image,
			                 int_diag, idh, idk, idl,
			                 ir_inn, ir_mid, ir_out,
			                 term_lock, masks);
			break;

			default :
			ERROR("Unrecognised integration method %i\n", meth);
			break;

		}

	}

	for ( i=0; i<image->detgeom->n_panels; i++ ) {
		free(masks[i]);
	}
}


void integrate_all_4(struct image *image, IntegrationMethod meth,
                     PartialityModel pmodel, double push_res,
                     double ir_inn, double ir_mid, double ir_out,
                     IntDiag int_diag,
                     signed int idh, signed int idk, signed int idl,
                     pthread_mutex_t *term_lock)
{
	integrate_all_5(image, meth, pmodel, 0.0, ir_inn, ir_mid, ir_out,
	                int_diag, idh, idk, idl, 0, 0);
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
	integrate_all_3(image, meth, PMODEL_XSPHERE, 0.0,
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


char *str_integration_method(IntegrationMethod m)
{
	char tmp[64];
	switch ( m & INTEGRATION_METHOD_MASK ) {
		case INTEGRATION_NONE :
		strcpy(tmp, "none");
		break;

		case INTEGRATION_RINGS :
		strcpy(tmp, "rings");
		break;

		case INTEGRATION_PROF2D :
		strcpy(tmp, "prof2d");
		break;

		default :
		strcpy(tmp, "unknown");
		break;
	}

	if ( m & INTEGRATION_SATURATED ) {
		strcat(tmp, "-sat");
	}

	if ( m & INTEGRATION_CENTER ) {
		strcat(tmp, "-cen");
	}

	if ( m & INTEGRATION_GRADIENTBG ) {
		strcat(tmp, "-grad");
	}

	return strdup(tmp);
}


IntegrationMethod integration_method(const char *str, int *err)
{
	int n, i;
	char **methods;
	IntegrationMethod meth = INTEGRATION_NONE;

	if ( err != NULL ) *err = 0;
	n = assplode(str, ",-", &methods, ASSPLODE_NONE);

	for ( i=0; i<n; i++ ) {

		if ( strcmp(methods[i], "rings") == 0 ) {
			meth = INTEGRATION_DEFAULTS_RINGS;

		} else if ( strcmp(methods[i], "prof2d") == 0 ) {
			meth = INTEGRATION_DEFAULTS_PROF2D;

		} else if ( strcmp(methods[i], "none") == 0 ) {
			return INTEGRATION_NONE;

		} else if ( strcmp(methods[i], "sat") == 0 ) {
			meth |= INTEGRATION_SATURATED;

		} else if ( strcmp(methods[i], "nosat") == 0 ) {
			meth &= ~INTEGRATION_SATURATED;

		} else if ( strcmp(methods[i], "cen") == 0 ) {
			meth |= INTEGRATION_CENTER;

		} else if ( strcmp(methods[i], "nocen") == 0 ) {
			meth &= ~INTEGRATION_CENTER;

		} else if ( strcmp(methods[i], "rescut") == 0 ) {
			ERROR("'rescut'/'norescut' in integration method is no "
			      "longer used.  Set --push-res instead.\n");
			if ( err != NULL ) *err = 1;
			return INTEGRATION_NONE;

		} else if ( strcmp(methods[i], "norescut") == 0 ) {
			ERROR("'rescut'/'norescut' in integration method is no "
			      "longer used.  Set --push-res instead.\n");
			if ( err != NULL ) *err = 1;
			return INTEGRATION_NONE;

		} else if ( strcmp(methods[i], "grad") == 0 ) {
			meth |= INTEGRATION_GRADIENTBG;

		} else if ( strcmp(methods[i], "nograd") == 0 ) {
			meth &= ~INTEGRATION_GRADIENTBG;

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
