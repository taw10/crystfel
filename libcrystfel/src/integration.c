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

#include "reflist.h"
#include "cell.h"
#include "crystal.h"
#include "cell-utils.h"
#include "geometry.h"
#include "image.h"
#include "peaks.h"
#include "integration.h"


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

	n = num_reflections(list);
	*np = 0;  /* For now */

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
		assert(i <= n);
	}

	qsort(il, n, sizeof(struct integr_ind), compare_resolution);

	*np = n;
	return il;
}


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
	BM_IG,
	BM_BG,
	BM_PK
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
};


struct peak_box
{
	int fid_fs;   /* Coordinates of corner */
	int fid_ss;

	int pn;           /* Panel number */
	struct panel *p;  /* The panel itself */

	/* Fitted background parameters */
	double a;
	double b;
	double c;

	/* Measured intensity (tentative, profile fitted or otherwise) */
	double intensity;

	int rp;   /* Reference profile number */

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

	fs = bx->fid_fs + p;
	ss = bx->fid_ss + q;

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


static void show_peak_box(struct intcontext *ic, struct peak_box *bx)
{
	int q;

	printf("Pixel values                                              ");
	printf("Box flags              ");
	printf("Fitted background\n");

	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {
			printf("%5.0f ", boxi(ic, bx, p, q));
		}

		printf("    ");

		for ( p=0; p<ic->w; p++ ) {
			printf("%i ", ic->bm[p+q*ic->w]);
		}

		printf("    ");

		for ( p=0; p<ic->w; p++ ) {
			printf("%5.0f ", bx->a*p + bx->b*q + bx->c);
		}

		printf("\n");
	}
	printf("-----------> p\n");
	printf("Reference profile number %i\n", bx->rp);

}


static void show_reference_profile(struct intcontext *ic, int i)
{
	int q;

	printf("Reference profile number %i:\n", i);
	printf("Pixel values\n");

	for ( q=ic->w-1; q>=0; q-- ) {

		int p;

		for ( p=0; p<ic->w; p++ ) {
			printf("%5.0f ", ic->reference_profiles[i][p+ic->w*q]);
		}

		printf("\n");
	}
	printf("-----------> p\n");

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

		if ( ic->bm[p + ic->w*q] == BM_BG ) {

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
		}
		}
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


static struct peak_box *add_box(struct intcontext *ic)
{
	int idx;

	if ( ic->n_boxes == ic->max_boxes ) {
		if ( alloc_boxes(ic, ic->max_boxes+32) ) {
			return NULL;
		}
	}

	idx = ic->n_boxes++;

	ic->boxes[idx].verbose = 0;

	return &ic->boxes[idx];
}


static double tentative_intensity(struct intcontext *ic, struct peak_box *bx)
{
	int p, q;
	double intensity = 0.0;

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {

		if ( ic->bm[p + ic->w*q] != BM_PK ) continue;
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
		if ( ic->bm[p + ic->w*q] != BM_PK ) continue;
		bi = boxi(ic, bx, p, q);

		num_p += bi*p;
		num_q += bi*q;
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

		if ( ic->bm[p + ic->w*q] == BM_IG ) continue;
		bi = boxi(ic, bx, p, q);

		val = bi*bx->intensity;
		val -= p*bx->intensity*bx->a + q*bx->intensity*bx->b
		       + bx->intensity*bx->c;

		ic->reference_profiles[bx->rp][p+ic->w*q] += val;
		ic->reference_den[bx->rp][p+ic->w*q] += pow(bx->intensity, 2.0);
	}
	}
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

		max /= 10000.0;

		for ( p=0; p<ic->w; p++ ) {
		for ( q=0; q<ic->w; q++ ) {

			ic->reference_profiles[i][p+ic->w*q] /= max;

		}
		}

	}

	show_reference_profile(ic, 2);
}


static void measure_all_intensities(RefList *list, struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;
	struct intcontext ic;

	ic.halfw = 4;
	ic.image = image;
	if ( init_intcontext(&ic) ) {
		ERROR("Failed to initialise integration.\n");
		return;
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pfs, pss;
		int pw, ph;
		double pos_p, pos_q;
		signed int h, k, l;
		struct peak_box *bx;

		bx = add_box(&ic);

		get_indices(refl, &h, &k, &l);
		if ( (h==-24) && (k==6) && (l==-12) ) {
			bx->verbose = 1;
		}

		get_detector_pos(refl, &pfs, &pss);
		bx->fid_fs = lrint(pfs);
		bx->fid_ss = lrint(pss);
		bx->pn = find_panel_number(image->det, bx->fid_fs, bx->fid_ss);
		bx->p = &image->det->panels[bx->pn];

		bx->fid_fs -= bx->p->min_fs;
		bx->fid_ss -= bx->p->min_ss;

		bx->fid_fs -= ic.halfw;
		bx->fid_ss -= ic.halfw;

		pw = bx->p->w;
		ph = bx->p->h;
		if ( (bx->fid_fs + ic.w >= pw) || (bx->fid_ss + ic.w >= ph ) ) {
			continue;
		}
		if ( (bx->fid_fs < 0) || (bx->fid_ss < 0 ) ) {
			continue;
		}

		fit_bg(&ic, bx);

		observed_position(&ic, bx, &pos_p, &pos_q);
		pos_p -= ic.halfw;
		pos_q -= ic.halfw;
		if ( bx->verbose ) {
			STATUS("%f %f\n", pos_p, pos_q);
		}

		bx->intensity = tentative_intensity(&ic, bx);
		set_intensity(refl, bx->intensity);

		/* Which reference profile? */
		bx->rp = bx->pn;
		add_to_reference_profile(&ic, bx);
	}

	calculate_reference_profiles(&ic);


}


static void estimate_mosaicity(Crystal *cr, struct image *image)
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
	measure_all_intensities(list, image);

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
		if ( !cutoff ) {
			set_redundancy(refl, 1);
		} else {
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


static void integrate_refine(Crystal *cr, struct image *image, int use_closer,
                             double min_snr,
                             double ir_inn, double ir_mid, double ir_out,
                             int integrate_saturated, int **bgMasks)
{
	RefList *reflections;

	/* Create initial list of reflections with nominal parameters */
	reflections = find_intersections(image, cr);
	measure_all_intensities(reflections, image);

	/* Find resolution limit of pattern using this list */
	estimate_resolution(reflections, cr, image);

	reflist_free(reflections);

	STATUS("Initial resolution estimate = %.2f nm^-1 or %.2f A\n",
	       crystal_get_resolution_limit(cr)/1e9,
	       1e9 / crystal_get_resolution_limit(cr));

	/* Estimate the mosaicity of the crystal using this resolution limit */
	estimate_mosaicity(cr, image);

	/* Create new list of reflections with refined mosaicity */
	reflections = find_intersections(image, cr);
	measure_all_intensities(reflections, image);

	estimate_resolution(reflections, cr, image);
}


static void integrate_rings(Crystal *cr, struct image *image, int use_closer,
                            double min_snr,
                            double ir_inn, double ir_mid, double ir_out,
                            int integrate_saturated, int **bgMasks)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	UnitCell *cell;
	int n_saturated = 0;
	double limit = 0.0;

	list = find_intersections(image, cr);
	if ( list == NULL ) return;

	if ( num_reflections(list) == 0 ) return;

	cell = crystal_get_cell(cr);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double fs, ss, intensity;
		double d;
		int idx;
		double sigma, snr;
		double pfs, pss;
		int r;
		signed int h, k, l;
		struct panel *p;
		int pnum, j, found;
		int saturated;
		double one_over_d;

		get_detector_pos(refl, &pfs, &pss);
		get_indices(refl, &h, &k, &l);

		/* Is there a really close feature which was detected? */
		if ( use_closer ) {

			struct imagefeature *f;

			if ( image->features != NULL ) {
				f = image_feature_closest(image->features,
					                  pfs, pss, &d, &idx);
			} else {
				f = NULL;
			}

			/* FIXME: Horrible hardcoded value */
			if ( (f != NULL) && (d < 10.0) ) {

				double exe;

				exe = get_excitation_error(refl);

				pfs = f->fs;
				pss = f->ss;

				set_detector_pos(refl, exe, pfs, pss);

			}

		}

		p = find_panel(image->det, pfs, pss);
		if ( p == NULL ) continue;  /* Next peak */
		found = 0;
		for ( j=0; j<image->det->n_panels; j++ ) {
			if ( &image->det->panels[j] == p ) {
				pnum = j;
				found = 1;
				break;
			}
		}
		if ( !found ) {
			ERROR("Couldn't find panel %p in list.\n", p);
			return;
		}

		r = integrate_peak(image, pfs, pss, &fs, &ss,
		                   &intensity, &sigma, ir_inn, ir_mid, ir_out,
		                   bgMasks[pnum], &saturated);

		if ( !r && saturated ) {
			n_saturated++;
			if ( !integrate_saturated ) r = 1;
		}

		/* I/sigma(I) cutoff
		 * Rejects reflections below --min-integration-snr, or if the
		 * SNR is clearly silly.  Silly indicates that the intensity
		 * was zero. */
		snr = fabs(intensity)/sigma;
		if ( !r && (isnan(snr) || (snr < min_snr)) ) {
			r = 1;
		}

		/* Record intensity and set redundancy to 1 on success */
		if ( !r ) {
			set_intensity(refl, intensity);
			set_esd_intensity(refl, sigma);
			set_redundancy(refl, 1);
		} else {
			set_redundancy(refl, 0);
		}

		one_over_d = resolution(cell, h, k, l);
		if ( one_over_d > limit ) limit = one_over_d;

	}

	crystal_set_num_saturated_reflections(cr, n_saturated);
	crystal_set_resolution_limit(cr, limit);
	crystal_set_reflections(cr, list);
}


void integrate_all(struct image *image, IntegrationMethod meth,
                   int use_closer, double min_snr,
                   double ir_inn, double ir_mid, double ir_out,
                   int integrate_saturated)
{
	int i;
	int **bgMasks;

	/* Make background masks for all panels */
	bgMasks = calloc(image->det->n_panels, sizeof(int *));
	if ( bgMasks == NULL ) {
		ERROR("Couldn't create list of background masks.\n");
		return;
	}
	for ( i=0; i<image->det->n_panels; i++ ) {
		int *mask;
		mask = make_BgMask(image, &image->det->panels[i], ir_inn);
		if ( mask == NULL ) {
			ERROR("Couldn't create background mask.\n");
			return;
		}
		bgMasks[i] = mask;
	}

	for ( i=0; i<image->n_crystals; i++ ) {

		switch ( meth & INTEGRATION_METHOD_MASK ) {

			case INTEGRATION_NONE :
			return;

			case INTEGRATION_RINGS :
			integrate_rings(image->crystals[i], image, use_closer,
			                min_snr, ir_inn, ir_mid, ir_out,
				        integrate_saturated, bgMasks);
			return;

			case INTEGRATION_REFINE :
			integrate_refine(image->crystals[i], image, use_closer,
			                 min_snr, ir_inn, ir_mid, ir_out,
				         integrate_saturated, bgMasks);
			return;

			default :
			ERROR("Unrecognised integration method %i\n", meth);
			return;

		}

	}

	for ( i=0; i<image->det->n_panels; i++ ) {
		free(bgMasks[i]);
	}
	free(bgMasks);
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

		} else if ( strcmp(methods[i], "refine") == 0) {
			meth = INTEGRATION_DEFAULTS_REFINE;

		} else if ( strcmp(methods[i], "none") == 0) {
			return INTEGRATION_NONE;

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
