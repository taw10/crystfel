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


struct intcontext
{
	int halfw;
	int w;
	struct image *image;
	gsl_matrix *bgm;
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


static double boxi(struct intcontext *ic,
                   int fid_fs, int fid_ss, int p, int q)
{
	int fs, ss;

	fs = fid_fs + p;
	ss = fid_ss + q;

	/* FIXME: Wrong panel? */
	return ic->image->data[fs + ic->image->width*ss];
}


static void fit_bg(struct intcontext *ic,
                   int fid_fs, int fid_ss,
                   double *pa, double *pb, double *pc)
{
	int p, q;
	gsl_vector *v;
	gsl_vector *ans;

	v = gsl_vector_calloc(3);

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {
		double bi;
		bi = boxi(ic, fid_fs, fid_ss, p, q);

		addv(v, 0, bi*p);
		addv(v, 1, bi*q);
		addv(v, 2, bi);
	}
	}

	ans = solve_svd(v, ic->bgm);
	gsl_vector_free(v);

	*pa = gsl_vector_get(ans, 0);
	*pb = gsl_vector_get(ans, 1);
	*pc = gsl_vector_get(ans, 2);

	gsl_vector_free(ans);
}


static void init_intcontext(struct intcontext *ic)
{
	int p, q;

	ic->w = 2*ic->halfw + 1;

	ic->bgm = gsl_matrix_calloc(3, 3);
	if ( ic->bgm == NULL ) {
		ERROR("Failed to initialise matrix.\n");
		return;
	}

	for ( p=0; p<ic->w; p++ ) {
	for ( q=0; q<ic->w; q++ ) {
		addm(ic->bgm, 0, 0, p*p);
		addm(ic->bgm, 0, 1, p*q);
		addm(ic->bgm, 0, 2, p);
		addm(ic->bgm, 1, 0, p*q);
		addm(ic->bgm, 1, 1, q*q);
		addm(ic->bgm, 1, 2, q);
		addm(ic->bgm, 2, 0, p);
		addm(ic->bgm, 2, 1, q);
		addm(ic->bgm, 2, 2, 1);
	}
	}
}


static void measure_all_intensities(RefList *list, struct image *image)
{
	Reflection *refl;
	RefListIterator *iter;
	struct intcontext ic;

	ic.halfw = 4;
	ic.image = image;
	init_intcontext(&ic);

	for ( refl = first_refl(list, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
	{
		double pfs, pss;
		int fid_fs, fid_ss;
		double a, b, c;

		get_detector_pos(refl, &pfs, &pss);
		fid_fs = lrint(pfs);
		fid_ss = lrint(pss);
		fit_bg(&ic, fid_fs, fid_ss, &a, &b, &c);

		//set_intensity(refl, intensity);
	}
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
