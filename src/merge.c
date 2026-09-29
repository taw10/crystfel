/*
 * merge.c
 *
 * Parallel weighted merging of intensities
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
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
#include <gsl/gsl_statistics.h>

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "utils.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "cell-utils.h"
#include "merge.h"
#include "error_model.h"


struct merge_queue_args
{
	RefList *full;
	pthread_rwlock_t full_lock;
	struct crystal_refls *crystals;
	int n_started;
	double push_res;
	int use_weak;
	long long int n_reflections;
	int ln_merge;
	int n_used;
	ErrorModel *emodel;
};


struct merge_worker_args
{
	struct merge_queue_args *qargs;
	Crystal *crystal;
	RefList *refls;
	int crystal_number;
	int n_reflections;
	int used;
	ErrorModel *emodel;
};


static void *create_merge_job(void *vqargs)
{
	struct merge_worker_args *wargs;
	struct merge_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct merge_worker_args));
	wargs->qargs = qargs;
	wargs->crystal_number = qargs->n_started;
	wargs->crystal = qargs->crystals[qargs->n_started].cr;
	wargs->refls = qargs->crystals[qargs->n_started].refls;
	wargs->used = 0;
	wargs->emodel = qargs->emodel;

	qargs->n_started++;

	return wargs;
}


static int alloc_contribs(struct reflection_contributions *c)
{
	c->contribs = realloc(c->contribs, c->max_contrib*sizeof(double));
	c->contrib_esds = realloc(c->contrib_esds, c->max_contrib*sizeof(double));
	c->contrib_legacy_weights = realloc(c->contrib_legacy_weights, c->max_contrib*sizeof(double));
	if ( c->contribs == NULL ) return 1;
	if ( c->contrib_esds == NULL ) return 1;
	if ( c->contrib_legacy_weights == NULL ) return 1;
	return 0;
}


/* Find reflection hkl in 'list', creating it if it's not there, under
 * protection of 'lock' and returning a locked reflection */
static Reflection *get_locked_reflection(RefList *list, pthread_rwlock_t *lock,
                                         signed int h, signed int k, signed  int l)
{
	Reflection *f;

	pthread_rwlock_rdlock(lock);
	f = find_refl(list, h, k, l);
	if ( f == NULL ) {

		/* Swap read lock for write lock */
		pthread_rwlock_unlock(lock);
		pthread_rwlock_wrlock(lock);

		/* In the gap between the unlock and the wrlock, the
		 * reflection might have been created by another thread.
		 * So, we must check again */
		f = find_refl(list, h, k, l);
		if ( f == NULL ) {

			struct reflection_contributions *c;

			f = add_refl(list, h, k, l);
			lock_reflection(f);
			pthread_rwlock_unlock(lock);
			set_intensity(f, 0.0);
			set_temp1(f, 0.0);
			set_temp2(f, 0.0);

			c = malloc(sizeof(struct reflection_contributions));
			if ( c != NULL ) {
				c->n_contrib = 0;
				c->max_contrib = 32;
				c->contribs = NULL;
				c->contrib_esds = NULL;
				c->contrib_legacy_weights = NULL;
				if ( alloc_contribs(c) ) {
					set_contributions(f, NULL);
				} else {
					set_contributions(f, c);
				}
			} else {
				set_contributions(f, NULL);
			}

		} else {
			/* Someone else created it */
			lock_reflection(f);
			pthread_rwlock_unlock(lock);
		}

	} else {

		lock_reflection(f);
		pthread_rwlock_unlock(lock);

	}

	return f;
}


static void run_merge_job(void *vwargs, int cookie)
{
	struct merge_worker_args *wargs = vwargs;
	Crystal *cr = wargs->crystal;
	RefList *full = wargs->qargs->full;
	double push_res = wargs->qargs->push_res;
	int ln_merge = wargs->qargs->ln_merge;
	Reflection *refl;
	RefListIterator *iter;
	double G, B;

	wargs->n_reflections = 0;

	/* If this crystal's scaling was dodgy, it doesn't contribute to the
	 * merged intensities */
	if ( crystal_get_user_flag(cr) != 0 ) return;

	wargs->used = 1;

	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);

	for ( refl = first_refl(wargs->refls, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		double res;
		double Ii, esdIi, legacy_weight;
		struct reflection_contributions *c;

		if ( get_partiality(refl) < MIN_PART_MERGE ) continue;
		if ( isnan(get_esd_intensity(refl)) ) continue;

		if ( !wargs->qargs->use_weak || ln_merge ) {

			if (get_intensity(refl) < 3.0*fabs(get_esd_intensity(refl))) {
				continue;
			}

			if ( get_flag(refl) ) continue;

		}

		get_indices(refl, &h, &k, &l);
		f = get_locked_reflection(full, &wargs->qargs->full_lock,
		                          h, k, l);

		res = resolution(crystal_get_cell(cr), h, k, l);
		if ( 2.0*res > crystal_get_resolution_limit(cr)+push_res ) {
			unlock_reflection(f);
			continue;
		}

		Ii = correct_reflection(get_intensity(refl), refl, G,  B, res);
		esdIi = correct_reflection(get_esd_intensity(refl), refl, G, B, res);
		legacy_weight = get_partiality(refl) / correct_reflection_nopart(1.0, refl, G, B, res);

		if ( ln_merge ) {
			Ii = log(Ii);
		}

		set_redundancy(f, get_redundancy(f)+1);
		c = get_contributions(f);
		if ( c == NULL ) {
			ERROR("Failed to get contributions for %i %i %i\n", h, k, l);
		} else {
			c->contribs[c->n_contrib] = Ii;
			c->contrib_esds[c->n_contrib] = esdIi;
			c->contrib_legacy_weights[c->n_contrib] = legacy_weight;
			c->n_contrib++;
			if ( c->n_contrib == c->max_contrib ) {
				c->max_contrib += 64;
				alloc_contribs(c);
			}
		}

		unlock_reflection(f);

		wargs->n_reflections++;

	}
}


static void finalise_merge_job(void *vqargs, void *vwargs)
{
	struct merge_queue_args *qargs = vqargs;
	struct merge_worker_args *wargs = vwargs;
	qargs->n_reflections += wargs->n_reflections;
	qargs->n_used += wargs->used;
	free(vwargs);
}


RefList *merge_intensities(struct crystal_refls *crystals, int n,
                           int n_threads, int min_meas,
                           double push_res, int use_weak, int ln_merge,
                           ErrorModel *emodel, int refine_emodel,
                           int *pn_used)
{
	RefList *full;
	RefList *out;
	struct merge_queue_args qargs;
	Reflection *refl_in;
	RefListIterator *iter;

	if ( n == 0 ) return NULL;

	full = reflist_new();

	qargs.full = full;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.push_res = push_res;
	qargs.use_weak = use_weak;
	qargs.n_reflections = 0;
	qargs.ln_merge = ln_merge;
	qargs.n_used = 0;
	qargs.emodel = emodel;
	pthread_rwlock_init(&qargs.full_lock, NULL);

	run_threads(n_threads, run_merge_job, create_merge_job,
	            finalise_merge_job, &qargs, n, 0, 0, 0);

	pthread_rwlock_destroy(&qargs.full_lock);

	out = reflist_new();
	for ( refl_in = first_refl(full, &iter);
	      refl_in != NULL;
	      refl_in = next_refl(refl_in, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		struct reflection_contributions *c;
		double *weights;
		double wmean;

		c = get_contributions(refl_in);

		if ( get_redundancy(refl_in) < min_meas ) {

			/* We do not need the contribution list any more */
			struct reflection_contributions *c;
			c = get_contributions(refl_in);
			free(c->contribs);
			free(c->contrib_esds);
			free(c->contrib_legacy_weights);
			free(c);

			continue;

		}

		get_indices(refl_in, &h, &k, &l);
		f = add_refl(out, h, k, l);

		set_contributions(f, c);
		set_redundancy(f, get_redundancy(refl_in));

		/* Unweighted mean (used for Ev06, Ev11 etc, also deltaCChalf) */
		set_unweighted_mean(f, gsl_stats_mean(c->contribs, 1, c->n_contrib));

		/* Max scaled intensity (used for Kh23) */
		set_max_measurement(f, gsl_stats_max(c->contribs, 1, c->n_contrib));

		weights = make_weights_array(c, f, emodel);

		wmean = gsl_stats_wmean(weights, 1, c->contribs, 1, c->n_contrib);
		set_intensity(f, wmean);

		set_esd_intensity(f, merged_esd(c, weights, wmean, emodel));

		/* Correct for averaging log of intensities */
		if ( ln_merge ) {
			set_unweighted_mean(f, exp(get_unweighted_mean(f)));
			set_intensity(f, exp(get_intensity(f)));
			set_esd_intensity(f, exp(get_esd_intensity(f)));
			set_max_measurement(f, exp(get_max_measurement(f)));
		}

		free(weights);

	}

	reflist_free(full);

	if ( refine_emodel ) {
		refine_error_model(out, emodel);
	}

	if ( pn_used != NULL ) {
		*pn_used = qargs.n_used;
	}

	return out;
}


/* Correct 'val' (probably an intensity from one pattern, maybe an e.s.d.)
 * for scaling and Lorentz factors but not partiality nor polarisation */
double correct_reflection_nopart(double val, Reflection *refl, double osf,
                                 double Bfac, double res)
{
	double corr = osf * exp(-Bfac*res*res);
	return (val / corr) / get_lorentz(refl);
}


/* Correct 'val' (probably an intensity from one pattern, maybe an e.s.d.)
 * for scaling, partiality and Lorentz factors but not polarisation */
double correct_reflection(double val, Reflection *refl, double osf, double Bfac,
                          double res)
{
	double Ipart = correct_reflection_nopart(val, refl, osf, Bfac, res);
	return Ipart / get_partiality(refl);
}


double residual(RefList *list, Crystal *cr, const RefList *full, int free,
                int *pn_used, const char *filename)
{
	Reflection *refl;
	RefListIterator *iter;
	int n_used = 0;
	double num = 0.0;
	double den = 0.0;
	double G = crystal_get_osf(cr);
	double B = crystal_get_Bfac(cr);
	UnitCell *cell = crystal_get_cell(cr);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double w, res;
		signed int h, k, l;
		Reflection *match;
		double I_full;
		double int1, pobs, pcalc;

		if ( free != get_flag(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);
		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;
		I_full = get_intensity(match);

		if ( get_redundancy(match) < 2 ) continue;

		int1 = correct_reflection_nopart(get_intensity(refl), refl, G, B, res);
		pobs = int1 / I_full;
		if ( pobs > 1.0 ) pobs = 1.0;
		if ( pobs < 0.0 ) pobs = 0.0;

		pcalc = get_partiality(refl);

		w = 1.0 / correct_reflection_nopart(1.0, refl, G, B, res);
		if ( isnan(w) ) {
			w = 0.0;
		}

		num += w*fabs(pobs-pcalc);
		den += w;

		n_used++;

	}

	if ( pn_used != NULL ) *pn_used = n_used;
	return num/den;
}


double log_residual(RefList *list, Crystal *cr,
                    const RefList *full, int free,
                    int *pn_used, const char *filename)
{
	double dev = 0.0;
	double G, B;
	Reflection *refl;
	RefListIterator *iter;
	int n_used = 0;
	FILE *fh = NULL;

	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);
	if ( filename != NULL ) {
		fh = fopen(filename, "a");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", filename);
		}
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double p, L, s, w;
		signed int h, k, l;
		Reflection *match;
		double esd, I_full, I_partial;
		double fx;

		if ( free && !get_flag(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;

		p = get_partiality(refl);
		L = get_lorentz(refl);
		I_partial = get_intensity(refl);
		I_full = get_intensity(match);
		esd = get_esd_intensity(refl);
		s = resolution(crystal_get_cell(cr), h, k, l);

		if ( I_partial <= 3.0*esd ) continue; /* Also because of log */
		if ( get_redundancy(match) < 2 ) continue;
		if ( I_full <= 0 ) continue;  /* Because log */
		if ( p <= 0.0 ) continue; /* Because of log */

		fx = log(G) - B*s*s + log(p) + log(I_full) - log(I_partial) - log(L);
		w = 1.0;
		dev += w*fx*fx;
		n_used++;

		if ( fh != NULL ) {
			fprintf(fh, "%4i %4i %4i %e %e\n",
			        h, k, l, s, dev);
		}

	}

	if ( fh != NULL ) fclose(fh);

	if ( pn_used != NULL ) *pn_used = n_used;
	return dev;
}


/* Has to match run_merge_job to be useful */
void write_unmerged(const char *fn,
                    struct crystal_refls *crystals,
                    struct image **images,
                    int n_crystals)
{
	FILE *fh;
	int i;

	fh = fopen(fn, "w");
	if ( fh == NULL ) {
		ERROR("Failed to open %s\n", fn);
		return;
	}

	for ( i=0; i<n_crystals; i++ ) {

		Reflection *refl;
		RefListIterator *iter;
		double G, B;
		UnitCell *cell;

		fprintf(fh, "Crystal %i\n", i);
		fprintf(fh, "Filename: %s %s\n", images[i]->filename, images[i]->ev);
		if ( crystal_get_user_flag(crystals[i].cr) ) {
			fprintf(fh, "Flagged: yes\n");
		} else {
			fprintf(fh, "Flagged: no\n");
		}

		G = crystal_get_osf(crystals[i].cr);
		B = crystal_get_Bfac(crystals[i].cr);
		cell = crystal_get_cell(crystals[i].cr);

		for ( refl = first_refl(crystals[i].refls, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			double res;

			get_indices(refl, &h, &k, &l);
			res = resolution(cell, h, k, l);
			fprintf(fh, "%4i %4i %4i %10.2f %10.2f", h, k, l,
			        correct_reflection(get_intensity(refl), refl, G,  B, res),
			        get_partiality(refl) / correct_reflection_nopart(1.0, refl, G, B, res));

			if ( get_partiality(refl) < MIN_PART_MERGE ) {
				fprintf(fh, " partiality_too_small");
			}
			if ( isnan(get_esd_intensity(refl)) ) {
				fprintf(fh, " nan_esd");
			}
			fprintf(fh, "\n");

		}

		fprintf(fh, "\n");

	}

	fclose(fh);
}


static UnitCell *first_accepted_cell(struct crystal_refls *crystals,
                                     int n_crystals)
{
	int i;
	for ( i=0; i<n_crystals; i++ ) {
		if ( crystal_get_user_flag(crystals[i].cr) ) continue;
		return crystal_get_cell(crystals[i].cr);
	}
	return NULL;
}


void average_unit_cell(struct crystal_refls *crystals,
                       int n_crystals,
                       const char *outcell_filename)
{
	int i;
	UnitCell *cmean;
	double a_sumw = 0.0, a_mean = 0.0, a_M2 = 0.0;
	double b_sumw = 0.0, b_mean = 0.0, b_M2 = 0.0;
	double c_sumw = 0.0, c_mean = 0.0, c_M2 = 0.0;
	double al_sumw = 0.0, al_mean = 0.0, al_M2 = 0.0;
	double be_sumw = 0.0, be_mean = 0.0, be_M2 = 0.0;
	double ga_sumw = 0.0, ga_mean = 0.0, ga_M2 = 0.0;
	char ua, cen;
	LatticeType lt;
	int n_nonmatch = 0;
	int n_acc = 0;

	UnitCell *cell = first_accepted_cell(crystals, n_crystals);
	if ( cell == NULL ) {
		ERROR("No accepted crystals - can't calculate average cell.\n");
		return;
	}
	lt = cell_get_lattice_type(cell);
	ua = cell_get_unique_axis(cell);
	cen = cell_get_centering(cell);

	for ( i=0; i<n_crystals; i++ ) {
		UnitCell *cell;
		double a, b, c, al, be, ga;
		if ( crystal_get_user_flag(crystals[i].cr) ) continue;
		cell = crystal_get_cell(crystals[i].cr);
		if ( (cell_get_lattice_type(cell) != lt)
		  || (cell_get_unique_axis(cell) != ua)
		  || (cell_get_centering(cell) != cen) )
		{
			n_nonmatch++;
			continue;
		}
		cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
		mean_variance(a, 1.0, &a_sumw, &a_mean, &a_M2);
		mean_variance(b, 1.0, &b_sumw, &b_mean, &b_M2);
		mean_variance(c, 1.0, &c_sumw, &c_mean, &c_M2);
		mean_variance(al, 1.0, &al_sumw, &al_mean, &al_M2);
		mean_variance(be, 1.0, &be_sumw, &be_mean, &be_M2);
		mean_variance(ga, 1.0, &ga_sumw, &ga_mean, &ga_M2);
		n_acc++;
	}

	cmean = cell_new_from_parameters(a_mean, b_mean, c_mean,
	                                 al_mean, be_mean, ga_mean);
	cell_set_unique_axis(cmean, ua);
	cell_set_lattice_type(cmean, lt);
	cell_set_centering(cmean, cen);

	STATUS("Average unit cell from %i accepted crystals:\n", n_acc);
	cell_print(cmean);

	if ( n_nonmatch > 0 ) {
		ERROR("WARNING: %i crystals were not included in the average "
		      "cell, because they had different lattice types, "
		      "centering or unique axes.\n", n_nonmatch);
	}

	if ( outcell_filename != NULL ) {
		write_cell_to_file(cmean, outcell_filename);
	}
	cell_free(cmean);
}
