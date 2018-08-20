/*
 * merge.c
 *
 * Parallel weighted merging of intensities
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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fit.h>

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "utils.h"
#include "reflist.h"
#include "cell-utils.h"


/* Minimum partiality of a reflection for it to be merged */
#define MIN_PART_MERGE (0.3)


struct merge_queue_args
{
	RefList *full;
	pthread_rwlock_t full_lock;
	Crystal **crystals;
	int n_started;
	double push_res;
	int use_weak;
	long long int n_reflections;
};


struct merge_worker_args
{
	struct merge_queue_args *qargs;
	Crystal *crystal;
	int crystal_number;
	int n_reflections;
};


static void *create_merge_job(void *vqargs)
{
	struct merge_worker_args *wargs;
	struct merge_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct merge_worker_args));
	wargs->qargs = qargs;
	wargs->crystal_number = qargs->n_started;
	wargs->crystal = qargs->crystals[qargs->n_started++];

	return wargs;
}


static int alloc_contribs(struct reflection_contributions *c)
{
	c->contribs = realloc(c->contribs, c->max_contrib*sizeof(Reflection *));
	c->contrib_crystals = realloc(c->contrib_crystals,
	                              c->max_contrib*sizeof(Crystal *));
	if ( c->contribs == NULL ) return 1;
	if ( c->contrib_crystals == NULL ) return 1;
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
				c->contrib_crystals = NULL;
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
	Reflection *refl;
	RefListIterator *iter;
	double G, B;

	wargs->n_reflections = 0;

	/* If this crystal's scaling was dodgy, it doesn't contribute to the
	 * merged intensities */
	if ( crystal_get_user_flag(cr) != 0 ) return;

	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		double mean, sumweight, M2, temp, delta, R;
		double corr, res, w;
		struct reflection_contributions *c;

		if ( get_partiality(refl) < MIN_PART_MERGE ) continue;

		if ( !wargs->qargs->use_weak ) {

			if (get_intensity(refl) < 3.0*get_esd_intensity(refl)) {
				continue;
			}

			if ( get_flag(refl) ) continue;

		}

		get_indices(refl, &h, &k, &l);
		f = get_locked_reflection(full, &wargs->qargs->full_lock,
		                          h, k, l);

		mean = get_intensity(f);
		sumweight = get_temp1(f);
		M2 = get_temp2(f);

		res = resolution(crystal_get_cell(cr), h, k, l);

		if ( 2.0*res > crystal_get_resolution_limit(cr)+push_res ) {
			unlock_reflection(f);
			continue;
		}

		/* Total (multiplicative) correction factor */
		corr = 1.0/G * exp(B*res*res) * get_lorentz(refl)
		        / get_partiality(refl);
		if ( isnan(corr) ) {
			ERROR("NaN corr:\n");
			ERROR("G = %f, B = %e\n", G, B);
			ERROR("res = %e\n", res);
			ERROR("p = %f\n", get_partiality(refl));
			unlock_reflection(f);
			continue;
		}

		/* Reflections count less the more they have to be scaled up */
		w = 1.0;

		/* Running mean and variance calculation */
		temp = w + sumweight;
		delta = get_intensity(refl)*corr - mean;
		R = delta * w / temp;
		set_intensity(f, mean + R);
		set_temp2(f, M2 + sumweight * delta * R);
		set_temp1(f, temp);
		set_redundancy(f, get_redundancy(f)+1);

		/* Record this contribution */
		c = get_contributions(f);
		if ( c != NULL ) {
			c->contribs[c->n_contrib] = refl;
			c->contrib_crystals[c->n_contrib++] = cr;
			if ( c->n_contrib == c->max_contrib ) {
				c->max_contrib += 64;
				alloc_contribs(c);
			}
		} /* else, too bad! */

		unlock_reflection(f);

		wargs->n_reflections++;

	}
}


static void finalise_merge_job(void *vqargs, void *vwargs)
{
	struct merge_queue_args *qargs = vqargs;
	struct merge_worker_args *wargs = vwargs;
	qargs->n_reflections += wargs->n_reflections;
	free(vwargs);
}


RefList *merge_intensities(Crystal **crystals, int n, int n_threads,
                           int min_meas,
                           double push_res, int use_weak)
{
	RefList *full;
	RefList *full2;
	struct merge_queue_args qargs;
	Reflection *refl;
	RefListIterator *iter;

	if ( n == 0 ) return NULL;

	full = reflist_new();

	qargs.full = full;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.push_res = push_res;
	qargs.use_weak = use_weak;
	qargs.n_reflections = 0;
	pthread_rwlock_init(&qargs.full_lock, NULL);

	run_threads(n_threads, run_merge_job, create_merge_job,
	            finalise_merge_job, &qargs, n, 0, 0, 0);

	pthread_rwlock_destroy(&qargs.full_lock);

	/* Calculate ESDs from variances, including only reflections with
	 * enough measurements */
	full2 = reflist_new();
	if ( full2 == NULL ) return NULL;
	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));

		if ( red >= min_meas ) {

			signed int h, k, l;
			Reflection *r2;

			get_indices(refl, &h, &k, &l);
			r2 =  add_refl(full2, h, k, l);
			copy_data(r2, refl);
		}
	}

	reflist_free(full);
	return full2;
}


double residual(Crystal *cr, const RefList *full, int free,
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

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double p, w, corr, res;
		signed int h, k, l;
		Reflection *match;
		double I_full;
		double int1, int2;

		if ( free && !get_flag(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);
		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;
		I_full = get_intensity(match);

		if ( get_redundancy(match) < 2 ) continue;

		p = get_partiality(refl);
		//if ( p < 0.2 ) continue;

		corr = get_lorentz(refl) / (G * exp(-B*res*res));
		int1 = get_intensity(refl) * corr;
		int2 = p*I_full;
		w = p;

		num += fabs(int1 - int2) * w;
		den += fabs(int1 + int2) * w/2.0;

		n_used++;

	}

	if ( pn_used != NULL ) *pn_used = n_used;
	return num/(den*sqrt(2));
}


double log_residual(Crystal *cr, const RefList *full, int free,
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

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
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

		if ( fh != NULL ) {
			fprintf(fh, "%4i %4i %4i %e %e\n",
			        h, k, l, s, dev);
		}

	}

	if ( fh != NULL ) fclose(fh);

	if ( pn_used != NULL ) *pn_used = n_used;
	return dev;
}
