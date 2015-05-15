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
#define MIN_PART_MERGE (0.05)


struct merge_queue_args
{
	RefList *full;
	pthread_rwlock_t full_lock;
	Crystal **crystals;
	int n_started;
	PartialityModel pmodel;
	double push_res;
};


struct merge_worker_args
{
	struct merge_queue_args *qargs;
	Crystal *crystal;
	int crystal_number;
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


static void run_merge_job(void *vwargs, int cookie)
{
	struct merge_worker_args *wargs = vwargs;
	Crystal *cr = wargs->crystal;
	RefList *full = wargs->qargs->full;
	double push_res = wargs->qargs->push_res;
	Reflection *refl;
	RefListIterator *iter;
	double G, B;

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
		double corr, res, w, esd;

		if ( get_partiality(refl) < MIN_PART_MERGE ) continue;

		get_indices(refl, &h, &k, &l);
		pthread_rwlock_rdlock(&wargs->qargs->full_lock);
		f = find_refl(full, h, k, l);
		if ( f == NULL ) {

			/* Swap read lock for write lock */
			pthread_rwlock_unlock(&wargs->qargs->full_lock);
			pthread_rwlock_wrlock(&wargs->qargs->full_lock);

			/* In the gap between the unlock and the wrlock, the
			 * reflection might have been created by another thread.
			 * So, we must check again */
			f = find_refl(full, h, k, l);
			if ( f == NULL ) {
				f = add_refl(full, h, k, l);
				lock_reflection(f);
				pthread_rwlock_unlock(&wargs->qargs->full_lock);
				set_intensity(f, 0.0);
				set_temp1(f, 0.0);
				set_temp2(f, 0.0);

			} else {

				/* Someone else created it */
				lock_reflection(f);
				pthread_rwlock_unlock(&wargs->qargs->full_lock);

			}

		} else {

			lock_reflection(f);
			pthread_rwlock_unlock(&wargs->qargs->full_lock);

		}

		mean = get_intensity(f);
		sumweight = get_temp1(f);
		M2 = get_temp2(f);

		res = resolution(crystal_get_cell(cr), h, k, l);

		if ( 2.0*res > crystal_get_resolution_limit(cr)+push_res ) {
			unlock_reflection(f);
			continue;
		}

		/* Total (multiplicative) correction factor */
		corr = exp(-G) * exp(B*res*res) * get_lorentz(refl)
		        / get_partiality(refl);

		esd = get_esd_intensity(refl) * corr;
		w = 1.0;

		/* Running mean and variance calculation */
		temp = w + sumweight;
		delta = get_intensity(refl)*corr - mean;
		R = delta * w / temp;
		set_intensity(f, mean + R);
		set_temp2(f, M2 + sumweight * delta * R);
		set_temp1(f, temp);
		set_redundancy(f, get_redundancy(f)+1);
		unlock_reflection(f);
	}
}


static void finalise_merge_job(void *vqargs, void *vwargs)
{
	free(vwargs);
}


RefList *lsq_intensities(Crystal **crystals, int n, int n_threads,
                         PartialityModel pmodel, int min_meas, double push_res)
{
	RefList *full;
	RefList *full2;
	struct merge_queue_args qargs;
	Reflection *refl;
	RefListIterator *iter;

	full = reflist_new();

	qargs.full = full;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.pmodel = pmodel;
	qargs.push_res = push_res;
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
