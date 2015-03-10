/*
 * hrs-scaling.c
 *
 * Intensity scaling using generalised HRS target function
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


/* Minimum partiality of a reflection for it to be used for scaling */
#define MIN_PART_SCALE (0.05)

/* Minimum partiality of a reflection for it to be merged */
#define MIN_PART_MERGE (0.05)

/* Maximum number of iterations of scaling per macrocycle. */
#define MAX_CYCLES (100)


struct scale_queue_args
{
	RefList *reference;
	Crystal **crystals;
	int n_started;
	PartialityModel pmodel;
	UnitCell *cell;
};


struct scale_worker_args
{
	Crystal *crystal;
	RefList *reference;
	PartialityModel pmodel;
	int crystal_number;
};


static void *create_scale_job(void *vqargs)
{
	struct scale_worker_args *wargs;
	struct scale_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct scale_worker_args));
	wargs->reference = qargs->reference;
	wargs->pmodel = qargs->pmodel;

	wargs->crystal_number = qargs->n_started;
	wargs->crystal = qargs->crystals[qargs->n_started++];

	return wargs;
}


static void run_scale_job(void *vwargs, int cookie)
{
	struct scale_worker_args *wargs = vwargs;
	Crystal *cr = wargs->crystal;
	RefList *reference = wargs->reference;
	Reflection *refl;
	RefListIterator *iter;
	int n = 0;
	double *x;
	double *y;
	double *w;
	int max_n = 256;
	double c0, c1, cov00, cov01, cov11, chisq;
	double G, B;
	int r;

	/* If this crystal's scaling was dodgy, it doesn't contribute to the
	 * merged intensities */
	if ( crystal_get_user_flag(cr) != 0 ) return;

	x = malloc(max_n*sizeof(double));
	y = malloc(max_n*sizeof(double));
	w = malloc(max_n*sizeof(double));
	if ( (x==NULL) || (y==NULL) || (w==NULL) ) {
		ERROR("Failed to allocate memory for scaling.\n");
		return;
	}

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double Ih, Ihl, corr;
		Reflection *r;
		double res;

		get_indices(refl, &h, &k, &l);
		res = resolution(crystal_get_cell(cr), h, k, l);

		if ( (get_partiality(refl) < MIN_PART_SCALE)
		  || (get_intensity(refl) < 5.0*get_esd_intensity(refl)) ) {
			continue;
		}

		/* Look up by asymmetric indices */
		r = find_refl(reference, h, k, l);
		if ( r == NULL ) {
			continue;
		}

		Ih = get_intensity(r);
		corr = get_partiality(refl) * get_lorentz(refl);
		Ihl = get_intensity(refl) / corr;

		if ( Ihl <= 0.0 ) continue;
		if ( Ih <= 0.0 ) continue;
		if ( isnan(Ih) || isinf(Ihl) ) continue;
		if ( isnan(Ih) || isinf(Ih) ) continue;

		if ( n == max_n ) {
			max_n *= 2;
			x = realloc(x, max_n*sizeof(double));
			y = realloc(y, max_n*sizeof(double));
			w = realloc(w, max_n*sizeof(double));
			if ( (x==NULL) || (y==NULL) || (w==NULL) ) {
				ERROR("Failed to allocate memory for scaling.\n");
				return;
			}
		}

		x[n] = res*res;
		y[n] = log(Ihl/Ih);
		w[n] = 1.0;//log(get_esd_intensity(refl)*corr);
		n++;

	}

	if ( n < 2 ) {
		crystal_set_user_flag(cr, 1);
		return;
	}

	r = gsl_fit_wlinear(x, 1, w, 1, y, 1, n, &c0, &c1,
	                     &cov00, &cov01, &cov11, &chisq);

	G = 1.0/exp(c0);
	B = -c1/2.0;

	free(x);
	free(y);
	free(w);
	if ( r || isnan(c0) ) {
		crystal_set_user_flag(cr, 1);
		return;
	}

	crystal_set_osf(cr, G);
	crystal_set_Bfac(cr, B);
}


static void finalise_scale_job(void *vqargs, void *vwargs)
{
	struct scale_worker_args *wargs = vwargs;
	free(wargs);
}


static void iterate_scale(Crystal **crystals, int n, RefList *reference,
                          int n_threads, PartialityModel pmodel)
{
	struct scale_queue_args qargs;

	assert(reference != NULL);

	qargs.reference = reference;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.pmodel = pmodel;

	run_threads(n_threads, run_scale_job, create_scale_job,
	            finalise_scale_job, &qargs, n, 0, 0, 0);
}


struct merge_queue_args
{
	RefList *full;
	pthread_rwlock_t full_lock;
	Crystal **crystals;
	int n_started;
	PartialityModel pmodel;
};


struct merge_worker_args
{
	Crystal *crystal;
	RefList *full;
	pthread_rwlock_t *full_lock;
	PartialityModel pmodel;
	int crystal_number;
};


static void *create_merge_job(void *vqargs)
{
	struct merge_worker_args *wargs;
	struct merge_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct merge_worker_args));
	wargs->full = qargs->full;
	wargs->full_lock = &qargs->full_lock;
	wargs->pmodel = qargs->pmodel;

	wargs->crystal_number = qargs->n_started;
	wargs->crystal = qargs->crystals[qargs->n_started++];

	return wargs;
}


static void run_merge_job(void *vwargs, int cookie)
{
	struct merge_worker_args *wargs = vwargs;
	Crystal *cr = wargs->crystal;
	RefList *full = wargs->full;
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
		pthread_rwlock_rdlock(wargs->full_lock);
		f = find_refl(full, h, k, l);
		if ( f == NULL ) {

			/* Swap read lock for write lock */
			pthread_rwlock_unlock(wargs->full_lock);
			pthread_rwlock_wrlock(wargs->full_lock);

			/* In the gap between the unlock and the wrlock, the
			 * reflection might have been created by another thread.
			 * So, we must check again */
			f = find_refl(full, h, k, l);
			if ( f == NULL ) {
				f = add_refl(full, h, k, l);
				lock_reflection(f);
				pthread_rwlock_unlock(wargs->full_lock);
				set_intensity(f, 0.0);
				set_temp1(f, 0.0);
				set_temp2(f, 0.0);

			} else {

				/* Someone else created it */
				lock_reflection(f);
				pthread_rwlock_unlock(wargs->full_lock);

			}

		} else {

			lock_reflection(f);
			pthread_rwlock_unlock(wargs->full_lock);

		}

		mean = get_intensity(f);
		sumweight = get_temp1(f);
		M2 = get_temp2(f);

		res = resolution(crystal_get_cell(cr), h, k, l);

		/* Total (multiplicative) correction factor */
		corr = G * exp(2.0*B*res*res) * get_lorentz(refl)
		        / get_partiality(refl);

		esd = get_esd_intensity(refl) * corr;
		w = 1.0 / (esd*esd);

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
                         PartialityModel pmodel)
{
	RefList *full;
	struct merge_queue_args qargs;
	Reflection *refl;
	RefListIterator *iter;

	full = reflist_new();

	qargs.full = full;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.pmodel = pmodel;
	pthread_rwlock_init(&qargs.full_lock, NULL);

	run_threads(n_threads, run_merge_job, create_merge_job,
	            finalise_merge_job, &qargs, n, 0, 0, 0);

	pthread_rwlock_destroy(&qargs.full_lock);

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));
	}

	return full;
}


static void reject_outliers(double *old_osfs, int n, Crystal **crystals)
{
	int i;

	for ( i=0; i<n; i++ ) {

		double osf = crystal_get_osf(crystals[i]);
		double Bfac = crystal_get_Bfac(crystals[i]);
		int bad = 0;

		if ( isnan(osf) || (osf < 0.0) || (osf > 10.0) ) bad = 1;
		if ( isnan(Bfac) || (Bfac<-40e-20) || (Bfac>40e-20) ) bad = 1;

		if ( bad ) {
			crystal_set_user_flag(crystals[i], 1);
			//STATUS("crystal %i bad: osf=%f B=%f A^2\n",
			//       i, osf, Bfac*1e20);
		} else {
			//STATUS("OK %i: osf=%f B=%f A^2\n",
			//       i, osf, Bfac*1e20);
		}

	}
}


static void reset_scaling_flag(Crystal *crystal)
{
	/* Every pattern gets another chance at being scaled, but stays flagged
	 * if it was flagged for any other reason */
	if ( crystal_get_user_flag(crystal) == 1 ) {
		crystal_set_user_flag(crystal, 0);
	}
}


static int test_convergence(double *old_osfs, int n, Crystal **crystals)
{
	int i;
	double total_change = 0.0;
	double mean_change;
	int n_change = 0;

	for ( i=0; i<n; i++ ) {
		if ( crystal_get_user_flag(crystals[i]) == 0 ) {
			double new_osf = crystal_get_osf(crystals[i]);
			total_change += fabs(new_osf - old_osfs[i]);
			n_change++;
		}
	}
	mean_change = total_change / n_change;

	STATUS("Mean OSF change = %f\n", mean_change);

	return mean_change < 0.01;
}


/* Scale the stack of images */
RefList *scale_intensities(Crystal **crystals, int n, int n_threads,
                           PartialityModel pmodel, int min_redundancy)
{
	int i;
	RefList *full = NULL;
	double *old_osfs;
	int done;

	for ( i=0; i<n; i++ ) {
		reset_scaling_flag(crystals[i]);
		crystal_set_osf(crystals[i], 1.0);
		crystal_set_Bfac(crystals[i], 0.0);
	}

	/* Create an initial list to refine against */
	full = lsq_intensities(crystals, n, n_threads, pmodel);

	old_osfs = malloc(n*sizeof(double));
	if ( old_osfs == NULL ) return NULL;

	/* Iterate */
	i = 0;
	do {

		double total_sf = 0.0;
		int n_sf = 0;
		double norm_sf;
		int j;

		for ( j=0; j<n; j++ ) {
			old_osfs[j] = crystal_get_osf(crystals[j]);
			reset_scaling_flag(crystals[j]);
		}

		iterate_scale(crystals, n, full, n_threads, pmodel);
		reject_outliers(old_osfs, n, crystals);

		/* Normalise the scale factors */
		for ( j=0; j<n; j++ ) {
			double osf =  crystal_get_osf(crystals[j]);
			if ( crystal_get_user_flag(crystals[j]) == 0 ) {
				total_sf += osf;
				n_sf++;
			}
		}
		norm_sf = total_sf / n_sf;
		for ( j=0; j<n; j++ ) {
			crystal_set_osf(crystals[j],
			                crystal_get_osf(crystals[j])/norm_sf);
		}

		done = test_convergence(old_osfs, n, crystals);

		/* Generate list for next iteration */
		reflist_free(full);
		full = lsq_intensities(crystals, n, n_threads, pmodel);

		i++;

	} while ( !done && (i < MAX_CYCLES) );

	if ( i == MAX_CYCLES ) {
		ERROR("WARNING: Scaling did not converge.\n");
	}

	free(old_osfs);
	return full;
}
