/*
 * hrs-scaling.c
 *
 * Intensity scaling using generalised HRS target function
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
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

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "utils.h"
#include "reflist.h"


/* Maximum number of iterations of scaling per macrocycle. */
#define MAX_CYCLES (10)


struct scale_queue_args
{
	RefList *reference;
	Crystal **crystals;
	int n_started;
	double max_shift;
	PartialityModel pmodel;
};


struct scale_worker_args
{
	Crystal *crystal;
	double shift;
	RefList *reference;
	PartialityModel pmodel;
};


static void *create_scale_job(void *vqargs)
{
	struct scale_worker_args *wargs;
	struct scale_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct scale_worker_args));
	wargs->reference = qargs->reference;
	wargs->pmodel = qargs->pmodel;

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
	double num = 0.0;
	double den = 0.0;
	double g;
	const double G = crystal_get_osf(cr);

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double Ih, Ihl, esd, corr;
		Reflection *r;

		if ( !get_scalable(refl) ) continue;

		/* Look up by asymmetric indices */
		get_indices(refl, &h, &k, &l);
		r = find_refl(reference, h, k, l);
		if ( r == NULL ) {
			ERROR("%3i %3i %3i isn't in the "
			      "reference list, so why is it "
			      "marked as scalable?\n", h, k, l);
			Ih = 0.0;
		} else {
			Ih = get_intensity(r);
		}

		/* If you change this, be sure to also change
		 * run_merge_job() and run_esd_job(). */
		switch ( wargs->pmodel ) {

			case PMODEL_UNITY :
			corr = 1.0;
			break;

			case PMODEL_SPHERE :
			corr = get_partiality(refl) * get_lorentz(refl);
			break;

			default :
			ERROR("Unrecognised partiality model!\n");
			abort();
			break;

		}

		Ihl = get_intensity(refl) / corr;
		esd = get_esd_intensity(refl) / corr;

		num += Ih * Ihl;
		den += Ihl * Ihl;

	}

	g = num / den;
	if ( !isnan(g) && !isinf(g) ) {
		crystal_set_osf(cr, g);
		wargs->shift = fabs((g/G)-1.0);
	} else {
		wargs->shift = 0.0;
	}
}


static void finalise_scale_job(void *vqargs, void *vwargs)
{
	struct scale_queue_args *qargs = vqargs;
	struct scale_worker_args *wargs = vwargs;

	if ( wargs->shift > qargs->max_shift ) qargs->max_shift = wargs->shift;
	free(wargs);
}


static double iterate_scale(Crystal **crystals, int n, RefList *reference,
                            int n_threads, PartialityModel pmodel)
{
	struct scale_queue_args qargs;

	assert(reference != NULL);

	qargs.reference = reference;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.max_shift = 0.0;
	qargs.pmodel = pmodel;

	run_threads(n_threads, run_scale_job, create_scale_job,
	            finalise_scale_job, &qargs, n, 0, 0, 0);

	return qargs.max_shift;
}


struct merge_queue_args
{
	RefList *full;
	pthread_mutex_t full_lock;
	Crystal **crystals;
	int n_started;
	PartialityModel pmodel;
};


struct merge_worker_args
{
	Crystal *crystal;
	RefList *full;
	pthread_mutex_t *full_lock;
	PartialityModel pmodel;
};


static void *create_merge_job(void *vqargs)
{
	struct merge_worker_args *wargs;
	struct merge_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct merge_worker_args));
	wargs->full = qargs->full;
	wargs->full_lock = &qargs->full_lock;
	wargs->pmodel = qargs->pmodel;

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
	double G;

	G = crystal_get_osf(cr);

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		double num, den;
		int red;
		double Ihl, esd, corr;

		if ( !get_scalable(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		/* FIXME (somehow): Huge contention on this lock */
		pthread_mutex_lock(wargs->full_lock);
		f = find_refl(full, h, k, l);
		if ( f == NULL ) {
			f = add_refl(full, h, k, l);
			lock_reflection(f);
			pthread_mutex_unlock(wargs->full_lock);
			num = 0.0;
			den = 0.0;
			red = 0;
		} else {
			lock_reflection(f);
			pthread_mutex_unlock(wargs->full_lock);
			num = get_temp1(f);
			den = get_temp2(f);
			red = get_redundancy(f);
		}

		/* If you change this, be sure to also change
		 * run_scale_job() and run_esd_job(). */
		switch ( wargs->pmodel ) {

			case PMODEL_UNITY :
			corr = 1.0;
			break;

			case PMODEL_SPHERE :
			corr = get_partiality(refl) * get_lorentz(refl);
			break;

			default :
			ERROR("Unrecognised partiality model!\n");
			abort();
			break;

		}

		Ihl = get_intensity(refl) / corr;
		esd = get_esd_intensity(refl) / corr;

		num += Ihl * G;
		den += 1.0;
		red++;

		set_temp1(f, num);
		set_temp2(f, den);
		set_redundancy(f, red);
		unlock_reflection(f);
	}
}


static void finalise_merge_job(void *vqargs, void *vwargs)
{
	free(vwargs);
}


static RefList *lsq_intensities(Crystal **crystals, int n, int n_threads,
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
	pthread_mutex_init(&qargs.full_lock, NULL);

	run_threads(n_threads, run_merge_job, create_merge_job,
	            finalise_merge_job, &qargs, n, 0, 0, 0);

	pthread_mutex_destroy(&qargs.full_lock);

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double Ih;
		Ih = get_temp1(refl) / get_temp2(refl);
		set_intensity(refl, Ih);
	}

	return full;
}



struct esd_queue_args
{
	RefList *full;
	Crystal **crystals;
	int n_started;
	PartialityModel pmodel;
};


struct esd_worker_args
{
	Crystal *crystal;
	RefList *full;
	PartialityModel pmodel;
};


static void *create_esd_job(void *vqargs)
{
	struct esd_worker_args *wargs;
	struct esd_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct esd_worker_args));
	wargs->full = qargs->full;
	wargs->pmodel = qargs->pmodel;

	wargs->crystal = qargs->crystals[qargs->n_started++];

	return wargs;
}


static void run_esd_job(void *vwargs, int cookie)
{
	struct esd_worker_args *wargs = vwargs;
	Crystal *cr = wargs->crystal;
	RefList *full = wargs->full;
	Reflection *refl;
	RefListIterator *iter;
	double G;

	G = crystal_get_osf(cr);

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		double num;
		double Ihl, Ih, corr;

		if ( !get_scalable(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		f = find_refl(full, h, k, l);
		assert(f != NULL);

		lock_reflection(f);

		num = get_temp1(f);

		/* If you change this, be sure to also change
		 * run_scale_job() and run_merge_job(). */
		switch ( wargs->pmodel ) {

			case PMODEL_UNITY :
			corr = 1.0;
			break;

			case PMODEL_SPHERE :
			corr = get_partiality(refl) * get_lorentz(refl);
			break;;

			default :
			ERROR("Unrecognised partiality model!\n");
			abort();
			break;

		}

		Ih = get_intensity(f);
		Ihl = get_intensity(refl) / (corr * G);

		num += pow(Ihl - Ih, 2.0);

		set_temp1(f, num);
		unlock_reflection(f);
	}
}


static void finalise_esd_job(void *vqargs, void *vwargs)
{
	free(vwargs);
}


static void calculate_esds(Crystal **crystals, int n, RefList *full,
                           int n_threads, int min_red, PartialityModel pmodel)
{
	struct esd_queue_args qargs;
	Reflection *refl;
	RefListIterator *iter;

	qargs.full = full;
	qargs.n_started = 0;
	qargs.crystals = crystals;
	qargs.pmodel = pmodel;

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		set_temp1(refl, 0.0);
		set_temp2(refl, 0.0);
	}

	run_threads(n_threads, run_esd_job, create_esd_job,
	            finalise_esd_job, &qargs, n, 0, 0, 0);

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double esd;
		int red = get_redundancy(refl);
		esd = sqrt(get_temp1(refl));
		esd /= (double)red;
		set_esd_intensity(refl, esd);

		if ( red < min_red ) {
			set_redundancy(refl, 0);
		}
	}
}


static double max_outlier_change(Crystal **crystals, int n, double *old_osfs)
{
	double rdtot = 0.0;
	double rdmax = 0.0;
	double rdmean;
	int j;

	for ( j=0; j<n; j++ ) {
		double r;
		r = crystal_get_osf(crystals[j]) / old_osfs[j];
		rdtot += r;
	}
	rdmean = rdtot / n;
	for ( j=0; j<n; j++ ) {
		double r;
		r = crystal_get_osf(crystals[j]) / old_osfs[j];
		old_osfs[j] = crystal_get_osf(crystals[j]);
		if ( fabs(r-rdmean) > rdmax ) {
			rdmax = fabs(r-rdmean);
		}
	}

	return rdmax;
}


/* Scale the stack of images */
RefList *scale_intensities(Crystal **crystals, int n, RefList *gref,
                           int n_threads, int noscale, PartialityModel pmodel,
                           int min_redundancy)
{
	int i;
	double max_corr;
	RefList *full = NULL;
	double *old_osfs;
	double rdval;

	old_osfs = malloc(n*sizeof(double));
	if ( old_osfs == NULL ) return NULL;

	for ( i=0; i<n; i++ ) crystal_set_osf(crystals[i], 1.0);

	if ( noscale ) {
		full = lsq_intensities(crystals, n, n_threads, pmodel);
		calculate_esds(crystals, n, full, n_threads, min_redundancy,
		               pmodel);
		free(old_osfs);
		return full;
	}

	/* No reference -> create an initial list to refine against */
	if ( gref == NULL ) {
		full = lsq_intensities(crystals, n, n_threads, pmodel);
	}

	for ( i=0; i<n; i++ ) {
		old_osfs[i] = crystal_get_osf(crystals[i]);
	}

	/* Iterate */
	i = 0;
	do {

		RefList *reference;

		/* Refine against reference or current "full" estimates */
		if ( gref != NULL ) {
			reference = gref;
		} else {
			reference = full;
		}

		max_corr = iterate_scale(crystals, n, reference, n_threads,
		                         pmodel);

		rdval = max_outlier_change(crystals, n, old_osfs);

		STATUS("Scaling iteration %2i: max correction = %5.2f "
		       "dev of worst outlier = %5.2f\n",
		       i+1, max_corr, rdval);

		/* No reference -> generate list for next iteration */
		if ( gref == NULL ) {
			reflist_free(full);
			full = lsq_intensities(crystals, n, n_threads, pmodel);
		}

		i++;

	} while ( (rdval > 0.05) && (i < MAX_CYCLES) );

	if ( i == MAX_CYCLES ) {
		ERROR("Warning: Scaling did not converge.\n");
	}

	if ( gref != NULL ) {
		full = lsq_intensities(crystals, n, n_threads, pmodel);
	} /* else we already did it */

	calculate_esds(crystals, n, full, n_threads, min_redundancy, pmodel);

	free(old_osfs);
	return full;
}
