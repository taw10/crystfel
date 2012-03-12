/*
 * hrs-scaling.c
 *
 * Intensity scaling using generalised HRS target function
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
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
#include <gsl/gsl_blas.h>

#include "image.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "utils.h"
#include "reflist.h"


/* Maximum number of iterations of scaling per macrocycle. */
#define MAX_CYCLES (50)

/* ESD of restraint driving scale factors to unity */
#define SCALING_RESTRAINT (1.0)


struct scale_queue_args
{
	RefList *reference;
	struct image *images;
	int n_started;
	double max_shift;
};


struct scale_worker_args
{
	struct image *image;
	double shift;
	RefList *reference;
};


static void *create_scale_job(void *vqargs)
{
	struct scale_worker_args *wargs;
	struct scale_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct scale_worker_args));
	wargs->reference = qargs->reference;

	wargs->image = &qargs->images[qargs->n_started++];

	return wargs;
}


static void run_scale_job(void *vwargs, int cookie)
{
	struct scale_worker_args *wargs = vwargs;
	struct image *image = wargs->image;
	RefList *reference = wargs->reference;
	Reflection *refl;
	RefListIterator *iter;
	double num = 0.0;
	double den = 0.0;
	double corr;

	if ( image->pr_dud ) {
		wargs->shift = 0.0;
		return;
	}

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double Ih, Ihl, esd;
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
			if ( get_redundancy(r) < 2 ) continue;
			Ih = get_intensity(r);
		}

		Ihl = get_intensity(refl) / get_partiality(refl);
		esd = get_esd_intensity(refl) / get_partiality(refl);

		num += Ih * (Ihl/image->osf) / pow(esd/image->osf, 2.0);
		den += pow(Ih, 2.0)/pow(esd/image->osf, 2.0);

	}

	//num += image->osf / pow(SCALING_RESTRAINT, 2.0);
	//den += pow(image->osf, 2.0)/pow(SCALING_RESTRAINT, 2.0);

	corr = num / den;
	if ( !isnan(corr) && !isinf(corr) ) {
		image->osf *= corr;
	}
	wargs->shift = fabs(corr-1.0);

}


static void finalise_scale_job(void *vqargs, void *vwargs)
{
	struct scale_queue_args *qargs = vqargs;
	struct scale_worker_args *wargs = vwargs;

	if ( wargs->shift > qargs->max_shift ) qargs->max_shift = wargs->shift;
	free(wargs);
}


static double iterate_scale(struct image *images, int n, RefList *reference,
                            int n_threads)
{
	struct scale_queue_args qargs;

	assert(reference != NULL);

	qargs.reference = reference;
	qargs.n_started = 0;
	qargs.images = images;
	qargs.max_shift = 0.0;

	run_threads(n_threads, run_scale_job, create_scale_job,
	            finalise_scale_job, &qargs, n, 0, 0, 0);

	return qargs.max_shift;
}


struct merge_queue_args
{
	RefList *full;
	pthread_mutex_t full_lock;
	struct image *images;
	int n_started;
};


struct merge_worker_args
{
	struct image *image;
	RefList *full;
	pthread_mutex_t *full_lock;
};


static void *create_merge_job(void *vqargs)
{
	struct merge_worker_args *wargs;
	struct merge_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct merge_worker_args));
	wargs->full = qargs->full;
	wargs->full_lock = &qargs->full_lock;

	wargs->image = &qargs->images[qargs->n_started++];

	return wargs;
}


static void run_merge_job(void *vwargs, int cookie)
{
	struct merge_worker_args *wargs = vwargs;
	struct image *image = wargs->image;
	RefList *full = wargs->full;
	Reflection *refl;
	RefListIterator *iter;
	double G;

	if ( image->pr_dud ) return;

	G = image->osf;

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		double num, den;
		int red;
		double Ihl, esd, pcalc;

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

		pcalc = get_partiality(refl);
		Ihl = get_intensity(refl) / pcalc;
		esd = get_esd_intensity(refl) / pcalc;

		num += (Ihl/G) / pow(esd/G, 2.0);
		den += 1.0 / pow(esd/G, 2.0);
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


static RefList *lsq_intensities(struct image *images, int n, int n_threads)
{
	RefList *full;
	struct merge_queue_args qargs;
	Reflection *refl;
	RefListIterator *iter;

	full = reflist_new();

	qargs.full = full;
	qargs.n_started = 0;
	qargs.images = images;
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
	struct image *images;
	int n_started;
};


struct esd_worker_args
{
	struct image *image;
	RefList *full;
};


static void *create_esd_job(void *vqargs)
{
	struct esd_worker_args *wargs;
	struct esd_queue_args *qargs = vqargs;

	wargs = malloc(sizeof(struct esd_worker_args));
	wargs->full = qargs->full;

	wargs->image = &qargs->images[qargs->n_started++];

	return wargs;
}


static void run_esd_job(void *vwargs, int cookie)
{
	struct esd_worker_args *wargs = vwargs;
	struct image *image = wargs->image;
	RefList *full = wargs->full;
	Reflection *refl;
	RefListIterator *iter;
	double G;

	if ( image->pr_dud ) return;

	G = image->osf;

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		Reflection *f;
		signed int h, k, l;
		double num;
		double Ihl, Ih;

		if ( !get_scalable(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		f = find_refl(full, h, k, l);
		assert(f != NULL);

		lock_reflection(f);

		num = get_temp1(f);

		Ih = get_intensity(f);
		Ihl = get_intensity(refl) / (get_partiality(refl) * G);

		num += pow(Ihl - Ih, 2.0);

		set_temp1(f, num);
		unlock_reflection(f);
	}
}


static void finalise_esd_job(void *vqargs, void *vwargs)
{
	free(vwargs);
}


static void calculate_esds(struct image *images, int n, RefList *full,
                           int n_threads, int min_red)
{
	struct esd_queue_args qargs;
	Reflection *refl;
	RefListIterator *iter;

	qargs.full = full;
	qargs.n_started = 0;
	qargs.images = images;

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


/* Scale the stack of images */
RefList *scale_intensities(struct image *images, int n, RefList *gref,
                           int n_threads, int noscale)
{
	int i;
	double max_corr;
	RefList *full = NULL;
	const int min_redundancy = 3;

	for ( i=0; i<n; i++ ) images[i].osf = 1.0;

	if ( noscale ) {
		full = lsq_intensities(images, n, n_threads);
		calculate_esds(images, n, full, n_threads, min_redundancy);
		return full;
	}

	/* No reference -> create an initial list to refine against */
	if ( gref == NULL ) {
		full = lsq_intensities(images, n, n_threads);
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

		max_corr = iterate_scale(images, n, reference, n_threads);
		//STATUS("Scaling iteration %2i: max correction = %5.2f\n",
		//       i+1, max_corr);

		/* No reference -> generate list for next iteration */
		if ( gref == NULL ) {
			reflist_free(full);
			full = lsq_intensities(images, n, n_threads);
		}

		//show_scale_factors(images, n);

		i++;

	} while ( (max_corr > 0.01) && (i < MAX_CYCLES) );

	if ( gref != NULL ) {
		full = lsq_intensities(images, n, n_threads);
	} /* else we already did it */

	calculate_esds(images, n, full, n_threads, min_redundancy);

	return full;
}
