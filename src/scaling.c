/*
 * scaling.c
 *
 * Scaling
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2017 Thomas White <taw@physics.org>
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
#include <gsl/gsl_statistics_double.h>

#include "merge.h"
#include "post-refinement.h"
#include "symmetry.h"
#include "cell.h"
#include "cell-utils.h"
#include "scaling.h"


struct scale_args
{
	RefList *full;
	Crystal *crystal;
	int flags;
};


struct queue_args
{
	int n_started;
	int n_done;
	Crystal **crystals;
	int n_crystals;
	struct scale_args task_defaults;
};


static void scale_crystal(void *task, int id)
{
	struct scale_args *pargs = task;
	scale_one_crystal(pargs->crystal, pargs->full, pargs->flags);
}


static void *get_crystal(void *vqargs)
{
	struct scale_args *task;
	struct queue_args *qargs = vqargs;

	task = malloc(sizeof(struct scale_args));
	memcpy(task, &qargs->task_defaults, sizeof(struct scale_args));

	task->crystal = qargs->crystals[qargs->n_started];

	qargs->n_started++;

	return task;
}


static void done_crystal(void *vqargs, void *task)
{
	struct queue_args *qa = vqargs;
	qa->n_done++;
	progress_bar(qa->n_done, qa->n_crystals, "Scaling");
	free(task);
}


static double total_log_r(Crystal **crystals, int n_crystals, RefList *full,
                          int *ninc)
{
	int i;
	double total = 0.0;
	int n = 0;

	for ( i=0; i<n_crystals; i++ ) {
		double r;
		if ( crystal_get_user_flag(crystals[i]) ) continue;
		r = log_residual(crystals[i], full, 0, NULL, NULL);
		if ( isnan(r) ) continue;
		total += r;
		n++;
	}

	if ( ninc != NULL ) *ninc = n;
	return total;
}


/* Perform iterative scaling, all the way to convergence */
void scale_all(Crystal **crystals, int n_crystals, int nthreads, int scaleflags)
{
	struct scale_args task_defaults;
	struct queue_args qargs;
	double old_res, new_res;
	int niter = 0;

	task_defaults.crystal = NULL;
	task_defaults.flags = scaleflags;

	qargs.task_defaults = task_defaults;
	qargs.n_crystals = n_crystals;
	qargs.crystals = crystals;

	/* Don't have threads which are doing nothing */
	if ( n_crystals < nthreads ) nthreads = n_crystals;

	new_res = INFINITY;
	do {
		RefList *full;
		int ninc;
		double bef_res;

		full = merge_intensities(crystals, n_crystals, nthreads,
		                         2, INFINITY, 0, 1);
		old_res = new_res;
		bef_res = total_log_r(crystals, n_crystals, full, NULL);

		qargs.task_defaults.full = full;
		qargs.n_started = 0;
		qargs.n_done = 0;
		run_threads(nthreads, scale_crystal, get_crystal, done_crystal,
		            &qargs, n_crystals, 0, 0, 0);

		new_res = total_log_r(crystals, n_crystals, full, &ninc);
		STATUS("Log residual went from %e to %e, %i crystals\n",
		       bef_res, new_res, ninc);

		int i;
		double meanB = 0.0;
		for ( i=0; i<n_crystals; i++ ) {
			meanB += crystal_get_Bfac(crystals[i]);
		}
		meanB /= n_crystals;
		STATUS("Mean B = %e\n", meanB);

		reflist_free(full);
		niter++;

	} while ( (fabs(new_res-old_res) >= 0.01*old_res) && (niter < 10) );

	if ( niter == 10 ) {
		ERROR("Too many iterations - giving up!\n");
	}
}


/* Calculates G and B, by which cr's reflections should be multiplied to fit reference */
int scale_one_crystal(Crystal *cr, const RefList *listR, int flags)
{
	const Reflection *reflS;
	RefListIterator *iter;
	int max_n = 256;
	int n = 0;
	double *x;
	double *y;
	double *w;
	int r;
	double cov00, cov01, cov11, chisq;
	int n_esdS = 0;
	int n_esdR = 0;
	int n_ihS = 0;
	int n_ihR = 0;
	int n_nanS = 0;
	int n_nanR = 0;
	int n_infS = 0;
	int n_infR = 0;
	int n_part = 0;
	int n_nom = 0;
	int n_red = 0;
	RefList *listS = crystal_get_reflections(cr);
	UnitCell *cell = crystal_get_cell(cr);
	double G, B;

	assert(cell != NULL);
	assert(listR != NULL);
	assert(listS != NULL);

	x = malloc(max_n*sizeof(double));
	w = malloc(max_n*sizeof(double));
	y = malloc(max_n*sizeof(double));
	if ( (x==NULL) || (y==NULL) || (w==NULL) ) {
		ERROR("Failed to allocate memory for scaling.\n");
		return 1;
	}

	int nb = 0;
	for ( reflS = first_refl_const(listS, &iter);
	      reflS != NULL;
	      reflS = next_refl_const(reflS, iter) )
	{
		signed int h, k, l;
		const Reflection *reflR;
		double IhR, IhS, esdS, pS, LS;
		double s;

		nb++;

		get_indices(reflS, &h, &k, &l);
		reflR = find_refl(listR, h, k, l);
		if ( reflR == NULL ) {
			n_nom++;
			continue;
		}

		s = resolution(cell, h, k, l);

		IhR = get_intensity(reflR);
		IhS = get_intensity(reflS);
		esdS = get_esd_intensity(reflS);
		pS = get_partiality(reflS);
		LS = get_lorentz(reflS);

		/* Problem cases in approximate descending order of severity */
		if ( isnan(IhR) ) { n_nanR++; continue; }
		if ( isinf(IhR) ) { n_infR++; continue; }
		if ( isnan(IhS) ) { n_nanS++; continue; }
		if ( isinf(IhS) ) { n_infS++; continue; }
		if ( pS < 0.3 ) { n_part++; continue; }
		if ( IhS <= 0.0 ) { n_ihS++; continue; }
		if ( IhS <= 3.0*esdS ) { n_esdS++; continue; }
		if ( IhR <= 0.0 ) { n_ihR++; continue; }
		if ( get_redundancy(reflR) < 2 ) { n_red++; continue; }

		if ( n == max_n ) {
			max_n *= 2;
			x = realloc(x, max_n*sizeof(double));
			y = realloc(y, max_n*sizeof(double));
			w = realloc(w, max_n*sizeof(double));
			if ( (x==NULL) || (y==NULL) || (w==NULL) ) {
				ERROR("Failed to allocate memory for scaling.\n");
				return 1;
			}
		}

		x[n] = s*s;
		y[n] = log(LS) + log(IhS) -log(pS) - log(IhR);
		w[n] = pS;
		n++;

	}

	if ( n < 2 ) {
		if ( flags & SCALE_VERBOSE_ERRORS ) {
			ERROR("Not enough reflections for scaling (had %i, but %i remain)\n", nb, n);
			if ( n_esdR ) ERROR("%i reference reflection esd\n", n_esdR);
			if ( n_esdS ) ERROR("%i subject reflection esd\n", n_esdS);
			if ( n_ihR ) ERROR("%i reference reflection intensity\n", n_ihR);
			if ( n_red ) ERROR("%i reference reflection redundancy\n", n_red);
			if ( n_ihS ) ERROR("%i subject reflection intensity\n", n_ihS);
			if ( n_nanR ) ERROR("%i reference reflection nan\n", n_nanR);
			if ( n_nanS ) ERROR("%i subject reflection nan\n", n_nanS);
			if ( n_infR ) ERROR("%i reference reflection inf\n", n_infR);
			if ( n_infS ) ERROR("%i subject reflection inf\n", n_infS);
			if ( n_part ) ERROR("%i subject reflection partiality\n", n_part);
			if ( n_nom ) ERROR("%i no match in reference list\n", n_nom);
		}
		free(x);
		free(y);
		free(w);
		return 1;
	}

	if ( flags & SCALE_NO_B ) {
		G = gsl_stats_wmean(w, 1, y, 1, n);
		B = 0.0;
		r = 0;
	} else {
		r = gsl_fit_wlinear(x, 1, w, 1, y, 1, n, &G, &B, &cov00, &cov01, &cov11, &chisq);
	}

	if ( r ) {
		ERROR("Scaling failed.\n");
		free(x);
		free(y);
		free(w);
		return 1;
	}

	if ( isnan(G) ) {

		if ( flags & SCALE_VERBOSE_ERRORS ) {
			ERROR("Scaling gave NaN (%i pairs)\n", n);
			if ( n < 10 ) {
				int i;
				for ( i=0; i<n; i++ ) {
					STATUS("%3i %e %e %e\n", i, x[i], y[i], w[i]);
				}
			}
		}

		free(x);
		free(y);
		free(w);
		return 1;
	}

	crystal_set_osf(cr, exp(G));
	crystal_set_Bfac(cr, -B);

	free(x);
	free(y);
	free(w);

	return 0;
}
