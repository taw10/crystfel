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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_fit.h>

#include "merge.h"
#include "post-refinement.h"
#include "symmetry.h"
#include "cell.h"
#include "cell-utils.h"
#include "scaling.h"


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
		double fx, dc;

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

		fx = -log(G) + log(p) - log(L) - B*s*s + log(I_full);
		dc = log(I_partial) - fx;
		w = 1.0;
		dev += w*dc*dc;

		if ( fh != NULL ) {
			fprintf(fh, "%4i %4i %4i %e %e\n",
			        h, k, l, s, dev);
		}

	}

	if ( fh != NULL ) fclose(fh);

	if ( pn_used != NULL ) *pn_used = n_used;
	return dev;
}


struct scale_args
{
	RefList *full;
	Crystal *crystal;
	PartialityModel pmodel;
	int n_reflections;
};


struct queue_args
{
	int n_started;
	int n_done;
	Crystal **crystals;
	int n_crystals;
	long long int n_reflections;
	struct scale_args task_defaults;
};


static void scale_crystal(void *task, int id)
{
	struct scale_args *pargs = task;
	int r;
	double G;

	/* Simple iterative algorithm */
	r = linear_scale(pargs->full, crystal_get_reflections(pargs->crystal), &G, 1);
	if ( r == 0 ) {
		crystal_set_osf(pargs->crystal, G);
	} /* else don't change it */
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
	struct scale_args *wargs = task;

	qa->n_done++;
	qa->n_reflections += wargs->n_reflections;

	progress_bar(qa->n_done, qa->n_crystals, "Scaling");
	free(task);
}


/* Calculates G, by which list2 should be multiplied to fit list1 */
int linear_scale(const RefList *list1, const RefList *list2, double *G,
                 int complain_loudly)
{
	const Reflection *refl1;
	const Reflection *refl2;
	RefListIterator *iter;
	int max_n = 256;
	int n = 0;
	double *x;
	double *y;
	double *w;
	int r;
	double cov11;
	double sumsq;
	int n_esd1 = 0;
	int n_esd2 = 0;
	int n_ih1 = 0;
	int n_ih2 = 0;
	int n_nan1 = 0;
	int n_nan2 = 0;
	int n_inf1 = 0;
	int n_inf2 = 0;
	int n_part = 0;
	int n_nom = 0;

	x = malloc(max_n*sizeof(double));
	w = malloc(max_n*sizeof(double));
	y = malloc(max_n*sizeof(double));
	if ( (x==NULL) || (y==NULL) || (w==NULL) ) {
		ERROR("Failed to allocate memory for scaling.\n");
		return 1;
	}

	int nb = 0;
	for ( refl2 = first_refl_const(list2, &iter);
	      refl2 != NULL;
	      refl2 = next_refl_const(refl2, iter) )
	{
		signed int h, k, l;
		double Ih1, Ih2;
		double esd1, esd2;
		nb++;

		get_indices(refl2, &h, &k, &l);
		refl1 = find_refl(list1, h, k, l);
		if ( refl1 == NULL ) {
			n_nom++;
			continue;
		}

		Ih1 = get_intensity(refl1);
		Ih2 = get_intensity(refl2);
		esd1 = get_esd_intensity(refl1);
		esd2 = get_esd_intensity(refl2);

		/* Problem cases in approximate descending order of severity */
		if ( isnan(Ih1) ) { n_nan1++; continue; }
		if ( isinf(Ih1) ) { n_inf1++; continue; }
		if ( isnan(Ih2) ) { n_nan2++; continue; }
		if ( isinf(Ih2) ) { n_inf2++; continue; }
		if ( get_partiality(refl2) < 0.1 ) { n_part++; continue; }
		if ( Ih1 <= 0.0 ) { n_ih1++; continue; }
		if ( Ih2 <= 0.0 ) { n_ih2++; continue; }
		if ( Ih1 <= 3.0*esd1 ) { n_esd1++;  continue; }
		if ( Ih2 <= 3.0*esd2 ) { n_esd2++; continue; }
		if ( get_redundancy(refl1) < 2 ) continue;

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

		x[n] = Ih2 / get_partiality(refl2);
		y[n] = Ih1;
		w[n] = get_partiality(refl2);
		n++;

	}

	if ( n < 2 ) {
		if ( complain_loudly ) {
			ERROR("Not enough reflections for scaling (had %i, but %i remain)\n", nb, n);
			if ( n_esd1 ) ERROR("%i reference reflection esd\n", n_esd1);
			if ( n_esd2 ) ERROR("%i subject reflection esd\n", n_esd2);
			if ( n_ih1 ) ERROR("%i reference reflection intensity\n", n_ih1);
			if ( n_ih2 ) ERROR("%i subject reflection intensity\n", n_ih2);
			if ( n_nan1 ) ERROR("%i reference reflection nan\n", n_nan1);
			if ( n_nan2 ) ERROR("%i subject reflection nan\n", n_nan2);
			if ( n_inf1 ) ERROR("%i reference reflection inf\n", n_inf1);
			if ( n_inf2 ) ERROR("%i subject reflection inf\n", n_inf2);
			if ( n_part ) ERROR("%i subject reflection partiality\n", n_part);
			if ( n_nom ) ERROR("%i no match in reference list\n", n_nom);
		}
		*G = 1.0;
		return 1;
	}

	r = gsl_fit_wmul(x, 1, w, 1, y, 1, n, G, &cov11, &sumsq);

	if ( r ) {
		ERROR("Scaling failed.\n");
		return 1;
	}

	if ( isnan(*G) ) {

		if ( complain_loudly ) {
			ERROR("Scaling gave NaN (%i pairs)\n", n);
			if ( n < 10 ) {
				int i;
				for ( i=0; i<n; i++ ) {
					STATUS("%i %e %e %e\n", i, x[i], y[i], w[n]);
					}
			}
		}

		*G = 1.0;
		return 1;
	}

	free(x);
	free(y);
	free(w);

	return 0;
}


void scale_all_to_reference(Crystal **crystals, int n_crystals,
                            RefList *reference, int nthreads)
{
	struct scale_args task_defaults;
	struct queue_args qargs;

	task_defaults.crystal = NULL;

	qargs.task_defaults = task_defaults;
	qargs.n_crystals = n_crystals;
	qargs.crystals = crystals;

	/* Don't have threads which are doing nothing */
	if ( n_crystals < nthreads ) nthreads = n_crystals;

	qargs.task_defaults.full = reference;
	qargs.n_started = 0;
	qargs.n_done = 0;
	qargs.n_reflections = 0;
	run_threads(nthreads, scale_crystal, get_crystal, done_crystal,
	            &qargs, n_crystals, 0, 0, 0);
}
