/*
 * post-refinement.c
 *
 * Post refinement
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
#include <gsl/gsl_multimin.h>

#include "image.h"
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "cell-utils.h"
#include "reflist-utils.h"


struct prdata
{
	int refined;
};

const char *str_prflag(enum prflag flag)
{
	switch ( flag ) {

		case PRFLAG_OK :
		return "OK";

		case PRFLAG_FEWREFL :
		return "not enough reflections";

		case PRFLAG_SOLVEFAIL :
		return "PR solve failed";

		case PRFLAG_EARLY :
		return "early rejection";

		case PRFLAG_CC :
		return "low CC";

		case PRFLAG_BIGB :
		return "B too big";

		default :
		return "Unknown flag";
	}
}


double residual(Crystal *cr, const RefList *full, int free,
                int *pn_used, const char *filename)
{
	double dev = 0.0;
	double G, B;
	Reflection *refl;
	RefListIterator *iter;
	FILE *fh = NULL;
	int n_used = 0;

	if ( filename != NULL ) {
		fh = fopen(filename, "w");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", filename);
		}
	}

	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);

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
		I_full = get_intensity(match);

		if ( get_redundancy(match) < 2 ) continue;

		p = get_partiality(refl);
		L = get_lorentz(refl);
		I_partial = get_intensity(refl);
		esd = get_esd_intensity(refl);
		s = resolution(crystal_get_cell(cr), h, k, l);

		if ( I_partial < 3.0*esd ) continue;

		fx = exp(G)*p*exp(-B*s*s)*I_full/L;
		dc = I_partial - fx;
		w = (s/1e9)*(s/1e9)/(esd*esd);
		dev += w*dc*dc;
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


static UnitCell *rotate_cell_xy(const UnitCell *cell, double ang1, double ang2)
{
	UnitCell *o;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xnew, ynew, znew;

	o = cell_new_from_cell(cell);

	cell_get_reciprocal(o, &asx, &asy, &asz,
	                       &bsx, &bsy, &bsz,
	                       &csx, &csy, &csz);

	/* "a" around x */
	xnew = asx;
	ynew = asy*cos(ang1) + asz*sin(ang1);
	znew = -asy*sin(ang1) + asz*cos(ang1);
	asx = xnew;  asy = ynew;  asz = znew;

	/* "b" around x */
	xnew = bsx;
	ynew = bsy*cos(ang1) + bsz*sin(ang1);
	znew = -bsy*sin(ang1) + bsz*cos(ang1);
	bsx = xnew;  bsy = ynew;  bsz = znew;

	/* "c" around x */
	xnew = csx;
	ynew = csy*cos(ang1) + csz*sin(ang1);
	znew = -csy*sin(ang1) + csz*cos(ang1);
	csx = xnew;  csy = ynew;  csz = znew;

	/* "a" around y */
	xnew = asx*cos(ang1) + asz*sin(ang2);
	ynew = asy;
	znew = -asx*sin(ang2) + asz*cos(ang2);
	asx = xnew;  asy = ynew;  asz = znew;

	/* "b" around y */
	xnew = bsx*cos(ang1) + bsz*sin(ang2);
	ynew = bsy;
	znew = -bsx*sin(ang2) + bsz*cos(ang2);
	bsx = xnew;  bsy = ynew;  bsz = znew;

	/* "c" around y */
	xnew = csx*cos(ang1) + csz*sin(ang2);
	ynew = csy;
	znew = -csx*sin(ang2) + csz*cos(ang2);
	csx = xnew;  csy = ynew;  csz = znew;

	cell_set_reciprocal(o, asx, asy, asy, bsx, bsy, bsz, csx, csy, csz);

	return o;
}


struct rf_priv {
	const Crystal *cr;
	const RefList *full;
	enum gparam *rv;
	int verbose;
};


static double residual_f(const gsl_vector *v, void *pp)
{
	struct rf_priv *pv = pp;
	int i;
	UnitCell *cell;
	RefList *list;
	Crystal *cr;
	double res;
	double ang1 = 0.0;
	double ang2 = 0.0;

	for ( i=0; i<v->size; i++ ) {
		switch ( pv->rv[i] ) {

			case GPARAM_ANG1 :
			ang1 = gsl_vector_get(v, i);
			break;

			case GPARAM_ANG2 :
			ang2 = gsl_vector_get(v, i);
			break;

			default :
			ERROR("Don't understand parameter %i\n", pv->rv[i]);
			break;

		}
	}

	cell = rotate_cell_xy(crystal_get_cell_const(pv->cr), ang1, ang2);
	cr = crystal_copy(pv->cr);
	list = copy_reflist(crystal_get_reflections(cr));
	crystal_set_reflections(cr, list);
	crystal_set_cell(cr, cell);

	update_predictions(cr);
	calculate_partialities(cr, PMODEL_XSPHERE);

	res = residual(cr, pv->full, 0, NULL, NULL);

	if ( isnan(res) ) {
		ERROR("NaN residual\n");
		ERROR("G=%e, B=%e\n", crystal_get_osf(cr), crystal_get_Bfac(cr));
		residual(cr, pv->full, 0, NULL, "nan-residual.dat");
		abort();
	}

	cell_free(cell);
	reflist_free(list);
	crystal_free(cr);

	return res;
}


static double get_initial_param(Crystal *cr, enum gparam p)
{
	switch ( p ) {

		case GPARAM_ANG1 : return 0.0;
		case GPARAM_ANG2 : return 0.0;

		default: return 0.0;

	}
}


static double get_stepsize(enum gparam p)
{
	switch ( p ) {

		case GPARAM_ANG1 : return deg2rad(0.01);
		case GPARAM_ANG2 : return deg2rad(0.01);

		default : return 0.0;

	}
}


static void do_pr_refine(Crystal *cr, const RefList *full,
                         PartialityModel pmodel, int verbose)
{
	int i;
	gsl_multimin_fminimizer *min;
	gsl_vector *v;
	gsl_vector *step;
	gsl_multimin_function f;
	enum gparam rv[32];
	struct rf_priv residual_f_priv;
	int n_params = 0;
	int n_iter = 0;
	int status;

	if ( verbose ) {
		STATUS("\nPR initial: dev = %10.5e, free dev = %10.5e\n",
		       residual(cr, full, 0, NULL, NULL),
		       residual(cr, full, 1, NULL, NULL));
	}

	/* The parameters to be refined */
	rv[n_params++] = GPARAM_ANG1;
	rv[n_params++] = GPARAM_ANG2;

	residual_f_priv.cr = cr;
	residual_f_priv.full = full;
	residual_f_priv.rv = rv;
	residual_f_priv.verbose = 1;
	f.f = residual_f;
	f.n = n_params;
	f.params = &residual_f_priv;

	v = gsl_vector_alloc(n_params);
	step = gsl_vector_alloc(n_params);

	for ( i=0; i<n_params; i++ ) {
		gsl_vector_set(v, i, get_initial_param(cr, rv[i]));
		gsl_vector_set(step, i, get_stepsize(rv[i]));
	}

	min = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,
	                                    n_params);
	gsl_multimin_fminimizer_set(min, &f, v, step);

	do {
		n_iter++;

		status = gsl_multimin_fminimizer_iterate(min);
		if ( status ) break;

		if ( verbose ) {
			double res = residual_f(min->x, &residual_f_priv);
			STATUS("%f %f residual = %e\n",
			       rad2deg(gsl_vector_get(min->x, 0)),
			       rad2deg(gsl_vector_get(min->x, 1)), res);
		}

		status = gsl_multimin_test_size(min->size, 1.0e-8);

	} while ( status == GSL_CONTINUE && n_iter < 100 );

	if ( verbose ) {
		STATUS("Done with refinement after %i iter\n", n_iter);
		STATUS("status = %i (%s)\n", status, gsl_strerror(status));
	}

	double tang;
	tang = fabs(gsl_vector_get(min->x, 0)) + fabs(gsl_vector_get(min->x, 1));
	if ( rad2deg(tang) > 5.0 ) {
		ERROR("More than 5 degrees total rotation!\n");
		residual_f_priv.verbose = 1;
		double res = residual_f(min->x, &residual_f_priv);
		STATUS("residual after rotation = %e\n", res);
		residual_f_priv.verbose = 2;
		res = residual_f(v, &residual_f_priv);
		STATUS("residual before rotation = %e\n", res);
		return;
	}

	if ( verbose ) {
		STATUS("PR final: dev = %10.5e, free dev = %10.5e\n",
		       residual(cr, full, 0, NULL, NULL),
		       residual(cr, full, 1, NULL, NULL));
	}

	gsl_multimin_fminimizer_free(min);
	gsl_vector_free(v);
	gsl_vector_free(step);
}


static struct prdata pr_refine(Crystal *cr, const RefList *full,
                               PartialityModel pmodel)
{
	int verbose = 0;
	struct prdata prdata;

	prdata.refined = 0;

	do_pr_refine(cr, full, pmodel, verbose);

	if ( crystal_get_user_flag(cr) == 0 ) {
		prdata.refined = 1;
	}

	return prdata;
}


struct refine_args
{
	RefList *full;
	Crystal *crystal;
	PartialityModel pmodel;
	struct prdata prdata;
};


struct queue_args
{
	int n_started;
	int n_done;
	Crystal **crystals;
	int n_crystals;
	struct refine_args task_defaults;
};


static void refine_image(void *task, int id)
{
	struct refine_args *pargs = task;
	Crystal *cr = pargs->crystal;

	pargs->prdata = pr_refine(cr, pargs->full, pargs->pmodel);
}


static void *get_image(void *vqargs)
{
	struct refine_args *task;
	struct queue_args *qargs = vqargs;

	task = malloc(sizeof(struct refine_args));
	memcpy(task, &qargs->task_defaults, sizeof(struct refine_args));

	task->crystal = qargs->crystals[qargs->n_started];

	qargs->n_started++;

	return task;
}


static void done_image(void *vqargs, void *task)
{
	struct queue_args *qa = vqargs;

	qa->n_done++;

	progress_bar(qa->n_done, qa->n_crystals, "Refining");
	free(task);
}


void refine_all(Crystal **crystals, int n_crystals,
                RefList *full, int nthreads, PartialityModel pmodel)
{
	struct refine_args task_defaults;
	struct queue_args qargs;

	task_defaults.full = full;
	task_defaults.crystal = NULL;
	task_defaults.pmodel = pmodel;
	task_defaults.prdata.refined = 0;

	qargs.task_defaults = task_defaults;
	qargs.n_started = 0;
	qargs.n_done = 0;
	qargs.n_crystals = n_crystals;
	qargs.crystals = crystals;

	/* Don't have threads which are doing nothing */
	if ( n_crystals < nthreads ) nthreads = n_crystals;

	run_threads(nthreads, refine_image, get_image, done_image,
	            &qargs, n_crystals, 0, 0, 0);
}
