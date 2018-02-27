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
#include <gsl/gsl_fit.h>

#include "image.h"
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "cell-utils.h"
#include "reflist-utils.h"
#include "scaling.h"


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
                int *pn_used, const char *filename, int complain)
{
	Reflection *refl;
	RefListIterator *iter;
	int n_used = 0;
	double num = 0.0;
	double den = 0.0;
	double G = crystal_get_osf(cr);
	double B = crystal_get_Bfac(cr);
	UnitCell *cell = crystal_get_cell(cr);

	if ( linear_scale(full, crystal_get_reflections(cr), &G, complain) ) {
		return GSL_NAN;
	}

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

		corr = G * exp(B*res*res) * get_lorentz(refl);
		int1 = get_intensity(refl) * corr;
		int2 = p*I_full;
		w = 1.0;

		num += fabs(int1 - int2) * w;
		den += fabs(int1 + int2) * w/2.0;

		n_used++;

	}

	if ( pn_used != NULL ) *pn_used = n_used;
	return num/(den*sqrt(2));
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
	xnew = asx*cos(ang2) + asz*sin(ang2);
	ynew = asy;
	znew = -asx*sin(ang2) + asz*cos(ang2);
	asx = xnew;  asy = ynew;  asz = znew;

	/* "b" around y */
	xnew = bsx*cos(ang2) + bsz*sin(ang2);
	ynew = bsy;
	znew = -bsx*sin(ang2) + bsz*cos(ang2);
	bsx = xnew;  bsy = ynew;  bsz = znew;

	/* "c" around y */
	xnew = csx*cos(ang2) + csz*sin(ang2);
	ynew = csy;
	znew = -csx*sin(ang2) + csz*cos(ang2);
	csx = xnew;  csy = ynew;  csz = znew;

	cell_set_reciprocal(o, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	return o;
}


/* We set all the step sizes to 1, then scale them.
 * This way, the size of the simplex stays meaningful and we possibly also
 *  avoid some roundoff errors */
static double get_scale(enum gparam p)
{
	switch ( p ) {

		case GPARAM_ANG1 : return deg2rad(0.01);
		case GPARAM_ANG2 : return deg2rad(0.01);
		case GPARAM_R : return 0.0005e9;
		case GPARAM_WAVELENGTH : return 0.001e-10;

		default : return 0.0;

	}
}


struct rf_priv {
	const Crystal *cr;
	const RefList *full;
	enum gparam *rv;
	int verbose;
	int serial;
	const gsl_vector *initial;
};


static double get_actual_val(const gsl_vector *v, const gsl_vector *initial,
                             enum gparam *rv, int i)
{
	return gsl_vector_get(v, i) * get_scale(rv[i])
	             + gsl_vector_get(initial, i);
}


static void apply_parameters(const gsl_vector *v, const gsl_vector *initial,
                             enum gparam *rv, Crystal *cr)
{
	int i;
	double ang1, ang2, R, lambda;
	UnitCell *cell;

	/* Default parameters if not used in refinement */
	ang1 = 0.0;
	ang2 = 0.0;
	R = crystal_get_profile_radius(cr);
	lambda = crystal_get_image(cr)->lambda;

	for ( i=0; i<v->size; i++ ) {

		double val;

		val = get_actual_val(v, initial, rv, i);

		switch ( rv[i] ) {

			case GPARAM_ANG1 :
			ang1 = val;
			break;

			case GPARAM_ANG2 :
			ang2 = val;
			break;

			case GPARAM_R :
			R = val;
			break;

			case GPARAM_WAVELENGTH :
			lambda = val;
			break;

			default :
			ERROR("Don't understand parameter %i\n", rv[i]);
			break;

		}
	}

	cell = rotate_cell_xy(crystal_get_cell_const(cr), ang1, ang2);
	crystal_set_cell(cr, cell);

	crystal_set_profile_radius(cr, R);
	crystal_get_image(cr)->lambda = lambda;
}


static double residual_f(const gsl_vector *v, void *pp)
{
	struct rf_priv *pv = pp;
	RefList *list;
	struct image im;
	Crystal *cr;
	double res;
	int i;

	for ( i=0; i<v->size; i++ ) {
		if ( gsl_vector_get(v, i) > 100.0 ) return GSL_NAN;
	}

	cr = crystal_copy(pv->cr);
	im = *crystal_get_image(cr);
	crystal_set_image(cr, &im);
	apply_parameters(v, pv->initial, pv->rv, cr);

	if ( crystal_get_profile_radius(cr) <= 0.0 ) {
		crystal_free(cr);
		if ( pv->verbose ) STATUS("R < 0\n");
		return GSL_NAN;
	}

	list = copy_reflist(crystal_get_reflections(cr));
	crystal_set_reflections(cr, list);

	update_predictions(cr);
	calculate_partialities(cr, PMODEL_XSPHERE);

	res = residual(cr, pv->full, 0, NULL, NULL, 0);

	cell_free(crystal_get_cell(cr));
	reflist_free(crystal_get_reflections(cr));
	crystal_free(cr);

	return res;
}


static double get_initial_param(Crystal *cr, enum gparam p)
{
	switch ( p ) {

		case GPARAM_ANG1 : return 0.0;
		case GPARAM_ANG2 : return 0.0;
		case GPARAM_R : return crystal_get_profile_radius(cr);
		case GPARAM_WAVELENGTH : return crystal_get_image(cr)->lambda;

		default: return 0.0;

	}
}


static int check_angle_shifts(gsl_vector *cur, gsl_vector *initial,
                              enum gparam *rv, int n_params,
                              struct rf_priv *residual_f_priv)
{
	int i;
	double ang = 0.0;

	for ( i=0; i<n_params; i++ ) {
		if ( (rv[i] == GPARAM_ANG1) || (rv[i] == GPARAM_ANG2) ) {
			ang += fabs(get_actual_val(cur, initial, rv, i));
		}
	}

	if ( rad2deg(ang) > 5.0 ) {
		ERROR("More than 5 degrees total rotation!\n");
		residual_f_priv->verbose = 1;
		double res = residual_f(cur, residual_f_priv);
		STATUS("residual after rotation = %e\n", res);
		residual_f_priv->verbose = 2;
		res = residual_f(initial, residual_f_priv);
		STATUS("residual before rotation = %e\n", res);
		return 1;
	}
	return 0;
}


static void do_pr_refine(Crystal *cr, const RefList *full,
                         PartialityModel pmodel, int verbose, int serial)
{
	int i;
	gsl_multimin_fminimizer *min;
	gsl_vector *initial;
	gsl_vector *vals;
	gsl_vector *step;
	gsl_multimin_function f;
	enum gparam rv[32];
	struct rf_priv residual_f_priv;
	int n_params = 0;
	int n_iter = 0;
	int status;
	int r;
	double G;
	double residual_start, residual_free_start;

	residual_start = residual(cr, full, 0, NULL, NULL, 1);
	residual_free_start = residual(cr, full, 1, NULL, NULL, 1);

	if ( verbose ) {
		STATUS("\nPR initial: dev = %10.5e, free dev = %10.5e\n",
		       residual_start, residual_free_start);
	}

	/* The parameters to be refined */
	rv[n_params++] = GPARAM_ANG1;
	rv[n_params++] = GPARAM_ANG2;
	rv[n_params++] = GPARAM_R;
	rv[n_params++] = GPARAM_WAVELENGTH;

	residual_f_priv.cr = cr;
	residual_f_priv.full = full;
	residual_f_priv.rv = rv;
	residual_f_priv.verbose = verbose;
	residual_f_priv.serial = serial;
	f.f = residual_f;
	f.n = n_params;
	f.params = &residual_f_priv;

	initial = gsl_vector_alloc(n_params);
	vals = gsl_vector_alloc(n_params);
	step = gsl_vector_alloc(n_params);

	for ( i=0; i<n_params; i++ ) {
		gsl_vector_set(initial, i, get_initial_param(cr, rv[i]));
		gsl_vector_set(vals, i, 0.0);
		gsl_vector_set(step, i, 1.0);
	}

	residual_f_priv.initial = initial;
	min = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,
	                                    n_params);
	gsl_multimin_fminimizer_set(min, &f, vals, step);

	if ( verbose ) {
		double res = residual_f(min->x, &residual_f_priv);
		double size = gsl_multimin_fminimizer_size(min);
		STATUS("At start: %f %f %f %f ----> %f %f %e %f residual = %e size %f\n",
		       gsl_vector_get(min->x, 0),
		       gsl_vector_get(min->x, 1),
		       gsl_vector_get(min->x, 2),
		       gsl_vector_get(min->x, 3),
		       rad2deg(get_actual_val(min->x, initial, rv, 0)),
		       rad2deg(get_actual_val(min->x, initial, rv, 1)),
		       get_actual_val(min->x, initial, rv, 2),
		       get_actual_val(min->x, initial, rv, 3)*1e10,
		       res, size);
	}

	do {
		n_iter++;

		status = gsl_multimin_fminimizer_iterate(min);
		if ( status ) break;

		if ( verbose ) {
			double res = residual_f(min->x, &residual_f_priv);
			double size = gsl_multimin_fminimizer_size(min);
			STATUS("%f %f %f %f ----> %f %f %e %f residual = %e size %f\n",
			       gsl_vector_get(min->x, 0),
			       gsl_vector_get(min->x, 1),
			       gsl_vector_get(min->x, 2),
			       gsl_vector_get(min->x, 3),
			       rad2deg(get_actual_val(min->x, initial, rv, 0)),
			       rad2deg(get_actual_val(min->x, initial, rv, 1)),
			       get_actual_val(min->x, initial, rv, 2),
			       get_actual_val(min->x, initial, rv, 3)*1e10,
			       res, size);
		}

		status = gsl_multimin_test_size(min->size, 0.005);

	} while ( status == GSL_CONTINUE && n_iter < 1000 );

	if ( verbose ) {
		STATUS("Done with refinement after %i iter\n", n_iter);
		STATUS("status = %i (%s)\n", status, gsl_strerror(status));
	}

	if ( check_angle_shifts(min->x, initial, rv, n_params, &residual_f_priv) ) return;

	if ( verbose ) {

		double res = residual_f(min->x, &residual_f_priv);
		double size = gsl_multimin_fminimizer_size(min);
		STATUS("At end: %f %f %f %f ----> %f %f %e %f residual = %e size %f\n",
		       gsl_vector_get(min->x, 0),
		       gsl_vector_get(min->x, 1),
		       gsl_vector_get(min->x, 2),
		       gsl_vector_get(min->x, 3),
		       rad2deg(get_actual_val(min->x, initial, rv, 0)),
		       rad2deg(get_actual_val(min->x, initial, rv, 1)),
		       get_actual_val(min->x, initial, rv, 2),
		       get_actual_val(min->x, initial, rv, 3)*1e10,
		       res, size);

	}

	/* Apply the final shifts */
	apply_parameters(min->x, initial, rv, cr);
	update_predictions(cr);
	calculate_partialities(cr, PMODEL_XSPHERE);
	r = linear_scale(full, crystal_get_reflections(cr), &G, 0);
	if ( r == 0 ) {
		crystal_set_osf(cr, G);
	}

	if ( verbose ) {

		STATUS("After applying final shifts:\n");
		STATUS("PR final: dev = %10.5e, free dev = %10.5e\n",
		       residual(cr, full, 0, NULL, NULL, 0),
		       residual(cr, full, 1, NULL, NULL, 0));
		STATUS("Final R = %e m^-1\n", crystal_get_profile_radius(cr));

	}

	gsl_multimin_fminimizer_free(min);
	gsl_vector_free(initial);
	gsl_vector_free(vals);
	gsl_vector_free(step);
}


static struct prdata pr_refine(Crystal *cr, const RefList *full,
                               PartialityModel pmodel, int verbose, int serial)
{
	struct prdata prdata;

	prdata.refined = 0;

	do_pr_refine(cr, full, pmodel, verbose, serial);

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
	int serial;
	struct prdata prdata;
	int verbose;
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

	pargs->prdata = pr_refine(cr, pargs->full, pargs->pmodel, pargs->verbose, pargs->serial);
}


static void *get_image(void *vqargs)
{
	struct refine_args *task;
	struct queue_args *qargs = vqargs;

	task = malloc(sizeof(struct refine_args));
	memcpy(task, &qargs->task_defaults, sizeof(struct refine_args));

	task->crystal = qargs->crystals[qargs->n_started];
	task->serial = qargs->n_started;

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
                RefList *full, int nthreads, PartialityModel pmodel, int verbose)
{
	struct refine_args task_defaults;
	struct queue_args qargs;

	task_defaults.full = full;
	task_defaults.crystal = NULL;
	task_defaults.pmodel = pmodel;
	task_defaults.prdata.refined = 0;
	task_defaults.verbose = verbose;

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
