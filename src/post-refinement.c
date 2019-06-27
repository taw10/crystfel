/*
 * post-refinement.c
 *
 * Post refinement
 *
 * Copyright © 2012-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2018 Thomas White <taw@physics.org>
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

#include "image.h"
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "cell-utils.h"
#include "reflist-utils.h"
#include "scaling.h"
#include "merge.h"

struct rf_alteration
{
	double rot_x;
	double rot_y;
	double delta_R;
	double delta_wave;
};


struct rf_priv
{
	const Crystal *cr;      /**< Original crystal (before any refinement) */
	const RefList *full;
	int serial;
	int scaleflags;

	Crystal *cr_tgt;        /**< Crystal to use for testing modifications */
	struct image image_tgt; /**< Image structure to go with cr_tgt */
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

		case PRFLAG_DELTACCHALF :
		return "negative delta CC½";

		case PRFLAG_BIGB :
		return "B too big";

		case PRFLAG_SCALEBAD :
		return "bad scaling";

		default :
		return "Unknown flag";
	}
}


static void rotate_cell_xy(const UnitCell *source, UnitCell *tgt,
                           double ang1, double ang2)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xnew, ynew, znew;

	cell_get_reciprocal(source, &asx, &asy, &asz,
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

	cell_set_reciprocal(tgt, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);
}


static void apply_parameters(const Crystal *cr_orig, Crystal *cr_tgt,
                             struct rf_alteration alter)
{
	double R, lambda;

	R = crystal_get_profile_radius(cr_orig);
	lambda = crystal_get_image_const(cr_orig)->lambda;

	rotate_cell_xy(crystal_get_cell_const(cr_orig), crystal_get_cell(cr_tgt),
	               alter.rot_x, alter.rot_y);
	crystal_set_profile_radius(cr_tgt, R+alter.delta_R);
	crystal_get_image(cr_tgt)->lambda = lambda+alter.delta_wave;
}


static double calc_residual(struct rf_priv *pv, struct rf_alteration alter,
                            int free)
{
	apply_parameters(pv->cr, pv->cr_tgt, alter);

	if ( fabs(crystal_get_profile_radius(pv->cr_tgt)) > 5e9 ) {
		STATUS("radius > 5e9\n");
		return NAN;
	}

	/* Can happen with grid scans and certain --force-radius values */
	if ( fabs(crystal_get_profile_radius(pv->cr_tgt)) < 0.0000001e9 ) {
		STATUS("radius very small\n");
		return NAN;
	}

	if ( crystal_get_image(pv->cr_tgt)->lambda <= 0.0 ) {
		STATUS("lambda < 0\n");
		return NAN;
	}

	update_predictions(pv->cr_tgt);
	calculate_partialities(pv->cr_tgt, PMODEL_XSPHERE);


	return residual(pv->cr_tgt, pv->full, free, NULL, NULL);
}


static RefList *reindex_reflections(RefList *input, SymOpList *amb,
                                    SymOpList *sym, int idx)
{
	RefList *n;
	Reflection *refl;
	RefListIterator *iter;

	n = reflist_new();

	for ( refl = first_refl(input, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		Reflection *rn;

		get_indices(refl, &h, &k, &l);
		get_equiv(amb, NULL, idx, h, k, l, &h, &k, &l);
		get_asymm(sym, h, k, l, &h, &k, &l);
		rn = add_refl(n, h, k, l);
		copy_data(rn, refl);

		get_symmetric_indices(rn, &h, &k, &l);
		get_equiv(amb, NULL, idx, h, k, l, &h, &k, &l);
		set_symmetric_indices(rn, h, k, l);
	}

	return n;
}


static void reindex_cell(UnitCell *cell, SymOpList *amb, int idx)
{
	signed int h, k, l;
	struct rvec na, nb, nc;
	struct rvec as, bs, cs;

	cell_get_reciprocal(cell, &as.u, &as.v, &as.w,
	                          &bs.u, &bs.v, &bs.w,
	                          &cs.u, &cs.v, &cs.w);

	get_equiv(amb, NULL, idx, 1, 0, 0, &h, &k, &l);
	na.u = as.u*h + bs.u*k + cs.u*l;
	na.v = as.v*h + bs.v*k + cs.v*l;
	na.w = as.w*h + bs.w*k + cs.w*l;

	get_equiv(amb, NULL, idx, 0, 1, 0, &h, &k, &l);
	nb.u = as.u*h + bs.u*k + cs.u*l;
	nb.v = as.v*h + bs.v*k + cs.v*l;
	nb.w = as.w*h + bs.w*k + cs.w*l;

	get_equiv(amb, NULL, idx, 0, 0, 1, &h, &k, &l);
	nc.u = as.u*h + bs.u*k + cs.u*l;
	nc.v = as.v*h + bs.v*k + cs.v*l;
	nc.w = as.w*h + bs.w*k + cs.w*l;

	cell_set_reciprocal(cell, na.u, na.v, na.w,
	                          nb.u, nb.v, nb.w,
	                          nc.u, nc.v, nc.w);
}


static void try_reindex(Crystal *crin, const RefList *full,
                        SymOpList *sym, SymOpList *amb, int scaleflags)
{
	RefList *list;
	Crystal *cr;
	UnitCell *cell;
	double residual_original, residual_flipped;
	int idx, n;

	if ( sym == NULL || amb == NULL ) return;

	if ( scale_one_crystal(crin, full, scaleflags) ) return;
	residual_original = residual(crin, full, 0, NULL, NULL);

	cr = crystal_copy(crin);

	n = num_equivs(amb, NULL);

	for ( idx=0; idx<n; idx++ ) {

		cell = cell_new_from_cell(crystal_get_cell(crin));
		if ( cell == NULL ) return;
		reindex_cell(cell, amb, idx);
		crystal_set_cell(cr, cell);

		list = reindex_reflections(crystal_get_reflections(crin),
		                           amb, sym, idx);
		crystal_set_reflections(cr, list);

		update_predictions(cr);
		calculate_partialities(cr, PMODEL_XSPHERE);

		if ( scale_one_crystal(cr, full, scaleflags) ) return;
		residual_flipped = residual(cr, full, 0, NULL, NULL);

		if ( residual_flipped < residual_original ) {
			crystal_set_cell(crin, cell);
			crystal_set_reflections(crin, list);
			residual_original = residual_flipped;
		} else {
			cell_free(crystal_get_cell(cr));
			reflist_free(crystal_get_reflections(cr));
		}
	}

	crystal_free(cr);
}


void write_test_logs(Crystal *crystal, const RefList *full,
                     signed int cycle, int serial)
{
	FILE *fh;
	struct image *image = crystal_get_image(crystal);
	char tmp[256];
	char ins[16];

	snprintf(tmp, 256, "pr-logs/parameters-crystal%i.dat", serial);

	if ( cycle == 0 ) {
		fh = fopen(tmp, "w");
	} else {
		fh = fopen(tmp, "a");
	}

	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", tmp);
		return;
	}

	if ( cycle == 0 ) {
		char *evstr = get_event_string(image->event);
		fprintf(fh, "Image: %s %s\n", image->filename, evstr);
		free(evstr);
	}

	if ( cycle >= 0 ) {
		snprintf(ins, 16, "%i", cycle);
	} else {
		ins[0] = 'F';
		ins[1] = '\0';
	}

	fprintf(fh, "%s rlp_size = %e\n", ins, crystal_get_profile_radius(crystal)/1e10);
	fprintf(fh, "%s mosaicity = %e\n", ins, crystal_get_mosaicity(crystal));
	fprintf(fh, "%s wavelength = %f A\n", ins, image->lambda*1e10);
	fprintf(fh, "%s bandwidth = %f\n", ins, image->bw);
	fprintf(fh, "%s my scale factor = %e\n", ins, crystal_get_osf(crystal));

	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	cell_get_reciprocal(crystal_get_cell(crystal), &asx, &asy, &asz,
	                                               &bsx, &bsy, &bsz,
	                                               &csx, &csy, &csz);
	fprintf(fh, "%s %f, %f, %f, %f, %f, %f, %f, %f, %f\n", ins,
	        asx/1e10, bsx/1e10, csx/1e10,
	        asy/1e10, bsy/1e10, csy/1e10,
	        asz/1e10, bsz/1e10, csz/1e10);
	fclose(fh);
}


void write_specgraph(Crystal *crystal, const RefList *full,
                     signed int cycle, int serial)
{
	FILE *fh;
	char tmp[256];
	Reflection *refl;
	RefListIterator *iter;
	double G = crystal_get_osf(crystal);
	double B = crystal_get_Bfac(crystal);
	UnitCell *cell;
	struct image *image = crystal_get_image(crystal);
	char ins[16];

	snprintf(tmp, 256, "pr-logs/specgraph-crystal%i.dat", serial);

	if ( cycle == 0 ) {
		fh = fopen(tmp, "w");
	} else {
		fh = fopen(tmp, "a");
	}

	if ( fh == NULL ) {
		ERROR("Failed to open '%s'\n", tmp);
		return;
	}

	if ( cycle == 0 ) {
		char *evstr = get_event_string(image->event);
		fprintf(fh, "Image: %s %s\n", image->filename, evstr);
		fprintf(fh, "khalf/m   1/d(m)  pcalc    pobs   iteration  h  k  l\n");
		free(evstr);
	}

	cell = crystal_get_cell(crystal);

	if ( cycle >= 0 ) {
		snprintf(ins, 16, "%i", cycle);
	} else {
		ins[0] = 'F';
		ins[1] = '\0';
	}

	for ( refl = first_refl(crystal_get_reflections(crystal), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double Ipart, Ifull, pobs, pcalc;
		double res;
		signed int h, k, l;
		Reflection *match;

		/* Strong reflections only */
		if ( get_intensity(refl) < 3.0*get_esd_intensity(refl) ) continue;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);

		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;

		/* Don't calculate pobs if reference reflection is weak */
		if ( fabs(get_intensity(match)) / get_esd_intensity(match) < 3.0 ) continue;

		Ipart = correct_reflection_nopart(get_intensity(refl), refl, G, B, res);
		Ifull = get_intensity(match);
		pobs = Ipart / Ifull;
		pcalc = get_partiality(refl);

		fprintf(fh, "%e   %e   %f   %f   %s  %4i %4i %4i\n",
		        get_khalf(refl), 2.0*res, pcalc, pobs, ins, h, k, l);

	}

	fclose(fh);
}


static void write_angle_grid(Crystal *cr, const RefList *full,
                             signed int cycle, int serial, int scaleflags)
{
	FILE *fh;
	char fn[64];
	char ins[16];
	struct rf_priv priv;
	RefList *list;
	UnitCell *cell;

	priv.cr = cr;
	priv.full = full;
	priv.serial = serial;
	priv.scaleflags = scaleflags;
	priv.cr_tgt = crystal_copy(cr);
	priv.image_tgt = *crystal_get_image(cr);
	crystal_set_image(priv.cr_tgt, &priv.image_tgt);
	list = copy_reflist(crystal_get_reflections(cr));
	crystal_set_reflections(priv.cr_tgt, list);
	cell = cell_new_from_cell(crystal_get_cell(cr));
	crystal_set_cell(priv.cr_tgt, cell);

	if ( cycle >= 0 ) {
		snprintf(ins, 16, "%i", cycle);
	} else {
		ins[0] = 'F';
		ins[1] = '\0';
	}

	snprintf(fn, 64, "pr-logs/grid-crystal%i-cycle%s-ang1-ang2.dat",
	         serial, ins);
	fh = fopen(fn, "w");
	if ( fh != NULL ) {
		double v1, v2;
		fprintf(fh, "%e %e %e %s\n", -5.0e-3, 5.0e-3, 0.0, "rot_x/rad");
		fprintf(fh, "%e %e %e %s\n", -5.0e-3, 5.0e-3, 0.0, "rot_y/rad");
		for ( v2=-5.0e-3; v2<=5.1e-3; v2+=0.25e-3 ) {
			int first=1;
			for ( v1=-5.0e-3; v1<=5.1e-3; v1+=0.25e-3 ) {
				double res;
				struct rf_alteration alter;
				alter.rot_x = v1;
				alter.rot_y = v2;
				alter.delta_R = 0.0;
				alter.delta_wave = 0.0;
				res = calc_residual(&priv, alter, 0);
				if ( !first ) fprintf(fh, " ");
				first = 0;
				fprintf(fh, "%e", res);
			}
			fprintf(fh, "\n");
		}
	}
	fclose(fh);

	reflist_free(crystal_get_reflections(priv.cr_tgt));
	crystal_free(priv.cr_tgt);
}


static void write_radius_grid(Crystal *cr, const RefList *full,
                              signed int cycle, int serial, int scaleflags)
{
	FILE *fh;
	char fn[64];
	char ins[16];
	struct rf_priv priv;
	RefList *list;
	UnitCell *cell;

	priv.cr = cr;
	priv.full = full;
	priv.serial = serial;
	priv.scaleflags = scaleflags;
	priv.cr_tgt = crystal_copy(cr);
	priv.image_tgt = *crystal_get_image(cr);
	crystal_set_image(priv.cr_tgt, &priv.image_tgt);
	list = copy_reflist(crystal_get_reflections(cr));
	crystal_set_reflections(priv.cr_tgt, list);
	cell = cell_new_from_cell(crystal_get_cell(cr));
	crystal_set_cell(priv.cr_tgt, cell);

	if ( cycle >= 0 ) {
		snprintf(ins, 16, "%i", cycle);
	} else {
		ins[0] = 'F';
		ins[1] = '\0';
	}

	snprintf(fn, 64, "pr-logs/grid-crystal%i-cycle%s-R-wave.dat",
	         serial, ins);
	fh = fopen(fn, "w");
	if ( fh != NULL ) {
		double v1, v2;
		fprintf(fh, "%e %e %e %s\n", -4e-13, 4e-13, 0.0, "wavelength change/m");
		fprintf(fh, "%e %e %e %s\n", -2e6, 2e6, 0.0, "radius change/m^-1");
		for ( v2=-2e6; v2<=2.001e6; v2+=100000 ) {
			int first=1;
			for ( v1=-4e-13; v1<=4.001e-13; v1+=2e-14 ) {
				double res;
				struct rf_alteration alter;
				alter.rot_x = 0.0;
				alter.rot_y = 0.0;
				alter.delta_R = v2;
				alter.delta_wave = v1;
				res = calc_residual(&priv, alter, 0);
				if ( !first ) fprintf(fh, " ");
				first = 0;
				fprintf(fh, "%e", res);
			}
			fprintf(fh, "\n");
		}
	}
	fclose(fh);

	reflist_free(crystal_get_reflections(priv.cr_tgt));
	crystal_free(priv.cr_tgt);
}


void write_gridscan(Crystal *cr, const RefList *full,
                    signed int cycle, int serial, int scaleflags)
{
	write_angle_grid(cr, full, cycle, serial, scaleflags);
	write_radius_grid(cr, full, cycle, serial, scaleflags);
}


static void do_pr_refine(Crystal *cr, const RefList *full,
                         PartialityModel pmodel, int serial,
                         int cycle, int write_logs,
                         SymOpList *sym, SymOpList *amb, int scaleflags)
{
	struct rf_priv priv;
	struct rf_alteration alter;
	int n_iter = 0;
	int status = 0;
	double fom, freefom;
	RefList *list;
	FILE *fh = NULL;
	UnitCell *cell;

	try_reindex(cr, full, sym, amb, scaleflags);

	if ( scale_one_crystal(cr, full, scaleflags | SCALE_VERBOSE_ERRORS) ) {
		ERROR("Bad scaling at start of refinement.\n");
		return;
	}

	alter.rot_x = 0.0;
	alter.rot_y = 0.0;
	alter.delta_R = 0.0;
	alter.delta_wave = 0.0;

	priv.cr = cr;
	priv.full = full;
	priv.serial = serial;
	priv.scaleflags = scaleflags;
	priv.cr_tgt = crystal_copy(cr);
	priv.image_tgt = *crystal_get_image(cr);
	crystal_set_image(priv.cr_tgt, &priv.image_tgt);
	list = copy_reflist(crystal_get_reflections(cr));
	crystal_set_reflections(priv.cr_tgt, list);
	cell = cell_new_from_cell(crystal_get_cell(cr));
	crystal_set_cell(priv.cr_tgt, cell);

	fom = calc_residual(&priv, alter, 0);
	freefom = calc_residual(&priv, alter, 1);

	if ( write_logs ) {

		char fn[64];

		snprintf(fn, 63, "pr-logs/crystal%i-cycle%i.log", serial, cycle);
		fh = fopen(fn, "w");
		if ( fh != NULL ) {
			fprintf(fh, "iter  FoM        FreeFoM     rotx/rad   "
			            "roty/rad    radius/m      wavelength/m\n");
			fprintf(fh, "%5i %10.8f  %10.8f %10.8f %10.8f  %e  %e\n",
			        n_iter, fom, freefom,
			        alter.rot_x, alter.rot_y,
			        crystal_get_profile_radius(cr)+alter.delta_R,
			        crystal_get_image(cr)->lambda+alter.delta_wave);
		}

	}

	do {
		n_iter++;

		/* FIXME: Actually minimise */

		fom = calc_residual(&priv, alter, 0);
		freefom = calc_residual(&priv, alter, 1);

		if ( fh != NULL ) {
			fprintf(fh, "%5i %10.8f  %10.8f %10.8f %10.8f  %e  %e\n",
			        n_iter, fom, freefom,
			        alter.rot_x, alter.rot_y,
			        crystal_get_profile_radius(cr)+alter.delta_R,
			        crystal_get_image(cr)->lambda+alter.delta_wave);
		}

	} while ( n_iter < 10 );  /* FIXME: Termination criteria */

	if ( status == 0 ) {

		/* Apply the final shifts */
		apply_parameters(cr, cr, alter);
		update_predictions(cr);
		calculate_partialities(cr, PMODEL_XSPHERE);
		scale_one_crystal(cr, full, scaleflags);

	} else {
		ERROR("Bad refinement: crystal %i after %i iterations\n",
		      serial, n_iter);
	}

	if ( write_logs ) {
		write_gridscan(cr, full, cycle, serial, scaleflags);
		write_specgraph(cr, full, cycle, serial);
		write_test_logs(cr, full, cycle, serial);
	}

	if ( crystal_get_profile_radius(cr) > 5e9 ) {
		ERROR("WARNING: Very large radius: crystal %i\n", serial);
	}

	if ( fh != NULL ) {
		fclose(fh);
	}

	reflist_free(crystal_get_reflections(priv.cr_tgt));
	crystal_free(priv.cr_tgt);
}


struct refine_args
{
	RefList *full;
	Crystal *crystal;
	PartialityModel pmodel;
	int serial;
	int cycle;
	int no_logs;
	SymOpList *sym;
	SymOpList *amb;
	int scaleflags;
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
	int write_logs = 0;

	write_logs = !pargs->no_logs && (pargs->serial % 20 == 0);

	do_pr_refine(cr, pargs->full, pargs->pmodel,
	             pargs->serial, pargs->cycle, write_logs,
	             pargs->sym, pargs->amb, pargs->scaleflags);
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
                RefList *full, int nthreads, PartialityModel pmodel,
                int cycle, int no_logs,
                SymOpList *sym, SymOpList *amb, int scaleflags)
{
	struct refine_args task_defaults;
	struct queue_args qargs;

	task_defaults.full = full;
	task_defaults.crystal = NULL;
	task_defaults.pmodel = pmodel;
	task_defaults.cycle = cycle;
	task_defaults.no_logs = no_logs;
	task_defaults.sym = sym;
	task_defaults.amb = amb;
	task_defaults.scaleflags = scaleflags;

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
