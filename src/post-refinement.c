/*
 * post-refinement.c
 *
 * Post refinement
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
#include <gsl/gsl_blas.h>

#include "image.h"
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"
#include "cell.h"
#include "cell-utils.h"


/* Maximum number of iterations of NLSq to do for each image per macrocycle. */
#define MAX_CYCLES (10)

struct prdata
{
	int refined;
	int n_filtered;
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


/* Returns dp(gauss)/dr at "r" */
static double gaussian_fraction_gradient(double r, double R)
{
	const double ng = 2.6;
	const double sig = R/ng;

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return exp(-pow(r/sig, 2.0)/2.0) / (sig*sqrt(2.0*M_PI));
}


/* Returns dp(sph)/dr at "r" */
static double sphere_fraction_gradient(double r, double pr)
{
	double q, dpdq, dqdr;

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -pr ) return 0.0;
	if ( r > +pr ) return 0.0;

	q = (r + pr)/(2.0*pr);
	dpdq = 6.0*(q - q*q);
	dqdr = 1.0 / (2.0*pr);
	return dpdq * dqdr;
}


/* Returns dp/dr at "r" */
static double partiality_gradient(double r, double pr,
                                  PartialityModel pmodel,
                                  double rlow, double rhigh)
{
	double A, D;

	D = rlow - rhigh;

	switch ( pmodel ) {

		default:
		case PMODEL_UNITY:
		return 0.0;

		case PMODEL_SCSPHERE:
		A = sphere_fraction_gradient(r, pr)/D;
		return 4.0*pr*A/3.0;

		case PMODEL_SCGAUSSIAN:
		A = gaussian_fraction_gradient(r, pr)/D;
		return 4.0*pr*A/3.0;

	}
}


static double sphere_fraction_rgradient(double r, double R)
{
	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return -(3.0*r/(4.0*R*R)) * (1.0 - r*r/(R*R));
}


static double gaussian_fraction_rgradient(double r, double R)
{
	const double ng = 2.6;
	const double sig = R/ng;

	/* If the Ewald sphere isn't within the profile, the gradient is zero */
	if ( r < -R ) return 0.0;
	if ( r > +R ) return 0.0;

	return -(ng*r/(sqrt(2.0*M_PI)*R*R))*exp(-r*r/(2.0*sig*sig));
}


static double volume_fraction_rgradient(double r, double pr,
                                       PartialityModel pmodel)
{
	switch ( pmodel )
	{
		case PMODEL_UNITY :
		return 1.0;

		case PMODEL_SCSPHERE :
		return sphere_fraction_rgradient(r, pr);

		case PMODEL_SCGAUSSIAN :
		return gaussian_fraction_rgradient(r, pr);

		default :
		ERROR("No pmodel in volume_fraction_rgradient!\n");
		return 1.0;
	}
}


static double volume_fraction(double rlow, double rhigh, double pr,
                              PartialityModel pmodel)
{
	switch ( pmodel )
	{
		case PMODEL_UNITY :
		return 1.0;

		case PMODEL_SCSPHERE :
		return sphere_fraction(rlow, rhigh, pr);

		case PMODEL_SCGAUSSIAN :
		return gaussian_fraction(rlow, rhigh, pr);

		default :
		ERROR("No pmodel in volume_fraction!\n");
		return 1.0;
	}
}


/* Return the gradient of "fx" wrt parameter 'k' given the current
 * status of the crystal. */
double gradient(Crystal *cr, int k, Reflection *refl, PartialityModel pmodel)
{
	double glow, ghigh;
	double rlow, rhigh, p;
	struct image *image = crystal_get_image(cr);
	double R = crystal_get_profile_radius(cr);
	double gr;

	get_partial(refl, &rlow, &rhigh, &p);

	if ( k == GPARAM_R ) {

		double Rglow, Rghigh;
		double D, psph;

		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);

		Rglow = volume_fraction_rgradient(rlow, R, pmodel);
		Rghigh = volume_fraction_rgradient(rhigh, R, pmodel);

		gr = 4.0*psph/(3.0*D) + (4.0*R/(3.0*D))*(Rglow - Rghigh);
		return gr;

	}

	/* Calculate the gradient of partiality wrt excitation error. */
	glow = partiality_gradient(rlow, R, pmodel, rlow, rhigh);
	ghigh = partiality_gradient(rhigh, R, pmodel, rlow, rhigh);

	if ( k == GPARAM_DIV ) {

		double D, psph, ds;
		signed int hs, ks, ls;

		D = rlow - rhigh;
		psph = volume_fraction(rlow, rhigh, R, pmodel);
		get_symmetric_indices(refl, &hs, &ks, &ls);
		ds = 2.0 * resolution(crystal_get_cell(cr), hs, ks, ls);

		gr = (ds/2.0)*(glow+ghigh) - 4.0*R*psph*ds/(3.0*D*D);
		return gr;

	}

	gr = r_gradient(crystal_get_cell(cr), k, refl, image) * (glow-ghigh);
	return gr;
}


static void apply_cell_shift(UnitCell *cell, int k, double shift)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double as, bs, cs;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	as = modulus(asx, asy, asz);
	bs = modulus(bsx, bsy, bsz);
	cs = modulus(csx, csy, csz);

	/* Forbid any step which looks too large */
	switch ( k )
	{
		case GPARAM_ASX :
		case GPARAM_ASY :
		case GPARAM_ASZ :
		if ( fabs(shift) > 0.1*as ) return;
		break;

		case GPARAM_BSX :
		case GPARAM_BSY :
		case GPARAM_BSZ :
		if ( fabs(shift) > 0.1*bs ) return;
		break;

		case GPARAM_CSX :
		case GPARAM_CSY :
		case GPARAM_CSZ :
		if ( fabs(shift) > 0.1*cs ) return;
		break;
	}

	switch ( k )
	{
		case GPARAM_ASX :  asx += shift;  break;
		case GPARAM_ASY :  asy += shift;  break;
		case GPARAM_ASZ :  asz += shift;  break;
		case GPARAM_BSX :  bsx += shift;  break;
		case GPARAM_BSY :  bsy += shift;  break;
		case GPARAM_BSZ :  bsz += shift;  break;
		case GPARAM_CSX :  csx += shift;  break;
		case GPARAM_CSY :  csy += shift;  break;
		case GPARAM_CSZ :  csz += shift;  break;
	}

	cell_set_reciprocal(cell, asx, asy, asz,
	                          bsx, bsy, bsz,
	                          csx, csy, csz);
}


/* Apply the given shift to the 'k'th parameter of 'image'. */
static void apply_shift(Crystal *cr, int k, double shift)
{
	double t;
	struct image *image = crystal_get_image(cr);

	switch ( k ) {

		case GPARAM_DIV :
		if ( shift > 0.1*image->div ) return;
		image->div += shift;
		if ( image->div < 0.0 ) image->div = 0.0;
		break;

		case GPARAM_R :
		t = crystal_get_profile_radius(cr);
		if ( shift > 0.1*t ) return;
		t += shift;
		crystal_set_profile_radius(cr, t);
		break;

		case GPARAM_ASX :
		case GPARAM_ASY :
		case GPARAM_ASZ :
		case GPARAM_BSX :
		case GPARAM_BSY :
		case GPARAM_BSZ :
		case GPARAM_CSX :
		case GPARAM_CSY :
		case GPARAM_CSZ :
		apply_cell_shift(crystal_get_cell(cr), k, shift);
		break;

		default :
		ERROR("No shift defined for parameter %i\n", k);
		abort();

	}
}


/* Perform one cycle of post refinement on 'image' against 'full' */
static double pr_iterate(Crystal *cr, const RefList *full,
                         PartialityModel pmodel,
                         int *n_filtered, int verbose)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int param;
	Reflection *refl;
	RefListIterator *iter;
	RefList *reflections;
	double max_shift;
	int nref = 0;
	int num_params = 0;
	enum gparam rv[32];
	double G, B;

	if ( n_filtered != NULL ) *n_filtered = 0;

	rv[num_params++] = GPARAM_ASX;
	rv[num_params++] = GPARAM_ASY;
	rv[num_params++] = GPARAM_ASZ;
	rv[num_params++] = GPARAM_BSX;
	rv[num_params++] = GPARAM_BSY;
	rv[num_params++] = GPARAM_BSZ;
	rv[num_params++] = GPARAM_CSX;
	rv[num_params++] = GPARAM_CSY;
	rv[num_params++] = GPARAM_CSZ;
//	rv[num_params++] = GPARAM_R;
//	rv[num_params++] = GPARAM_DIV;

	M = gsl_matrix_calloc(num_params, num_params);
	v = gsl_vector_calloc(num_params);

	reflections = crystal_get_reflections(cr);
	G = crystal_get_osf(cr);
	B = crystal_get_Bfac(cr);

	/* Post-refinement terms */
	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int ha, ka, la;
		double I_full, delta_I, esd;
		double w;
		double I_partial;
		int k;
		double p, L, s;
		double fx;
		Reflection *match;
		double gradients[num_params];

		if ( get_flag(refl) ) continue;

		get_indices(refl, &ha, &ka, &la);
		match = find_refl(full, ha, ka, la);
		if ( match == NULL ) continue;
		I_full = get_intensity(match);

		if ( get_redundancy(match) < 2 ) continue;

		p = get_partiality(refl);
		L = get_lorentz(refl);
		I_partial = get_intensity(refl);
		esd = get_esd_intensity(refl);
		s = resolution(crystal_get_cell(cr), ha, ka, la);

		if ( I_partial < 3.0*esd ) continue;

		/* Calculate the weight for this reflection */
		w = (s/1e9)*(s/1e9) / (esd*esd);

		/* Calculate all gradients for this reflection */
		for ( k=0; k<num_params; k++ ) {
			gradients[k] = gradient(cr, rv[k], refl, pmodel);
			gradients[k] *= exp(G)*exp(-B*s*s)*I_full/L;
		}

		for ( k=0; k<num_params; k++ ) {

			int g;
			double v_c, v_curr;

			for ( g=0; g<num_params; g++ ) {

				double M_c, M_curr;

				/* Matrix is symmetric */
				if ( g > k ) continue;

				M_c = w * gradients[g] * gradients[k];

				M_curr = gsl_matrix_get(M, k, g);
				gsl_matrix_set(M, k, g, M_curr + M_c);
				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			fx = exp(G)*p*exp(-B*s*s)*I_full/L;
			delta_I = I_partial - fx;
			v_c = w * delta_I * gradients[k];
			v_curr = gsl_vector_get(v, k);
			gsl_vector_set(v, k, v_curr + v_c);

		}

		nref++;
	}
	if ( verbose ) {
		//STATUS("The original equation:\n");
		//show_matrix_eqn(M, v);
		STATUS("%i reflections went into the equations.\n", nref);
	}

	if ( nref < num_params ) {
		crystal_set_user_flag(cr, PRFLAG_FEWREFL);
		gsl_matrix_free(M);
		gsl_vector_free(v);
		return 0.0;
	}

	max_shift = 0.0;
	shifts = solve_svd(v, M, n_filtered, 0);
	if ( shifts != NULL ) {

		for ( param=0; param<num_params; param++ ) {
			double shift = gsl_vector_get(shifts, param);
			if ( verbose ) STATUS("Shift %i: %e\n", param, shift);
			if ( isnan(shift) ) {
				//ERROR("NaN shift parameter %i (image ser %i), "
				//       "%i reflections.\n", rv[param],
				//       crystal_get_image(cr)->serial,
				//       nref);
			} else {
				apply_shift(cr, rv[param], shift);
			}
			if ( fabs(shift) > max_shift ) max_shift = fabs(shift);
		}

	} else {
		crystal_set_user_flag(cr, PRFLAG_SOLVEFAIL);
	}

	gsl_matrix_free(M);
	gsl_vector_free(v);
	gsl_vector_free(shifts);

	return max_shift;
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
		fh = fopen(filename, "a");
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


static void write_residual_graph(Crystal *cr, const RefList *full)
{
	int i;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a;
	double step;
	double orig_asx;
	FILE *fh;
	UnitCell *cell;

	cell = crystal_get_cell(cr);

	fh = fopen("residual-plot.dat", "w");

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	a = modulus(asx, asy, asz);
	orig_asx = asx;

	step = a/100.0;

	for ( i=-50; i<=50; i++ ) {

		double res;
		int n;
		asx = orig_asx + (i*step);
		cell_set_reciprocal(cell, asx, asy, asz,
		                          bsx, bsy, bsz,
		                          csx, csy, csz);
		update_partialities(cr, PMODEL_SCSPHERE);
		res = residual(cr, full, 0, &n, NULL);
		fprintf(fh, "%i %e %e %i\n", i, asx, res, n);
	}

	cell_set_reciprocal(cell, orig_asx, asy, asz,
	                          bsx, bsy, bsz,
	                          csx, csy, csz);
	update_partialities(cr, PMODEL_SCSPHERE);
	fclose(fh);
}


static void do_pr_refine(Crystal *cr, const RefList *full,
                         PartialityModel pmodel, int verbose)
{
	int i, done;
	double old_dev;
	UnitCell *cell = crystal_get_cell(cr);

	old_dev = residual(cr, full, 0, 0, NULL);

	if ( verbose ) {
		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		cell_get_reciprocal(cell, &asx, &asy, &asz,
		                          &bsx, &bsy, &bsz,
		                          &csx, &csy, &csz);
		STATUS("Initial asx = %e\n", asx);
		STATUS("PR initial  dev =  %10.5e, free dev = %10.5e\n",
		       old_dev, residual(cr, full, 1, NULL, NULL));
	}

	i = 0;
	done = 0;
	do {

		double dev;

		pr_iterate(cr, full, pmodel, NULL, verbose);

		update_partialities(cr, pmodel);

		dev = residual(cr, full, 0, 0, NULL);
		if ( fabs(dev - old_dev) < dev*0.0001 ) done = 1;

		if ( verbose ) {
			STATUS("PR iter %2i: dev = %10.5e, free dev = %10.5e\n",
			       i+1, dev, residual(cr, full, 1, NULL, NULL));
		}

		i++;
		old_dev = dev;

	} while ( i < 30 && !done );

	if ( verbose ) {
		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		cell_get_reciprocal(cell, &asx, &asy, &asz,
		                          &bsx, &bsy, &bsz,
		                          &csx, &csy, &csz);
		STATUS("Final asx = %e\n", asx);
	}
}


static struct prdata pr_refine(Crystal *cr, const RefList *full,
                               PartialityModel pmodel)
{
	int verbose = 0;
	struct prdata prdata;

	prdata.refined = 0;
	prdata.n_filtered = 0;

	if ( verbose ) {
		write_residual_graph(cr, full);
	}

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
	task_defaults.prdata.n_filtered = 0;

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
