/*
 * error_model.c
 *
 * Error modelling
 *
 * Copyright © 2026 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2026 Thomas White <taw@physics.org>
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_multimin.h>

#include "utils.h"
#include "reflist.h"
#include "reflist-utils.h"
#include "merge.h"
#include "cell.h"
#include "cell-utils.h"
#include "error_model.h"


struct error_model
{
	ErrorModelType type;

	double sdfac;
	double sdb;
	double sdadd;
};


ErrorModel *error_model_new(ErrorModelType t)
{
	ErrorModel *e = malloc(sizeof(struct error_model));
	if ( e == NULL ) return NULL;
	e->type = t;

	e->sdfac = 1.0;
	e->sdb = 0.0;
	e->sdadd = 0.0;

	return e;
}


double merged_esd(struct reflection_contributions *c, double *weights,
                  double wmean, ErrorModel *emodel)
{
	if ( emodel->type == EMODEL_EQUIVS ) {
		double var;
		var = gsl_stats_wvariance_m(c->contrib_legacy_weights, 1,
		                            c->contribs, 1, c->n_contrib, wmean);
		return sqrt(var)/sqrt(c->n_contrib);
	} else {
		int i;
		double total = 0.0;
		for ( i=0; i<c->n_contrib; i++ ) {
			total += weights[i];
		}
		return sqrt(1.0/total);
	}
}


static double corr_esd(double sigij, double Ih, ErrorModel *emodel)
{
	double c;

	switch ( emodel->type ) {

		case EMODEL_EQUIVS:
		return sigij;

		case EMODEL_EV11:
		case EMODEL_KH23:
		c = sigij*sigij + emodel->sdb*emodel->sdb*Ih + emodel->sdadd*emodel->sdadd*Ih*Ih;
		return sqrt(emodel->sdfac*emodel->sdfac*c);

		case EMODEL_EV06:
		c = sigij*sigij + emodel->sdadd*emodel->sdadd*Ih*Ih;
		return sqrt(emodel->sdfac*emodel->sdfac*c);

	}
	abort();
}


double *make_weights_array(struct reflection_contributions *c, Reflection *refl, ErrorModel *emodel)
{
	int i;
	double *w;

	w = malloc(c->n_contrib * sizeof(double));

	for ( i=0; i<c->n_contrib; i++ ) {
		if ( emodel->type == EMODEL_EQUIVS ) {
			/* Preserve old behaviour for 'equivs'
			 * (no option provided for sigma-weighting without error modelling) */
			w[i] = c->contrib_legacy_weights[i];;
		} else if ( emodel->type == EMODEL_KH23 ) {
			double  esdIiC = corr_esd(c->contrib_esds[i], get_max_measurement(refl), emodel);
			w[i] = 1.0/(esdIiC*esdIiC);
		} else {
			double  esdIiC = corr_esd(c->contrib_esds[i], get_unweighted_mean(refl), emodel);
			w[i] = 1.0/(esdIiC*esdIiC);
		}
	}

	return w;
}


#define NQUANT (20)

static gsl_rstat_quantile_workspace **fill_quantiles(RefList *full, ErrorModel *emodel,
		                                     double *pminv, double *pmaxv)
{
	int i;
	Reflection *refl;
	RefListIterator *iter;
	gsl_rstat_quantile_workspace **quantiles;
	double minv = +INFINITY;
	double maxv = -INFINITY;

	quantiles = malloc(NQUANT*sizeof(gsl_rstat_quantile_workspace *));
	if ( quantiles == NULL ) return NULL;

	for ( i=0; i<NQUANT; i++ ) {
		double plotpos = ((double)i+1)/(NQUANT+1);
		quantiles[i] = gsl_rstat_quantile_alloc(plotpos);
		if ( quantiles[i] == NULL ) return NULL;
	}

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		int j;
		double Ih;
		struct reflection_contributions *c = get_contributions(refl);

		if ( c->n_contrib < 2 ) continue;

		if ( emodel->type == EMODEL_KH23 ) {
			/* Kh23 uses the highest intensity contribution in the
			 * error model equation, instead of the mean. */
			Ih = -INFINITY;
			for ( j=0; j<c->n_contrib; j++ ) {
				if ( c->contribs[j] > Ih ) Ih = c->contribs[j];
			}
		} else {
			Ih = get_intensity(refl);
		}

		for ( j=0; j<c->n_contrib; j++ ) {

			/* Mean (not max, for Kh23) without contribution j */
			double mIhj = (get_intensity(refl)*c->n_contrib - c->contribs[j])/(c->n_contrib-1);

			double bcorr = sqrt(((double)c->n_contrib-1)/c->n_contrib);
			double norm_dev = bcorr * (c->contribs[j] - mIhj)
			                   / corr_esd(c->contrib_esds[j], Ih, emodel);

			if ( norm_dev < -10 ) continue;
			if ( norm_dev > 10 ) continue;
			for ( i=0; i<NQUANT; i++ ) {
				gsl_rstat_quantile_add(norm_dev, quantiles[i]);
			}
			if ( norm_dev > maxv ) maxv = norm_dev;
			if ( norm_dev < minv ) minv = norm_dev;

		}
	}

	*pminv = minv;
	*pmaxv = maxv;
	return quantiles;
}


static double norm_res(ErrorModel *emodel, RefList *full)
{
	gsl_rstat_quantile_workspace **quantiles;
	double minv, maxv;
	int i;

	quantiles = fill_quantiles(full, emodel, &minv, &maxv);

	double total = 0.0;
	for ( i=0; i<NQUANT; i++ ) {
		double plotpos = ((double)i+1)/(NQUANT+1);
		total += pow(gsl_rstat_quantile_get(quantiles[i]) - gsl_cdf_gaussian_Pinv(plotpos, 1.0), 2.0);
		gsl_rstat_quantile_free(quantiles[i]);
	}
	free(quantiles);
	return total;
}


static void error_model_params_set_from_vector(ErrorModel *emodel, const gsl_vector *sdparams)
{
	switch ( emodel->type ) {

		case EMODEL_EQUIVS:
		break;

		case EMODEL_EV11:
		emodel->sdfac = fabs(gsl_vector_get(sdparams, 0));
		emodel->sdb   = fabs(gsl_vector_get(sdparams, 1));
		emodel->sdadd = fabs(gsl_vector_get(sdparams, 2));
		break;

		case EMODEL_EV06:
		emodel->sdfac = fabs(gsl_vector_get(sdparams, 0));
		emodel->sdadd = fabs(gsl_vector_get(sdparams, 1));
		break;

		case EMODEL_KH23:
		emodel->sdfac = fabs(gsl_vector_get(sdparams, 0));
		emodel->sdb   = fabs(gsl_vector_get(sdparams, 1));
		emodel->sdadd = fabs(gsl_vector_get(sdparams, 2));
		break;

	}
}


static double norm_res_equivs(const gsl_vector *sdparams, void *vp)
{
	ErrorModel emodel;
	RefList *full = vp;
	emodel.type = EMODEL_EQUIVS;
	return norm_res(&emodel, full);
}


static double norm_res_ev11(const gsl_vector *sdparams, void *vp)
{
	ErrorModel emodel;
	RefList *full = vp;
	emodel.type = EMODEL_EV11;
	error_model_params_set_from_vector(&emodel, sdparams);
	return norm_res(&emodel, full);
}


static double norm_res_ev06(const gsl_vector *sdparams, void *vp)
{
	ErrorModel emodel;
	RefList *full = vp;
	emodel.type = EMODEL_EV06;
	error_model_params_set_from_vector(&emodel, sdparams);
	return norm_res(&emodel, full);
}



static double norm_res_kh23(const gsl_vector *sdparams, void *vp)
{
	ErrorModel emodel;
	RefList *full = vp;
	emodel.type = EMODEL_KH23;
	error_model_params_set_from_vector(&emodel, sdparams);
	return norm_res(&emodel, full);
}


/* Return a wrapper function that can be used by GSL for minimisation */
static double (*error_model_norm_res_func(ErrorModelType t))(const gsl_vector *, void *)
{
	switch ( t ) {
		case EMODEL_EQUIVS: return norm_res_equivs;
		case EMODEL_EV11:   return norm_res_ev11;
		case EMODEL_EV06:   return norm_res_ev06;
		case EMODEL_KH23:   return norm_res_kh23;
	}
	abort();
}


void normal_probability_plot(RefList *full, ErrorModel *emodel)
{
	gsl_rstat_quantile_workspace **quantiles;
	double hstart;
	int i;
	double minv, maxv;

	quantiles = fill_quantiles(full, emodel, &minv, &maxv);
	printf("Bin start    Bin middle     Bin end        Density   Theoretical\n");
	printf("                        (=Sample quantile)            quantile  \n");
	printf("------------------------------------------------------------------\n");
	hstart = minv;
	printf("       -              -    %8.5f   %e             -\n", hstart, 0.0);
	for ( i=0; i<NQUANT; i++ ) {
		double plotpos = ((double)i+1)/(NQUANT+1);
		double hend = gsl_rstat_quantile_get(quantiles[i]);
		printf("%8.5f       %8.5f    %8.5f   %e      %8.5f\n",
		       hstart, hstart+(hend-hstart)/2.0, hend,
		       (1.0/(NQUANT+1))/(hend-hstart),
		       gsl_cdf_gaussian_Pinv(plotpos, 1.0));
		hstart = hend;
		gsl_rstat_quantile_free(quantiles[i]);
	}
	printf("%8.5f       %8.5f    %8.5f   %e             -\n",
	       hstart, hstart+(maxv-hstart)/2.0, maxv,
	       (1.0/(NQUANT+1))/(maxv-hstart));
	printf("\n\n");
	free(quantiles);
}


static int error_model_num_params(ErrorModel *emodel)
{
	switch ( emodel->type ) {
		case EMODEL_EQUIVS: return 0;
		case EMODEL_EV11:   return 3;
		case EMODEL_EV06:   return 2;
		case EMODEL_KH23:   return 3;
	}
	abort();
}


static gsl_vector *error_model_params_vector(ErrorModel *emodel)
{
	gsl_vector *sdparams = gsl_vector_alloc(error_model_num_params(emodel));

	switch ( emodel->type ) {

		case EMODEL_EQUIVS:
		break;

		case EMODEL_EV11:
		gsl_vector_set(sdparams, 0, emodel->sdfac);
		gsl_vector_set(sdparams, 1, emodel->sdb);
		gsl_vector_set(sdparams, 2, emodel->sdadd);
		break;

		case EMODEL_EV06:
		gsl_vector_set(sdparams, 0, emodel->sdfac);
		gsl_vector_set(sdparams, 1, emodel->sdadd);
		break;

		case EMODEL_KH23:
		gsl_vector_set(sdparams, 0, emodel->sdfac);
		gsl_vector_set(sdparams, 1, emodel->sdb);
		gsl_vector_set(sdparams, 2, emodel->sdadd);
		break;

	}

	return sdparams;
}


static gsl_vector *error_model_step_vector(ErrorModel *emodel)
{
	gsl_vector *stepsize = gsl_vector_alloc(error_model_num_params(emodel));

	switch ( emodel->type ) {

		case EMODEL_EQUIVS:
		break;

		case EMODEL_EV11:
		gsl_vector_set(stepsize, 0, 0.25);
		gsl_vector_set(stepsize, 1, 0.01);
		gsl_vector_set(stepsize, 2, 0.01);
		break;

		case EMODEL_EV06:
		gsl_vector_set(stepsize, 0, 0.25);
		gsl_vector_set(stepsize, 1, 0.01);
		break;

		case EMODEL_KH23:
		gsl_vector_set(stepsize, 0, 0.25);
		gsl_vector_set(stepsize, 1, 0.01);
		gsl_vector_set(stepsize, 2, 0.01);
		break;

	}

	return stepsize;
}


void refine_error_model(RefList *full, ErrorModel *emodel)
{
	gsl_multimin_fminimizer *mini;
	gsl_multimin_function myfunc;
	gsl_vector *sdparams;
	gsl_vector *stepsize;
	int r;
	int niter;

	if ( emodel->type == EMODEL_EQUIVS ) {
		STATUS("Not refining equivs model\n");
		return;
	}

	sdparams = error_model_params_vector(emodel);
	stepsize = error_model_step_vector(emodel);

	myfunc.n = error_model_num_params(emodel);
	myfunc.f = error_model_norm_res_func(emodel->type);
	myfunc.params = full;

	mini = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, myfunc.n);
	gsl_multimin_fminimizer_set(mini, &myfunc, sdparams, stepsize);
	gsl_vector_free(sdparams);
	gsl_vector_free(stepsize);

	STATUS("Refining error model...\n");
	niter = 0;
	do {
		niter++;
		r = gsl_multimin_fminimizer_iterate(mini);
		if ( r ) break;
		r = gsl_multimin_test_size(mini->size, 0.1);
		STATUS("%2i  |   ", niter);
		show_vector_oneline(mini->x);
	} while ( r == GSL_CONTINUE && niter < 20 );
	STATUS("Done.\n");

	error_model_params_set_from_vector(emodel, mini->x);

	gsl_multimin_fminimizer_free(mini);
}


ErrorModelType parse_error_model(const char *str, int *err)
{
	*err = 0;

	if ( strcmp(str, "equivs") == 0 ) {
		return EMODEL_EQUIVS;

	} else if ( strcmp(str, "ev11") == 0 ) {
		return EMODEL_EV11;

	} else if ( strcmp(str, "xscale") == 0 ) {
		return EMODEL_EV06;

	} else if ( strcmp(str, "ev06") == 0 ) {
		return EMODEL_EV06;

	} else if ( strcmp(str, "kh23") == 0 ) {
		return EMODEL_KH23;

	} else {
		*err = 1;
		return EMODEL_EQUIVS;
	}
}


static double isigi_asymptotic(ErrorModel *emodel)
{
	switch ( emodel->type ) {

		case EMODEL_EQUIVS:
		return NAN;

		case EMODEL_EV11:
		case EMODEL_KH23:
		case EMODEL_EV06:
		return 1.0/(emodel->sdfac*emodel->sdadd);
		/* The SdB term in Ev11 and Kh23 does not affect ISa */

	}
	abort();
}


void print_error_model(ErrorModel *emodel)
{
	switch ( emodel->type ) {

		case EMODEL_EQUIVS:
		STATUS("No error modelling was performed.");
		break;

		case EMODEL_EV11:
		STATUS("Error model parameters (Ev11): sdFac=%.3f, sdB=%.3f, sdAdd=%.3f\n",
		        emodel->sdfac, emodel->sdb, emodel->sdadd);
		STATUS("Overall (I/sigI)_asymptotic = %.3f\n", isigi_asymptotic(emodel));
		break;

		case EMODEL_KH23:
		STATUS("Error model parameters (Kh23): sdFac=%.3f, sdB=%.3f, sdAdd=%.3f\n",
		        emodel->sdfac, emodel->sdb, emodel->sdadd);
		STATUS("Overall (I/sigI)_asymptotic = %.3f\n", isigi_asymptotic(emodel));
		break;

		case EMODEL_EV06:
		STATUS("Error model parameters (Ev06): sdFac=%.3f, sdAdd=%.3f\n",
		        emodel->sdfac, emodel->sdadd);
		STATUS("Overall (I/sigI)_asymptotic = %.3f\n", isigi_asymptotic(emodel));
		break;

	}

}


void error_model_free(ErrorModel *emodel)
{
	free(emodel);
}
