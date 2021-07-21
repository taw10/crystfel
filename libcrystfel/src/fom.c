/*
 * fom.c
 *
 * Figure of merit calculation
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
 *   2013      Lorenzo Galli <lorenzo.galli@desy.de>
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

#include <libcrystfel-config.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>
#include <assert.h>

#include "utils.h"
#include "fom.h"
#include "cell.h"
#include "cell-utils.h"
#include "reflist.h"
#include "reflist-utils.h"

/**
 * \file fom.h
 */

struct fom_context
{
	enum fom_type fom;
	int nshells;
	int *cts;

	/* For R-factors */
	double *num;
	double *den;

	/* For "double" R-factors */
	double *num2;
	double *den2;

	/* For CCs */
	double **vec1;
	double **vec2;
	int *n;
	int nmax;

	/* For "counting" things e.g. d1sig or d2sig */
	int *n_within;

	long int *n_meas;
	long int *possible;
};


static struct fom_context *init_fom(enum fom_type fom, int nmax, int nshells)
{
	struct fom_context *fctx;
	int i;

	fctx = malloc(sizeof(struct fom_context));
	if ( fctx == NULL ) return NULL;

	fctx->fom = fom;
	fctx->nshells = nshells;
	fctx->cts = malloc(nshells*sizeof(int));
	for ( i=0; i<nshells; i++ ) {
		fctx->cts[i] = 0;
	}

	fctx->num2 = NULL;
	fctx->den2 = NULL;
	fctx->num = NULL;
	fctx->den = NULL;
	fctx->n_meas = NULL;
	fctx->vec1 = NULL;
	fctx->vec2 = NULL;
	fctx->n = NULL;
	fctx->n_within = NULL;
	fctx->possible = NULL;

	switch ( fctx->fom ) {

		case FOM_RANORSPLIT :
		fctx->num2 = malloc(nshells*sizeof(double));
		fctx->den2 = malloc(nshells*sizeof(double));
		if ( (fctx->num2 == NULL) || (fctx->den2 == NULL) ) goto out;
		for ( i=0; i<nshells; i++ ) {
			fctx->num2[i] = 0.0;
			fctx->den2[i] = 0.0;
		}
		/* Intentional fall-through (no break) */

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
		case FOM_RANO :
		case FOM_MEAN_INTENSITY :
		case FOM_SNR :
		case FOM_REDUNDANCY :
		fctx->num = malloc(nshells*sizeof(double));
		fctx->den = malloc(nshells*sizeof(double));
		if ( (fctx->num == NULL) || (fctx->den == NULL) ) goto out;
		for ( i=0; i<nshells; i++ ) {
			fctx->num[i] = 0.0;
			fctx->den[i] = 0.0;
		}
		break;

		case FOM_COMPLETENESS :
		/* Uses 'cts' and 'possible' only - see calculate_possible() */
		break;

		case FOM_NUM_MEASUREMENTS :
		fctx->n_meas = calloc(nshells, sizeof(long int));
		if ( fctx->n_meas == NULL ) goto out;
		break;

		case FOM_CC :
		case FOM_CCSTAR :
		case FOM_CCANO :
		case FOM_CRDANO :
		fctx->vec1 = malloc(nshells*sizeof(double *));
		fctx->vec2 = malloc(nshells*sizeof(double *));
		if ( (fctx->vec1 == NULL) || (fctx->vec2 == NULL) ) goto out;
		for ( i=0; i<nshells; i++ ) {
			fctx->vec1[i] = NULL;
			fctx->vec2[i] = NULL;
		}
		for ( i=0; i<nshells; i++ ) {
			fctx->vec1[i] = malloc(nmax*sizeof(double));
			if ( fctx->vec1[i] == NULL ) goto out;
			fctx->vec2[i] = malloc(nmax*sizeof(double));
			if ( fctx->vec2[i] == NULL ) goto out;
		}
		fctx->n = malloc(nshells*sizeof(int));
		if ( fctx->n == NULL ) goto out;
		for ( i=0; i<nshells; i++ ) {
			fctx->n[i] = 0;
		}
		fctx->nmax = nmax;
		break;

		case FOM_D1SIG :
		case FOM_D2SIG :
		fctx->n_within = malloc(nshells*sizeof(int));
		if ( fctx->n_within == NULL ) goto out;
		for ( i=0; i<nshells; i++ ) {
			fctx->n_within[i] = 0;
		}
		break;

	}

	return fctx;

out:
	free(fctx->num2);
	free(fctx->den2);
	free(fctx->num);
	free(fctx->den);
	free(fctx->n_meas);
	if ( fctx->vec1 != NULL ) {
		for ( i=0; i<nshells; i++ ) {
			free(fctx->vec1[i]);
		}
		free(fctx->vec1);
	}
	if ( fctx->vec2 != NULL ) {
		for ( i=0; i<nshells; i++ ) {
			free(fctx->vec2[i]);
		}
		free(fctx->vec2);
	}
	free(fctx->n);
	free(fctx->n_within);
	free(fctx);
	return NULL;
}


static int add_to_fom(struct fom_context *fctx,
                      Reflection *refl1,
                      Reflection *refl2,
                      Reflection *refl1bij,
                      Reflection *refl2bij,
                      int bin)
{
	double i1, i2, i1bij, i2bij, sig1, sig2;
	double im, imbij;
	int bad = 0;

	fctx->cts[bin]++;

	switch ( fctx->fom ) {

		case FOM_R1I :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		fctx->num[bin] += fabs(i1 - i2);
		fctx->den[bin] += i1;
		break;

		case FOM_R1F :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		fctx->num[bin] += fabs(sqrt(i1) - sqrt(i2));
		fctx->den[bin] += sqrt(i1);
		break;

		case FOM_R2 :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		fctx->num[bin] += pow(i1 - i2, 2.0);
		fctx->den[bin] += pow(i1, 2.0);
		break;

		case FOM_RSPLIT :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		fctx->num[bin] += fabs(i1 - i2);
		fctx->den[bin] += i1 + i2;
		break;

		case FOM_CC :
		case FOM_CCSTAR :
		assert(fctx->n[bin] < fctx->nmax);
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		fctx->vec1[bin][fctx->n[bin]] = i1;
		fctx->vec2[bin][fctx->n[bin]] = i2;
		fctx->n[bin]++;
		break;

		case FOM_CCANO :
		case FOM_CRDANO :
		assert(fctx->n[bin] < fctx->nmax);
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		i1bij = get_intensity(refl1bij);
		i2bij = get_intensity(refl2bij);
		fctx->vec1[bin][fctx->n[bin]] = i1 - i1bij;
		fctx->vec2[bin][fctx->n[bin]] = i2 - i2bij;
		fctx->n[bin]++;
		break;

		case FOM_RANORSPLIT :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		fctx->num2[bin] += fabs(i1 - i2);
		fctx->den2[bin] += i1 + i2;
		/* Intentional fall-through (no break) */

		case FOM_RANO :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		i1bij = get_intensity(refl1bij);
		i2bij = get_intensity(refl2bij);
		im = (i1 + i2)/2.0;
		imbij = (i1bij + i2bij)/2.0;
		fctx->num[bin] += fabs(im - imbij);
		fctx->den[bin] += im + imbij;
		break;

		case FOM_D1SIG :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		sig1 = get_esd_intensity(refl1);
		sig2 = get_esd_intensity(refl2);
		if ( fabs(i1-i2) < sqrt(sig1*sig1 + sig2*sig2) ) {
			fctx->n_within[bin]++;
		}
		break;

		case FOM_D2SIG :
		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		sig1 = get_esd_intensity(refl1);
		sig2 = get_esd_intensity(refl2);
		if ( fabs(i1-i2) < 2.0*sqrt(sig1*sig1 + sig2*sig2) ) {
			fctx->n_within[bin]++;
		}
		break;

		case FOM_NUM_MEASUREMENTS :
		fctx->n_meas[bin] += get_redundancy(refl1);
		break;

		case FOM_REDUNDANCY :
		fctx->num[bin] += get_redundancy(refl1);
		fctx->den[bin] += 1.0;
		break;

		case FOM_SNR :
		i1 = get_intensity(refl1);
		sig1 = get_esd_intensity(refl1);
		if ( isfinite(i1/sig1) ) {
			fctx->num[bin] += i1/sig1;
			fctx->den[bin] += 1.0;
		} else {
			bad = 1;
		}
		break;

		case FOM_MEAN_INTENSITY :
		i1 = get_intensity(refl1);
		fctx->num[bin] += i1;
		fctx->den[bin] += 1.0;
		break;

		case FOM_COMPLETENESS :
		/* fctx->cts already incremented, as needed.
		 * Will calculate possible reflections later */
		break;

	}

	return bad;
}


/**
 * Calculates the overall value for the %fom_context
 *
 * You must have previously called fom_calculate()
 */
double fom_overall_value(struct fom_context *fctx)
{
	double overall_num = INFINITY;
	double overall_den = 0.0;
	double overall_num2 = INFINITY;
	double overall_den2 = 0.0;
	int i;
	double *overall_vec1;
	double *overall_vec2;
	int overall_n;
	double *overall_along_diagonal;
	double *overall_perpend_diagonal;
	double variance_signal;
	double variance_error;
	double cc = INFINITY;
	long int total_meas = 0;
	long int overall_cts = 0;
	long int overall_possible = 0;

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
		case FOM_RANO :
		case FOM_REDUNDANCY :
		case FOM_SNR :
		case FOM_MEAN_INTENSITY :
		overall_num = 0.0;
		overall_den = 0.0;
		for ( i=0; i<fctx->nshells; i++ ) {
			overall_num += fctx->num[i];
			overall_den += fctx->den[i];
		}
		break;

		case FOM_RANORSPLIT :
		overall_num = 0.0;
		overall_den = 0.0;
		for ( i=0; i<fctx->nshells; i++ ) {
			overall_num += fctx->num[i];
			overall_den += fctx->den[i];
		}
		overall_num2 = 0.0;
		overall_den2 = 0.0;
		for ( i=0; i<fctx->nshells; i++ ) {
			overall_num2 += fctx->num2[i];
			overall_den2 += fctx->den2[i];
		}
		break;

		case FOM_CC :
		case FOM_CCSTAR :
		case FOM_CCANO :
		overall_vec1 = malloc(fctx->nmax*sizeof(double));
		overall_vec2 = malloc(fctx->nmax*sizeof(double));
		overall_n = 0;
		for ( i=0; i<fctx->nshells; i++ ) {
			int j;
			for ( j=0; j<fctx->n[i]; j++ ) {
				overall_vec1[overall_n] = fctx->vec1[i][j];
				overall_vec2[overall_n] = fctx->vec2[i][j];
				overall_n++;
			}
		}
		cc = gsl_stats_correlation(overall_vec1, 1, overall_vec2, 1,
		                           overall_n);
		free(overall_vec1);
		free(overall_vec2);
		break;

		case FOM_CRDANO :
		overall_along_diagonal = malloc(fctx->nmax*sizeof(double));
		overall_perpend_diagonal = malloc(fctx->nmax*sizeof(double));
		overall_n = 0;
		for ( i=0; i<fctx->nshells; i++ ) {
			int j;
			for ( j=0; j<fctx->n[i]; j++ ) {
				overall_along_diagonal[overall_n] =
					 ( fctx->vec1[i][j] + fctx->vec2[i][j] )
					 / sqrt(2.0);
				overall_perpend_diagonal[overall_n] =
					 ( fctx->vec1[i][j] - fctx->vec2[i][j] )
					 / sqrt(2.0);
				overall_n++;
			}
		}
		variance_signal = gsl_stats_variance_m(overall_along_diagonal,
		                                       1, overall_n, 0.0);
		variance_error = gsl_stats_variance_m(overall_perpend_diagonal,
		                                      1, overall_n, 0.0);
		cc = sqrt(variance_signal / variance_error );

		free(overall_along_diagonal);
		free(overall_perpend_diagonal);
		break;

		case FOM_D1SIG :
		case FOM_D2SIG :
		overall_num = 0.0;
		overall_den = 0.0;
		for ( i=0; i<fctx->nshells; i++ ) {
			overall_num += fctx->n_within[i];
			overall_den += fctx->cts[i];
		}
		break;

		case FOM_NUM_MEASUREMENTS :
		total_meas = 0;
		for ( i=0; i<fctx->nshells; i++ ) {
			total_meas += fctx->n_meas[i];
		}
		break;

		case FOM_COMPLETENESS :
		for ( i=0; i<fctx->nshells; i++ ) {
			overall_cts += fctx->cts[i];
			overall_possible += fctx->possible[i];
		}
		break;

	}

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_REDUNDANCY :
		case FOM_SNR :
		case FOM_MEAN_INTENSITY :
		return overall_num/overall_den;

		case FOM_COMPLETENESS :
		return (double)overall_cts / overall_possible;

		case FOM_NUM_MEASUREMENTS :
		return total_meas;

		case FOM_R2 :
		return sqrt(overall_num/overall_den);

		case FOM_RSPLIT :
		return 2.0*(overall_num/overall_den) / sqrt(2.0);

		case FOM_CC :
		case FOM_CCANO :
		case FOM_CRDANO :
		return cc;

		case FOM_CCSTAR :
		return sqrt((2.0*cc)/(1.0+cc));

		case FOM_RANO :
		return 2.0*(overall_num/overall_den);

		case FOM_RANORSPLIT :
		return (2.0*(overall_num/overall_den)) /
		       (2.0*(overall_num2/overall_den2) / sqrt(2.0));

		case FOM_D1SIG :
		case FOM_D2SIG :
		return overall_num/overall_den;

	}

	ERROR("This point is never reached.\n");
	abort();
}


/**
 * Calculates the figure of merit for the specified shell number.
 * You must have previously called fom_calculate()
 */
double fom_shell_value(struct fom_context *fctx, int i)
{
	double cc;
	int j;
	double variance_signal;
	double variance_error;
	double *along_diagonal;
	double *perpend_diagonal;

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_REDUNDANCY :
		case FOM_SNR :
		case FOM_MEAN_INTENSITY :
		return fctx->num[i]/fctx->den[i];

		case FOM_R2 :
		return sqrt(fctx->num[i]/fctx->den[i]);

		case FOM_RSPLIT :
		return 2.0*(fctx->num[i]/fctx->den[i]) / sqrt(2.0);

		case FOM_CC :
		case FOM_CCANO :
		return gsl_stats_correlation(fctx->vec1[i], 1, fctx->vec2[i], 1,
		                             fctx->n[i]);

		case FOM_CCSTAR :
		cc = gsl_stats_correlation(fctx->vec1[i], 1, fctx->vec2[i], 1,
		                           fctx->n[i]);
		return sqrt((2.0*cc)/(1.0+cc));

		case FOM_RANO :
		return 2.0 * fctx->num[i]/fctx->den[i];

		case FOM_RANORSPLIT :
		return (2.0*fctx->num[i]/fctx->den[i]) /
		       (2.0*(fctx->num2[i]/fctx->den2[i]) / sqrt(2.0));

		case FOM_CRDANO :
		along_diagonal = malloc(fctx->n[i] * sizeof(double));
		perpend_diagonal = malloc(fctx->n[i] * sizeof(double));
		for ( j=0; j<fctx->n[i]; j++ ) {
			along_diagonal[j] = ( fctx->vec1[i][j] +
						 fctx->vec2[i][j] ) / sqrt(2.0);
			perpend_diagonal[j] = ( fctx->vec1[i][j] -
						 fctx->vec2[i][j] ) / sqrt(2.0);
		}
		variance_signal = gsl_stats_variance_m(along_diagonal, 1,
		                                       fctx->n[i], 0.0);
		variance_error = gsl_stats_variance_m(perpend_diagonal, 1,
		                                      fctx->n[i], 0.0);
		free(along_diagonal);
		free(perpend_diagonal);
		return sqrt(variance_signal / variance_error);

		case FOM_D1SIG :
		case FOM_D2SIG :
		return (double)fctx->n_within[i] / fctx->cts[i];

		case FOM_NUM_MEASUREMENTS :
		return fctx->n_meas[i];

		case FOM_COMPLETENESS :
		return (double)fctx->cts[i] / fctx->possible[i];

	}

	ERROR("This point is never reached.\n");
	abort();
}


/**
 * \param rmin: The minimum value of 1/d, in m^-1
 * \param rmax: The maximum value of 1/d, in m^-1
 * \param nshells: The number of shells to use
 *
 * Create a %fom_shells structure for the specified minimum and maximum
 * resolution limits
 *
 * Returns the %fom_shells structure, or NULL on error.
 */
struct fom_shells *fom_make_resolution_shells(double rmin, double rmax,
                                              int nshells)
{
	struct fom_shells *s;
	double total_vol, vol_per_shell;
	int i;

	s = malloc(sizeof(struct fom_shells));
	if ( s == NULL ) return NULL;

	s->rmins = malloc(nshells*sizeof(double));
	s->rmaxs = malloc(nshells*sizeof(double));

	if ( (s->rmins==NULL) || (s->rmaxs==NULL) ) {
		ERROR("Couldn't allocate memory for resolution shells.\n");
		free(s->rmins);
		free(s->rmaxs);
		free(s);
		return NULL;
	}

	s->nshells = nshells;

	total_vol = pow(rmax, 3.0) - pow(rmin, 3.0);
	vol_per_shell = total_vol / nshells;
	s->rmins[0] = rmin;
	for ( i=1; i<nshells; i++ ) {

		double r;

		r = vol_per_shell + pow(s->rmins[i-1], 3.0);
		r = pow(r, 1.0/3.0);

		/* Shells of constant volume */
		s->rmaxs[i-1] = r;
		s->rmins[i] = r;

	}
	s->rmaxs[nshells-1] = rmax;

	return s;
}


/**
 * \param s: A %fom_shells structure
 * \param i: The shell number
 *
 * Returns the value of 1/d at the middle of the shell,
 *  i.e. the mean of the minimum and maximum 1/d values for the shell
 */
double fom_shell_centre(struct fom_shells *s, int i)
{
	return s->rmins[i] + (s->rmaxs[i] - s->rmins[i])/2.0;
}


static int get_bin(struct fom_shells *s, Reflection *refl, UnitCell *cell)
{
	double d;
	int bin, j;
	signed int h, k, l;

	get_indices(refl, &h, &k, &l);
	d = 2.0 * resolution(cell, h, k, l);

	bin = -1;
	for ( j=0; j<s->nshells; j++ ) {
		if ( (d>s->rmins[j]) && (d<=s->rmaxs[j]) ) {
			bin = j;
			break;
		}
	}

	/* Allow for slight rounding errors */
	if ( (bin == -1) && (d <= s->rmins[0]) ) bin = 0;
	if ( (bin == -1) && (d >= s->rmaxs[s->nshells-1]) ) bin = 0;
	assert(bin != -1);

	return bin;
}


static int wilson_scale(RefList *list1, RefList *list2, UnitCell *cell)
{
	Reflection *refl1;
	Reflection *refl2;
	RefListIterator *iter;
	int max_n = 256;
	int n = 0;
	double *x;
	double *y;
	int r;
	double G, B;
	double c0, c1, cov00, cov01, cov11, chisq;

	x = malloc(max_n*sizeof(double));
	y = malloc(max_n*sizeof(double));
	if ( (x==NULL) || (y==NULL) ) {
		ERROR("Failed to allocate memory for scaling.\n");
		return 1;
	}

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		double Ih1, Ih2;
		double res;

		get_indices(refl1, &h, &k, &l);
		res = resolution(cell, h, k, l);

		refl2 = find_refl(list2, h, k, l);
		assert(refl2 != NULL);

		Ih1 = get_intensity(refl1);
		Ih2 = get_intensity(refl2);

		if ( (Ih1 <= 0.0) || (Ih2 <= 0.0) ) continue;
		if ( isnan(Ih1) || isinf(Ih1) ) continue;
		if ( isnan(Ih2) || isinf(Ih2) ) continue;

		if ( n == max_n ) {
			max_n *= 2;
			x = realloc(x, max_n*sizeof(double));
			y = realloc(y, max_n*sizeof(double));
			if ( (x==NULL) || (y==NULL) ) {
				ERROR("Failed to allocate memory for scaling.\n");
				return 1;
			}
		}

		x[n] = res*res;
		y[n] = log(Ih1/Ih2);
		n++;

	}

	if ( n < 2 ) {
		ERROR("Not enough reflections for scaling\n");
		return 1;
	}

	r = gsl_fit_linear(x, 1, y, 1, n, &c0, &c1,
	                     &cov00, &cov01, &cov11, &chisq);

	if ( r ) {
		ERROR("Scaling failed.\n");
		return 1;
	}

	G = exp(c0);
	B = c1/2.0;

	STATUS("Relative scale factor = %f, relative B factor = %f A^2\n",
	       G, B*1e20);
	STATUS("A scale factor greater than 1 means that the second reflection "
	       "list is weaker than the first.\n");
	STATUS("A positive relative B factor means that the second reflection "
	       "list falls off with resolution more quickly than the first.\n");

	free(x);
	free(y);

	/* Apply the scaling factor */
	for ( refl2 = first_refl(list2, &iter);
	      refl2 != NULL;
	      refl2 = next_refl(refl2, iter) )
	{
		signed int h, k, l;
		double res;
		double corr;

		get_indices(refl2, &h, &k, &l);
		res = resolution(cell, h, k, l);

		corr = G * exp(2.0*B*res*res);
		set_intensity(refl2, get_intensity(refl2)*corr);
		set_esd_intensity(refl2, get_esd_intensity(refl2)*corr);

	}
	return 0;
}


static int calculate_possible(struct fom_context *fctx,
                              struct fom_shells *shells,
                              UnitCell *cell,
                              const SymOpList *sym)
{
	RefList *counted;
	int hmax, kmax, lmax;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	signed int h, k, l;

	fctx->possible = calloc(fctx->nshells, sizeof(long int));
	if ( fctx->possible == NULL ) return 1;

	counted = reflist_new();
	if ( counted == NULL ) {
		free(fctx->possible);
		return 1;
	}

	cell_get_cartesian(cell, &ax, &ay, &az,
	                         &bx, &by, &bz,
	                         &cx, &cy, &cz);
	hmax = shells->rmaxs[fctx->nshells-1] * modulus(ax, ay, az);
	kmax = shells->rmaxs[fctx->nshells-1] * modulus(bx, by, bz);
	lmax = shells->rmaxs[fctx->nshells-1] * modulus(cx, cy, cz);
	for ( h=-hmax; h<=hmax; h++ ) {
	for ( k=-kmax; k<=kmax; k++ ) {
	for ( l=-lmax; l<=lmax; l++ ) {

		double d;
		signed int hs, ks, ls;
		int bin;
		int i;

		get_asymm(sym, h, k, l, &hs, &ks, &ls);
		d = 2.0 * resolution(cell, hs, ks, ls);

		if ( forbidden_reflection(cell, h, k, l) ) continue;

		bin = -1;
		for ( i=0; i<fctx->nshells; i++ ) {
			if ( (d>shells->rmins[i]) && (d<=shells->rmaxs[i]) ) {
				bin = i;
				break;
			}
		}
		if ( bin == -1 ) continue;

		if ( find_refl(counted, hs, ks, ls) != NULL ) continue;
		add_refl(counted, hs, ks, ls);

		fctx->possible[bin]++;

	}
	}
	}
	reflist_free(counted);

	return 0;
}


int fom_is_anomalous(enum fom_type fom)
{
	switch ( fom ) {

		case FOM_CCANO:
		case FOM_RANO:
		case FOM_CRDANO:
		case FOM_RANORSPLIT:
		return 1;

		case FOM_R1I:
		case FOM_R1F:
		case FOM_R2:
		case FOM_RSPLIT:
		case FOM_CC:
		case FOM_CCSTAR:
		case FOM_D1SIG:
		case FOM_D2SIG:
		case FOM_NUM_MEASUREMENTS:
		case FOM_REDUNDANCY:
		case FOM_SNR:
		case FOM_MEAN_INTENSITY:
		case FOM_COMPLETENESS:
		return 0;
	}

	ERROR("This point never reached\n");
	abort();
}


int fom_is_comparison(enum fom_type fom)
{
	switch ( fom ) {

		case FOM_CCANO:
		case FOM_RANO:
		case FOM_CRDANO:
		case FOM_RANORSPLIT:
		case FOM_R1I:
		case FOM_R1F:
		case FOM_R2:
		case FOM_RSPLIT:
		case FOM_CC:
		case FOM_CCSTAR:
		case FOM_D1SIG:
		case FOM_D2SIG:
		return 1;

		case FOM_NUM_MEASUREMENTS:
		case FOM_REDUNDANCY:
		case FOM_SNR:
		case FOM_MEAN_INTENSITY:
		case FOM_COMPLETENESS:
		return 0;
	}

	ERROR("This point never reached\n");
	abort();
}


static int is_single_list(enum fom_type fom)
{
	switch ( fom ) {

		case FOM_CCANO:
		case FOM_RANO:
		case FOM_CRDANO:
		case FOM_RANORSPLIT:
		case FOM_R1I:
		case FOM_R1F:
		case FOM_R2:
		case FOM_RSPLIT:
		case FOM_CC:
		case FOM_CCSTAR:
		case FOM_D1SIG:
		case FOM_D2SIG:
		return 0;

		case FOM_NUM_MEASUREMENTS:
		case FOM_REDUNDANCY:
		case FOM_SNR:
		case FOM_MEAN_INTENSITY:
		case FOM_COMPLETENESS:
		return 1;
	}

	ERROR("This point never reached\n");
	abort();
}


/**
 * \param list1: A %RefList
 * \param list2: A %RefList
 * \param cell: A %UnitCell
 * \param shells: A %fom_shells structure
 * \param fom: The figure of merit to calculate
 * \param noscale: Non-zero to disable scaline of reflection lists
 * \param sym: The symmetry of \p list1 and \p list2.
 *
 * Calculates the specified figure of merit, comparing the two reflection lists.
 *
 * The \p cell and \p sym must match both reflection lists.  You should also have
 * called fom_select_reflection_pairs() to pre-process the lists.
 *
 * If the figure of merit does not involve comparison (e.g. %FOM_SNR),
 * then \p list1 will be used.  In this case, \p list2 and \p noscale will be
 * ignored.  Use fom_select_reflections() instead of fom_select_reflection_pairs()
 * in this case.
 *
 * \returns a %fom_context structure.  Use fom_shell_value() et al., to
 *  extract the actual figure of merit values.
 */
struct fom_context *fom_calculate(RefList *list1, RefList *list2, UnitCell *cell,
                                  struct fom_shells *shells, enum fom_type fom,
                                  int noscale, const SymOpList *sym)
{
	Reflection *refl1;
	RefListIterator *iter;
	struct fom_context *fctx;
	long int n_out = 0;
	long int n_rej = 0;

	fctx = init_fom(fom, num_reflections(list1), shells->nshells);

	if ( fctx == NULL ) {
		ERROR("Couldn't allocate memory for resolution shells.\n");
		return NULL;
	}

	if ( !is_single_list(fom) ) {
		if ( !noscale && wilson_scale(list1, list2, cell) ) {
			ERROR("Error with scaling.\n");
			return NULL;
		}

		for ( refl1 = first_refl(list1, &iter);
		      refl1 != NULL;
		      refl1 = next_refl(refl1, iter) )
		{
			Reflection *refl2;
			signed int h, k, l;
			set_flag(refl1, 0);
			get_indices(refl1, &h, &k, &l);
			refl2 = find_refl(list2, h, k, l);
			assert(refl2 != NULL);
			set_flag(refl2, 0);
		}
	}

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		int bin;
		Reflection *refl2;
		Reflection *refl1_bij = NULL;
		Reflection *refl2_bij = NULL;

		get_indices(refl1, &h, &k, &l);

		if ( is_single_list(fom) ) {
			refl2 = NULL;
		} else {
			refl2 = find_refl(list2, h, k, l);
			if ( refl2 == NULL ) continue;
		}

		bin = get_bin(shells, refl1, cell);
		if ( bin == -1 ) {
			n_out++;
			continue;
		}

		if ( fom_is_anomalous(fom) ) {

			signed int hb, kb, lb;

			if ( find_equiv_in_list(list1, -h, -k, -l, sym,
			                        &hb, &kb, &lb) )
			{
				refl1_bij = find_refl(list1, hb, kb, lb);
			}

			if ( find_equiv_in_list(list2, -h, -k, -l, sym,
			                        &hb, &kb, &lb) )
			{
				refl2_bij = find_refl(list2, hb, kb, lb);
			}

			/* Each reflection must only be counted once, whether
			 * we are visiting it now as "normal" or "bij" */
			if ( get_flag(refl1) ) continue;
			assert(!get_flag(refl2));
			set_flag(refl1, 1);
			set_flag(refl1_bij, 1);
			set_flag(refl2, 1);
			set_flag(refl2_bij, 1);

			assert(refl1_bij != NULL);
			assert(refl2_bij != NULL);

		}

		n_rej += add_to_fom(fctx, refl1, refl2, refl1_bij, refl2_bij, bin);

	}
	if ( n_out )  {
		ERROR("WARNING: %i reflection pairs outside range.\n", n_out);
	}
	if ( n_rej ) {
		if ( fom == FOM_SNR ) {
			ERROR("WARNING: %li reflections had infinite or "
			      "invalid values of I/sigma(I).\n", n_rej);
		} else {
			ERROR("WARNING: %li reflections rejected by add_to_fom\n",
			      n_rej);
		}
	}

	if ( fom == FOM_COMPLETENESS ) {
		calculate_possible(fctx, shells, cell, sym);
	}

	return fctx;
}


/**
 * \param list1: The first input %RefList
 * \param list2: The second input %RefList
 * \param plist1_acc: Pointer to location for accepted list
 * \param plist2_acc: Pointer to location for accepted list
 * \param cell: A %UnitCell
 * \param sym: The symmetry of \p raw_list
 * \param anom: Non-zero if you will calculate a FoM for anomalous signal
 * \param rmin_fix: If positive, minimum resolution to use
 * \param rmax_fix: If positive, maximum resolution to use
 * \param sigma_cutoff: Minimum I/sigI value
 * \param ignore_negs: Non-zero to filter out negative intensities
 * \param zero_negs: Non-zero to set negative intensities to zero
 * \param mul_cutoff: Minimum number of measurements per reflection
 *
 * Selects reflections suitable for use with fom_calculate().
 *
 * Use -INFINITY for \p sigma_cutoff to disable the check.
 * Set \p mul_cutoff to zero to disable the check.
 *
 * \returns a %fom_rejections structure with the counts of reflections.
 */
struct fom_rejections fom_select_reflection_pairs(RefList *list1, RefList *list2,
                                                  RefList **plist1_acc,
                                                  RefList **plist2_acc,
                                                  UnitCell *cell, SymOpList *sym,
                                                  int anom, double rmin_fix, double rmax_fix,
                                                  double sigma_cutoff, int ignore_negs,
                                                  int zero_negs, int mul_cutoff)
{
	Reflection *refl1;
	RefListIterator *iter;
	struct fom_rejections rej;
	RefList *list1_acc;
	RefList *list2_acc;

	rej.common = 0;
	rej.low_snr = 0;
	rej.negative_deleted = 0;
	rej.negative_zeroed = 0;
	rej.few_measurements = 0;
	rej.outside_resolution_range = 0;
	rej.no_bijvoet = 0;
	rej.centric = 0;
	rej.nan_inf_value = 0;

	list1_acc = reflist_new();
	list2_acc = reflist_new();

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		double val1, val2;
		double esd1, esd2;
		int mul1, mul2;
		Reflection *refl2;
		Reflection *refl1_acc;
		Reflection *refl2_acc;

		get_indices(refl1, &h, &k, &l);

		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;

		val1 = get_intensity(refl1);
		val2 = get_intensity(refl2);

		esd1 = get_esd_intensity(refl1);
		esd2 = get_esd_intensity(refl2);

		mul1 = get_redundancy(refl1);
		mul2 = get_redundancy(refl2);

		if ( !isfinite(val1) || !isfinite(val2)
		  || !isfinite(esd1) || !isfinite(esd2) )
		{
			rej.nan_inf_value++;
			continue;
		}

		if ( (val1 < sigma_cutoff * esd1)
		  || (val2 < sigma_cutoff * esd2) )
		{
			rej.low_snr++;
			continue;
		}

		if ( ignore_negs && ((val1 < 0.0) || (val2 < 0.0)) ) {
			rej.negative_deleted++;
			continue;
		}

		if ( (mul1 < mul_cutoff) || (mul2 < mul_cutoff) ) {
			rej.few_measurements++;
			continue;
		}

		if ( zero_negs ) {
			int d = 0;
			if ( val1 < 0.0 ) {
				val1 = 0.0;
				d = 1;
			}
			if ( val2 < 0.0 ) {
				val2 = 0.0;
				d = 1;
			}
			if ( d ) rej.negative_zeroed++;
			continue;
		}

		if ( rmin_fix > 0.0 ) {
			double res = 2.0*resolution(cell, h, k, l);
			if ( res < rmin_fix ) {
				rej.outside_resolution_range++;
				continue;
			}
		}

		if ( rmax_fix > 0.0 ) {
			double res = 2.0*resolution(cell, h, k, l);
			if ( res > rmax_fix ) {
				rej.outside_resolution_range++;
				continue;
			}
		}

		refl1_acc = add_refl(list1_acc, h, k, l);
		copy_data(refl1_acc, refl1);
		set_intensity(refl1_acc, val1);

		refl2_acc = add_refl(list2_acc, h, k, l);
		copy_data(refl2_acc, refl2);
		set_intensity(refl2_acc, val2);

		rej.common++;

	}

	/* For anomalous figures of merit, we additionally require that we have
	 * all the Bijvoet pairs after the above rejection tests */
	if ( anom ) {

		list1 = list1_acc;
		list2 = list2_acc;
		list1_acc = reflist_new();
		list2_acc = reflist_new();

		rej.common = 0;

		for ( refl1 = first_refl(list1, &iter);
		      refl1 != NULL;
		      refl1 = next_refl(refl1, iter) )
		{
			Reflection *refl1_bij = NULL;
			Reflection *refl2_bij = NULL;
			signed int h, k, l;
			signed int hb, kb, lb;
			Reflection *refl1_acc;
			Reflection *refl2_acc;
			Reflection *refl2;
			double val1, val2;

			get_indices(refl1, &h, &k, &l);

			refl2 = find_refl(list2, h, k, l);
			assert(refl2 != NULL);

			val1 = get_intensity(refl1);
			val2 = get_intensity(refl2);

			if ( is_centric(h, k, l, sym) ) {
				rej.centric++;
				continue;
			}

			if ( find_equiv_in_list(list1, -h, -k, -l, sym,
			                        &hb, &kb, &lb) )
			{
				refl1_bij = find_refl(list1, hb, kb, lb);
			}

			if ( find_equiv_in_list(list2, -h, -k, -l, sym,
			                        &hb, &kb, &lb) )
			{
				refl2_bij = find_refl(list2, hb, kb, lb);
			}

			if ( (refl1_bij == NULL) || (refl2_bij == NULL) ) {
				rej.no_bijvoet++;
				continue;
			}

			refl1_acc = add_refl(list1_acc, h, k, l);
			copy_data(refl1_acc, refl1);
			set_intensity(refl1_acc, val1);

			refl2_acc = add_refl(list2_acc, h, k, l);
			copy_data(refl2_acc, refl2);
			set_intensity(refl2_acc, val2);

			rej.common++;
		}
	}

	*plist1_acc = list1_acc;
	*plist2_acc = list2_acc;
	return rej;
}


/**
 * \param raw_list: The input %RefList
 * \param plist_acc: Pointer to location for accepted list
 * \param cell: A %UnitCell
 * \param sym: The symmetry of \p raw_list
 * \param rmin_fix: If positive, minimum resolution to use
 * \param rmax_fix: If positive, maximum resolution to use
 * \param sigma_cutoff: Minimum I/sigI value
 * \param ignore_negs: Non-zero to filter out negative intensities
 * \param zero_negs: Non-zero to set negative intensities to zero
 * \param mul_cutoff: Minimum number of measurements per reflection
 *
 * Use -INFINITY for \p sigma_cutoff to disable the check.
 * Set \p mul_cutoff to zero to disable the check.
 *
 * \returns a %fom_rejections structure with the counts of reflections.
 */
struct fom_rejections fom_select_reflections(RefList *raw_list,
                                             RefList **plist_acc,
                                             UnitCell *cell, SymOpList *sym,
                                             double rmin_fix, double rmax_fix,
                                             double sigma_cutoff, int ignore_negs,
                                             int zero_negs, int mul_cutoff)
{
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	struct fom_rejections rej;

	*plist_acc = NULL;

	rej.common = 0;
	rej.low_snr = 0;
	rej.negative_deleted = 0;
	rej.negative_zeroed = 0;
	rej.few_measurements = 0;
	rej.outside_resolution_range = 0;
	rej.no_bijvoet = 0;
	rej.centric = 0;
	rej.nan_inf_value = 0;

	list = reflist_new();
	if ( list == NULL ) return rej;

	for ( refl = first_refl(raw_list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double val, sig;
		int ig = 0;
		Reflection *new;

		get_indices(refl, &h, &k, &l);

		val = get_intensity(refl);
		sig = get_esd_intensity(refl);

		if ( !isfinite(val) || !isfinite(sig) ) {
			rej.nan_inf_value++;
			continue;
		}

		if ( val < sigma_cutoff * sig ) {
			rej.low_snr++;
			ig = 1;
		}

		if ( ignore_negs && (val < 0.0) ) {
			rej.negative_deleted++;
			ig = 1;
		}

		if ( zero_negs && (val < 0.0) ) {
			set_intensity(refl, 0.0);
			rej.negative_zeroed++;
		}

		if ( rmin_fix > 0.0 ) {
			double res = 2.0*resolution(cell, h, k, l);
			if ( res < rmin_fix ) {
				rej.outside_resolution_range++;
				continue;
			}
		}

		if ( rmax_fix > 0.0 ) {
			double res = 2.0*resolution(cell, h, k, l);
			if ( res > rmax_fix ) {
				rej.outside_resolution_range++;
				continue;
			}
		}

		if ( ig ) continue;

		new = add_refl(list, h, k, l);
		copy_data(new, refl);
	}

	*plist_acc = list;
	return rej;
}


/**
 * \param fctx: A %fom_context structure
 *
 * \returns the total number of unique reflections
 */
int fom_overall_num_reflections(struct fom_context *fctx)
{
	int i;
	long int n = 0;

	for ( i=0; i<fctx->nshells; i++ ) {
		n += fctx->cts[i];
	}
	return n;
}


/**
 * \param fctx: A %fom_context structure
 * \param i: Shell number
 *
 * \returns the number of unique reflections in the shell
 */
int fom_shell_num_reflections(struct fom_context *fctx, int i)
{
	return fctx->cts[i];
}


/**
 * \param fctx: A %fom_context structure
 *
 * This must only be called on a %fom_context for %FOM_COMPLETENESS.
 *
 * \returns the total number of reflections possible in all shells, taking into
 * account symmetry and lattice absences, but not screw axis/glide place absences.
 */
int fom_overall_num_possible(struct fom_context *fctx)
{
	int i;
	long int n = 0;

	assert(fctx->fom == FOM_COMPLETENESS);

	for ( i=0; i<fctx->nshells; i++ ) {
		n += fctx->possible[i];
	}
	return n;
}


/**
 * \param fctx: A %fom_context structure
 * \param i: Shell number
 *
 * This must only be called on a %fom_context for %FOM_COMPLETENESS.
 *
 * \returns the number of reflections possible in the shell, taking into account
 * symmetry and lattice absences, but not screw axis/glide place absences.
 */
int fom_shell_num_possible(struct fom_context *fctx, int i)
{
	assert(fctx->fom == FOM_COMPLETENESS);
	return fctx->possible[i];
}


const char *fom_name(enum fom_type f)
{
	switch ( f ) {
		case FOM_R1I : return "R1(I)";
		case FOM_R1F : return "R1(F)";
		case FOM_R2 : return "R2";
		case FOM_RSPLIT : return "Rsplit";
		case FOM_CC : return "CC";
		case FOM_CCSTAR : return "CC*";
		case FOM_CCANO : return "CCano";
		case FOM_CRDANO : return "CRDano";
		case FOM_RANO : return "Rano";
		case FOM_RANORSPLIT : return "Rano/Rsplit";
		case FOM_D1SIG : return "D<1sigma";
		case FOM_D2SIG : return "D<2sigma";
		case FOM_NUM_MEASUREMENTS : return "nMeas";
		case FOM_REDUNDANCY : return "Redundancy";
		case FOM_SNR : return "I/sigI";
		case FOM_MEAN_INTENSITY : return "mean I";
		case FOM_COMPLETENESS : return "Completeness";
		default : return "unknown FoM";
	}
}
