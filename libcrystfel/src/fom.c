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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

enum fom get_fom(const char *s)
{
	if ( strcasecmp(s, "r1i") == 0 ) return FOM_R1I;
	if ( strcasecmp(s, "r1f") == 0 ) return FOM_R1F;
	if ( strcasecmp(s, "r2") == 0 ) return FOM_R2;
	if ( strcasecmp(s, "rsplit") == 0 ) return FOM_RSPLIT;
	if ( strcasecmp(s, "cc") == 0 ) return FOM_CC;
	if ( strcasecmp(s, "ccstar") == 0 ) return FOM_CCSTAR;
	if ( strcasecmp(s, "ccano") == 0 ) return FOM_CCANO;
	if ( strcasecmp(s, "crdano") == 0 ) return FOM_CRDANO;
	if ( strcasecmp(s, "rano") == 0 ) return FOM_RANO;
	if ( strcasecmp(s, "rano/rsplit") == 0 ) return FOM_RANORSPLIT;
	if ( strcasecmp(s, "d1sig") == 0 ) return FOM_D1SIG;
	if ( strcasecmp(s, "d2sig") == 0 ) return FOM_D2SIG;

	ERROR("Unknown figure of merit '%s'.\n", s);
	exit(1);
}


static struct fom_context *init_fom(enum fom fom, int nmax, int nshells)
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

	switch ( fctx->fom ) {

		case FOM_RANORSPLIT :
		fctx->num2 = malloc(nshells*sizeof(double));
		fctx->den2 = malloc(nshells*sizeof(double));
		if ( (fctx->num2 == NULL) || (fctx->den2 == NULL) ) return NULL;
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
		fctx->num = malloc(nshells*sizeof(double));
		fctx->den = malloc(nshells*sizeof(double));
		if ( (fctx->num == NULL) || (fctx->den == NULL) ) return NULL;
		for ( i=0; i<nshells; i++ ) {
			fctx->num[i] = 0.0;
			fctx->den[i] = 0.0;
		}
		break;

		case FOM_CC :
		case FOM_CCSTAR :
		case FOM_CCANO :
		case FOM_CRDANO :
		fctx->vec1 = malloc(nshells*sizeof(double *));
		fctx->vec2 = malloc(nshells*sizeof(double *));
		if ( (fctx->vec1 == NULL) || (fctx->vec2 == NULL) ) return NULL;
		for ( i=0; i<nshells; i++ ) {
			fctx->vec1[i] = malloc(nmax*sizeof(double));
			if ( fctx->vec1[i] == NULL ) return NULL;
			fctx->vec2[i] = malloc(nmax*sizeof(double));
			if ( fctx->vec2[i] == NULL ) return NULL;
			fctx->n = malloc(nshells*sizeof(int));
			if ( fctx->n == NULL ) return NULL;
		}
		for ( i=0; i<nshells; i++ ) {
			fctx->n[i] = 0;
		}
		fctx->nmax = nmax;
		break;

		case FOM_D1SIG :
		case FOM_D2SIG :
		fctx->n_within = malloc(nshells*sizeof(int));
		if ( fctx->n_within == NULL ) return NULL;
		for ( i=0; i<nshells; i++ ) {
			fctx->n_within[i] = 0;
		}
		break;

	}

	return fctx;
}


static void add_to_fom(struct fom_context *fctx, double i1, double i2,
                       double i1bij, double i2bij, double sig1, double sig2,
                       int bin)
{
	double f1, f2;
	double im, imbij;

	fctx->cts[bin]++;

	/* Negative intensities have already been weeded out. */
	f1 = sqrt(i1);
	f2 = sqrt(i2);

	switch ( fctx->fom ) {

		case FOM_R1I :
		fctx->num[bin] += fabs(i1 - i2);
		fctx->den[bin] += i1;
		break;

		case FOM_R1F :
		fctx->num[bin] += fabs(f1 - f2);
		fctx->den[bin] += f1;
		break;

		case FOM_R2 :
		fctx->num[bin] += pow(i1 - i2, 2.0);
		fctx->den[bin] += pow(i1, 2.0);
		break;

		case FOM_RSPLIT :
		fctx->num[bin] += fabs(i1 - i2);
		fctx->den[bin] += i1 + i2;
		break;

		case FOM_CC :
		case FOM_CCSTAR :
		assert(fctx->n[bin] < fctx->nmax);
		fctx->vec1[bin][fctx->n[bin]] = i1;
		fctx->vec2[bin][fctx->n[bin]] = i2;
		fctx->n[bin]++;
		break;

		case FOM_CCANO :
		case FOM_CRDANO :
		assert(fctx->n[bin] < fctx->nmax);
		fctx->vec1[bin][fctx->n[bin]] = i1 - i1bij;
		fctx->vec2[bin][fctx->n[bin]] = i2 - i2bij;
		fctx->n[bin]++;
		break;

		case FOM_RANORSPLIT :
		fctx->num2[bin] += fabs(i1 - i2);
		fctx->den2[bin] += i1 + i2;
		/* Intentional fall-through (no break) */

		case FOM_RANO :
		im = (i1 + i2)/2.0;
		imbij = (i1bij + i2bij)/2.0;
		fctx->num[bin] += fabs(im - imbij);
		fctx->den[bin] += im + imbij;
		break;

		case FOM_D1SIG :
		if ( fabs(i1-i2) < sqrt(sig1*sig1 + sig2*sig2) ) {
			fctx->n_within[bin]++;
		}
		break;

		case FOM_D2SIG :
		if ( fabs(i1-i2) < 2.0*sqrt(sig1*sig1 + sig2*sig2) ) {
			fctx->n_within[bin]++;
		}
		break;

	}
}


double fom_overall(struct fom_context *fctx)
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

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
		case FOM_RANO :
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

	}

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		return overall_num/overall_den;

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


double fom_shell(struct fom_context *fctx, int i)
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

	}

	ERROR("This point is never reached.\n");
	abort();
}


struct shells *make_intensity_shells(double min_I, double max_I, int nshells)
{
	struct shells *s;
	int i;

	if ( min_I >= max_I ) {
		ERROR("Invalid intensity range.\n");
		return NULL;
	}

	/* Adjust minimum and maximum intensities to get the most densely
	 * populated part of the reflections */
	max_I = min_I + (max_I-min_I)/5000.0;

	s = malloc(sizeof(struct shells));
	if ( s == NULL ) return NULL;

	s->rmins = malloc(nshells*sizeof(double));
	s->rmaxs = malloc(nshells*sizeof(double));

	if ( (s->rmins==NULL) || (s->rmaxs==NULL) ) {
		ERROR("Couldn't allocate memory for shells.\n");
		free(s);
		return NULL;
	}

	s->config_intshells = 1;
	s->nshells = nshells;

	for ( i=0; i<nshells; i++ ) {

		s->rmins[i] = min_I + i*(max_I - min_I)/nshells;;
		s->rmaxs[i] = min_I + (i+1)*(max_I - min_I)/nshells;;

	}

	return s;
}


struct shells *make_resolution_shells(double rmin, double rmax, int nshells)
{
	struct shells *s;
	double total_vol, vol_per_shell;
	int i;

	s = malloc(sizeof(struct shells));
	if ( s == NULL ) return NULL;

	s->rmins = malloc(nshells*sizeof(double));
	s->rmaxs = malloc(nshells*sizeof(double));

	if ( (s->rmins==NULL) || (s->rmaxs==NULL) ) {
		ERROR("Couldn't allocate memory for resolution shells.\n");
		free(s);
		return NULL;
	}

	s->config_intshells = 0;
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


double shell_label(struct shells *s, int i)
{
	if ( s->config_intshells ) {
		return (i+0.5) / s->nshells;
	} else {
		return s->rmins[i] + (s->rmaxs[i] - s->rmins[i])/2.0;
	}
}


static int get_bin(struct shells *s, Reflection *refl, UnitCell *cell)
{
	if ( s->config_intshells ) {

		double intensity;
		int bin, j;

		intensity = get_intensity(refl);

		bin = -1;
		for ( j=0; j<s->nshells; j++ ) {
			if ( (intensity>s->rmins[j])
			  && (intensity<=s->rmaxs[j]) )
			{
				bin = j;
				break;
			}
		}

		return bin;

	} else {

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



struct fom_context *fom_calculate(RefList *list1, RefList *list2, UnitCell *cell,
                                  struct shells *shells, enum fom fom,
                                  int noscale, SymOpList *sym)
{
	Reflection *refl1;
	RefListIterator *iter;
	struct fom_context *fctx;
	int n_out;

	fctx = init_fom(fom, num_reflections(list1), shells->nshells);

	if ( fctx == NULL ) {
		ERROR("Couldn't allocate memory for resolution shells.\n");
		return NULL;
	}

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

	n_out = 0;
	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		int bin;
		double i1, i2;
		double i1bij, i2bij;
		double sig1, sig2;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);

		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;

		bin = get_bin(shells, refl1, cell);
		if ( bin == -1 ) {
			n_out++;
			continue;
		}

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		sig1 = get_esd_intensity(refl1);
		sig2 = get_esd_intensity(refl2);

		if ( (fom == FOM_CCANO) || (fom == FOM_CRDANO)
		  || (fom == FOM_RANO) || (fom == FOM_RANORSPLIT) )
		{

			Reflection *refl1_bij = NULL;
			Reflection *refl2_bij = NULL;
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

			i1bij = get_intensity(refl1_bij);
			i2bij = get_intensity(refl2_bij);

		} else {

			/* Make it obvious if these get used by mistake */
			i1bij = +INFINITY;
			i2bij = +INFINITY;

		}

		add_to_fom(fctx, i1, i2, i1bij, i2bij, sig1, sig2, bin);
	}
	if ( n_out)  {
		ERROR("WARNING: %i reflection pairs outside range.\n", n_out);
	}

	return fctx;
}

