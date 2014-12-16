/*
 * compare_hkl.c
 *
 * Compare reflection lists
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>

#include "version.h"
#include "utils.h"
#include "statistics.h"
#include "symmetry.h"
#include "reflist-utils.h"
#include "cell-utils.h"

enum fom
{
	FOM_R1I,
	FOM_R1F,
	FOM_R2,
	FOM_RSPLIT,
	FOM_CC,
	FOM_CCSTAR,
	FOM_CCANO,
	FOM_CRDANO,
	FOM_RANO,
	FOM_RANORSPLIT
};


static enum fom get_fom(const char *s)
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

	ERROR("Unknown figure of merit '%s'.\n", s);
	exit(1);
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file1.hkl> <file2.hkl>\n\n", s);
	printf(
"Compare intensity lists.\n"
"\n"
"  -y, --symmetry=<sym>       The symmetry of both the input files.\n"
"  -p, --pdb=<filename>       Unit cell file to use.\n"
"      --fom=<FoM>            Calculate this figure of merit  Choose from:\n"
"                              R1I, R1F, R2, Rsplit, CC, CCstar,\n"
"			       CCano, CRDano, Rano, Rano/Rsplit\n"
"      --nshells=<n>          Use <n> resolution shells.\n"
"  -u                         Force scale factor to 1.\n"
"      --shell-file=<file>    Write resolution shells to <file>.\n"
"\n"
"You can control which reflections are included in the calculation:\n"
"\n"
"      --ignore-negs          Ignore reflections with negative intensities.\n"
"      --zero-negs            Set negative intensities to zero.\n"
"      --sigma-cutoff=<n>     Discard reflections with I/sigma(I) < n.\n"
"      --rmin=<res>           Low resolution cutoff (1/d in m^-1).\n"
"      --rmax=<res>           High resolution cutoff (1/d in m^-1).\n"
"      --lowres=<n>           Low resolution cutoff in (d in A).\n"
"      --highres=<n>          High resolution cutoff in (d in A).\n"
"      --intensity-shells     Use shells of intensity instead of resolution.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
);
}


struct fom_context
{
	enum fom fom;
	int nshells;

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
};


static struct fom_context *init_fom(enum fom fom, int nmax, int nshells)
{
	struct fom_context *fctx;
	int i;

	fctx = malloc(sizeof(struct fom_context));
	if ( fctx == NULL ) return NULL;

	fctx->fom = fom;
	fctx->nshells = nshells;

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

	}

	return fctx;
}


static void add_to_fom(struct fom_context *fctx, double i1, double i2,
                       double i1bij, double i2bij, int bin)
{
	double f1, f2;
	double im, imbij;

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


	}
}


static double fom_overall(struct fom_context *fctx)
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

	}

	ERROR("This point is never reached.\n");
	abort();
}


static double fom_shell(struct fom_context *fctx, int i)
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

	}

	ERROR("This point is never reached.\n");
	abort();
}


struct shells
{
	int config_intshells;
	int nshells;
	double *rmins;
	double *rmaxs;
};


static struct shells *set_intensity_shells(double min_I, double max_I,
                                           int nshells)
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


static struct shells *set_resolution_shells(double rmin, double rmax,
                                            int nshells)
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


static double shell_label(struct shells *s, int i)
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


static void do_fom(RefList *list1, RefList *list2, UnitCell *cell,
                   double rmin, double rmax, enum fom fom,
                   int config_unity, int nshells, const char *filename,
                   int config_intshells, double min_I, double max_I,
                   SymOpList *sym)
{
	int *cts;
	int i;
	Reflection *refl1;
	RefListIterator *iter;
	FILE *fh;
	double scale;
	struct fom_context *fctx;
	struct shells *shells;
	const char *t1, *t2;
	int n_out;

	cts = malloc(nshells*sizeof(int));
	fctx = init_fom(fom, num_reflections(list1), nshells);

	if ( (fctx==NULL) || (cts==NULL) ) {
		ERROR("Couldn't allocate memory for resolution shells.\n");
		return;
	}

	for ( i=0; i<nshells; i++ ) {
		cts[i] = 0;
	}

	if ( config_unity ) {
		scale = 1.0;
	} else {
		scale = stat_scale_intensity(list1, list2);
	}

	/* Calculate the bins */
	if ( config_intshells ) {
		shells = set_intensity_shells(min_I, max_I, nshells);
	} else {
		shells = set_resolution_shells(rmin, rmax, nshells);
	}

	if ( shells == NULL ) {
		ERROR("Failed to set up shells.\n");
		return;
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
		i2 = scale * get_intensity(refl2);

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

			assert(refl1_bij != NULL);
			assert(refl2_bij != NULL);

			i1bij = get_intensity(refl1_bij);
			i2bij = scale * get_intensity(refl2_bij);

		} else {

			/* Make it obvious if these get used by mistake */
			i1bij = +INFINITY;
			i2bij = +INFINITY;

		}

		add_to_fom(fctx, i1, i2, i1bij, i2bij, bin);
		cts[bin]++;
	}
	if ( n_out)  {
		ERROR("WARNING: %i reflection pairs outside range.\n", n_out);
	}

	switch ( fom ) {

		case FOM_R1I :
		STATUS("Overall R1(I) = %.2f %%\n", 100.0*fom_overall(fctx));
		break;

		case FOM_R1F :
		STATUS("Overall R1(F) = %.2f %%\n", 100.0*fom_overall(fctx));
		break;

		case FOM_R2 :
		STATUS("Overall R(2) = %.2f %%\n", 100.0*fom_overall(fctx));
		break;

		case FOM_RSPLIT :
		STATUS("Overall Rsplit = %.2f %%\n", 100.0*fom_overall(fctx));
		break;

		case FOM_CC :
		STATUS("Overall CC = %.7f\n", fom_overall(fctx));
		break;

		case FOM_CCSTAR :
		STATUS("Overall CC* = %.7f\n", fom_overall(fctx));
		break;

		case FOM_CCANO :
		STATUS("Overall CCano = %.7f\n", fom_overall(fctx));
		break;

		case FOM_CRDANO :
		STATUS("Overall CRDano = %.7f\n", fom_overall(fctx));
		break;

		case FOM_RANO :
		STATUS("Overall Rano =  %.2f %%\n", 100.0*fom_overall(fctx));
		break;

		case FOM_RANORSPLIT :
		STATUS("Overall Rano/Rsplit =  %.7f\n", fom_overall(fctx));
		break;

	}

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		return;
	}

	if ( config_intshells ) {
		t1 = "Relative I  ";
		t2 = "";
	} else {
		t1 = "  1/d centre";
		t2 = "      d / A   Min 1/nm    Max 1/nm";
	}

	switch ( fom ) {

		case FOM_R1I :
		fprintf(fh, "%s  R1(I)/%%       nref%s\n", t1, t2);
		break;

		case FOM_R1F :
		fprintf(fh, "%s  R1(F)/%%       nref%s\n", t1, t2);
		break;

		case FOM_R2 :
		fprintf(fh, "%s     R2/%%       nref%s\n", t1, t2);
		break;

		case FOM_RSPLIT :
		fprintf(fh, "%s Rsplit/%%       nref%s\n", t1, t2);
		break;

		case FOM_CC :
		fprintf(fh, "%s       CC       nref%s\n", t1, t2);
		break;

		case FOM_CCSTAR :
		fprintf(fh, "%s      CC*       nref%s\n", t1, t2);
		break;

		case FOM_CCANO :
		fprintf(fh, "%s    CCano       nref%s\n", t1, t2);
		break;

		case FOM_CRDANO :
		fprintf(fh, "%s    CRDano       nref%s\n", t1, t2);
		break;

		case FOM_RANO :
		fprintf(fh, "%s   Rano/%%       nref%s\n", t1, t2);
		break;

		case FOM_RANORSPLIT :
		fprintf(fh, "%s Rano/Rsplit       nref%s\n", t1, t2);
		break;

	}

	for ( i=0; i<nshells; i++ ) {

		double r, cen;

		cen = shell_label(shells, i);
		r = fom_shell(fctx, i);

		switch ( fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
		case FOM_RANO :
		if ( config_intshells ) {
			fprintf(fh, "%10.3f %10.2f %10i\n",
				cen, r*100.0, cts[i]);
		} else {
			fprintf(fh, "%10.3f %10.2f %10i %10.2f "
			            "%10.3f  %10.3f\n",
				cen*1.0e-9, r*100.0, cts[i], (1.0/cen)*1e10,
				shells->rmins[i]*1.0e-9,
				shells->rmaxs[i]*1.0e-9);
		}
		break;

		case FOM_CC :
		case FOM_CCSTAR :
		case FOM_CCANO :
		case FOM_CRDANO :
		if ( config_intshells ) {
			fprintf(fh, "%10.3f %10.7f %10i\n",
				cen, r, cts[i]);
		} else {
			fprintf(fh, "%10.3f %10.7f %10i %10.2f "
			            "%10.3f  %10.3f\n",
				cen*1.0e-9, r, cts[i], (1.0/cen)*1e10,
				shells->rmins[i]*1.0e-9,
				shells->rmaxs[i]*1.0e-9);
		}
		break;

		case FOM_RANORSPLIT :
		if ( config_intshells ) {
			fprintf(fh, "%10.3f    %10.7f %10i\n",
				cen, r, cts[i]);
		} else {
			fprintf(fh, "%10.3f    %10.7f %10i %10.2f "
			            "%10.3f  %10.3f\n",
				cen*1.0e-9, r, cts[i], (1.0/cen)*1e10,
				shells->rmins[i]*1.0e-9, shells->rmaxs[i]*1.0e-9);
		}
		break;


	}

	}

	fclose(fh);
}


static void check_highres()
{
	static int have = 0;
	if ( have ) {
		ERROR("You cannot use --rmax and --highres at the same time.\n");
		exit(1);
	}
	have = 1;
}


static void check_lowres()
{
	static int have = 0;
	if ( have ) {
		ERROR("You cannot use --rmin and --lowres at the same time.\n");
		exit(1);
	}
	have = 1;
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *afile = NULL;
	char *bfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	int ncom, nrej, nneg, nres, nbij, ncen;
	RefList *list1_acc;
	RefList *list2_acc;
	RefList *list1;
	RefList *list2;
	RefList *list1_raw;
	RefList *list2_raw;
	enum fom fom = FOM_R1I;
	char *cellfile = NULL;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	double rmin, rmax;
	Reflection *refl1;
	RefListIterator *iter;
	float sigma_cutoff = -INFINITY;
	int config_ignorenegs = 0;
	int config_zeronegs = 0;
	int config_unity = 0;
	int config_intshells = 0;
	int nshells = 10;
	char *shell_file = NULL;
	double min_I = +INFINITY;
	double max_I = -INFINITY;
	float highres, lowres;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,               10 },
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},
		{"rmin",               1, NULL,                2},
		{"rmax",               1, NULL,                3},
		{"fom",                1, NULL,                4},
		{"sigma-cutoff",       1, NULL,                5},
		{"nshells",            1, NULL,                6},
		{"shell-file",         1, NULL,                7},
		{"highres",            1, NULL,                8},
		{"lowres",             1, NULL,                9},
		{"ignore-negs",        0, &config_ignorenegs,  1},
		{"zero-negs",          0, &config_zeronegs,    1},
		{"intensity-shells",   0, &config_intshells,   1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hy:p:u",
	                        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 10 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'p' :
			cellfile = strdup(optarg);
			break;

			case 'u' :
			config_unity = 1;
			break;

			case 0 :
			break;

			case 2 :
			check_lowres();
			if ( sscanf(optarg, "%e", &rmin_fix) != 1 ) {
				ERROR("Invalid value for --rmin\n");
				return 1;
			}
			break;

			case 3 :
			check_highres();
			if ( sscanf(optarg, "%e", &rmax_fix) != 1 ) {
				ERROR("Invalid value for --rmax\n");
				return 1;
			}
			break;

			case 4 :
			fom = get_fom(optarg);
			break;

			case 5 :
			if ( sscanf(optarg, "%f", &sigma_cutoff) != 1 ) {
				ERROR("Invalid value for --sigma-cutoff\n");
				return 1;
			}
			break;

			case 6 :
			if ( sscanf(optarg, "%i", &nshells) != 1 ) {
				ERROR("Invalid value for --nshells\n");
				return 1;
			}
			break;

			case 7 :
			shell_file = strdup(optarg);
			break;

			case 8 :
			check_highres();
			if ( sscanf(optarg, "%e", &highres) != 1 ) {
				ERROR("Invalid value for --highres\n");
				return 1;
			}
			rmax_fix = 1.0 / (highres/1e10);
			break;

			case 9 :
			check_lowres();
			if ( sscanf(optarg, "%e", &lowres) != 1 ) {
				ERROR("Invalid value for --lowres\n");
				return 1;
			}
			rmin_fix = 1.0 / (lowres/1e10);
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( argc != (optind+2) ) {
		ERROR("Please provide exactly two HKL files to compare.\n");
		return 1;
	}

	if ( !config_ignorenegs && !config_zeronegs ) {
		switch ( fom )
		{
			case FOM_R1F :
			ERROR("Your chosen figure of merit involves converting"
			      " intensities to structure factors, but you have"
			      " not specified how to handle negative"
			      " intensities.\n");
			ERROR("Please try again with --ignore-negs or"
			      " --zero-negs.\n");
			exit(1);

			case FOM_R2 :
			case FOM_R1I :
			case FOM_RSPLIT :
			case FOM_CC :
			case FOM_CCSTAR :
			case FOM_CCANO :
			case FOM_CRDANO :
			case FOM_RANO :
			case FOM_RANORSPLIT :
			break;
		}
	}

	if ( sym_str == NULL ) {
		sym_str = strdup("1");
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	if ( is_centrosymmetric(sym) ) {
		switch ( fom )
		{
			case FOM_R1F :
			case FOM_R2 :
			case FOM_R1I :
			case FOM_RSPLIT :
			case FOM_CC :
			case FOM_CCSTAR :
			break;

			case FOM_CCANO :
			case FOM_CRDANO :
			case FOM_RANO :
			case FOM_RANORSPLIT :
			ERROR("You are trying to measure an anomalous signal in"
			      " a centrosymmetric point group.\n");
			ERROR("This is a silly thing to do, and I'm refusing to"
			      " help you do it.\n");
			ERROR("Please review your earlier processing steps and"
			      " try again using a non-centrosymmetric point"
			      " group for '-y'.\n");
			return 1;
		}
	}

	afile = strdup(argv[optind++]);
	bfile = strdup(argv[optind]);

	if ( cellfile == NULL ) {
		ERROR("You must provide a unit cell.\n");
		exit(1);
	}

	if ( shell_file == NULL ) shell_file = strdup("shells.dat");

	cell = load_cell_from_file(cellfile);
	free(cellfile);

	list1_raw = read_reflections(afile);
	if ( list1_raw == NULL ) {
		ERROR("Couldn't read file '%s'\n", afile);
		return 1;
	}

	list2_raw = read_reflections(bfile);
	if ( list2_raw == NULL ) {
		ERROR("Couldn't read file '%s'\n", bfile);
		return 1;
	}

	/* Check that the intensities have the correct symmetry */
	if ( check_list_symmetry(list1_raw, sym) ) {
		ERROR("The first input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(sym));
		return 1;
	}
	if ( check_list_symmetry(list2_raw, sym) ) {
		ERROR("The second input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(sym));
		return 1;
	}

	resolution_limits(list1_raw, cell, &rmin, &rmax);
	STATUS("%s: %i reflections, resolution range %.2f to %.2f Angstroms"
	       " (%.5f to %.5f nm^-1).\n", afile,
	       num_reflections(list1_raw),
	       1e10/rmin, 1e10/rmax, rmin/1e9, rmax/1e9);

	resolution_limits(list2_raw, cell, &rmin, &rmax);
	STATUS("%s: %i reflections, resolution range %.2f to %.2f Angstroms"
	       " (%.5f to %.5f nm^-1).\n", bfile,
	       num_reflections(list2_raw),
	       1e10/rmin, 1e10/rmax, rmin/1e9, rmax/1e9);

	list1 = asymmetric_indices(list1_raw, sym);
	list2 = asymmetric_indices(list2_raw, sym);

	/* Select reflections to be used */
	ncom = 0;
	nrej = 0;
	nneg = 0;
	nres = 0;
	nbij = 0;
	ncen = 0;
	list1_acc = reflist_new();
	list2_acc = reflist_new();
	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		double val1, val2;
		double esd1, esd2;
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

		if ( (val1 < sigma_cutoff * esd1)
		  || (val2 < sigma_cutoff * esd2) )
		{
			nrej++;
			continue;
		}

		if ( config_ignorenegs && ((val1 < 0.0) || (val2 < 0.0)) ) {
			nneg++;
			continue;
		}

		if ( config_zeronegs ) {
			int d = 0;
			if ( val1 < 0.0 ) {
				val1 = 0.0;
				d = 1;
			}
			if ( val2 < 0.0 ) {
				val2 = 0.0;
				d = 1;
			}
			if ( d ) nneg++;
		}

		if ( rmin_fix > 0.0 ) {
			double res = 2.0*resolution(cell, h, k, l);
			if ( res < rmin_fix ) {
				nres++;
				continue;
			}
		}

		if ( rmax_fix > 0.0 ) {
			double res = 2.0*resolution(cell, h, k, l);
			if ( res > rmax_fix ) {
				nres++;
				continue;
			}
		}

		if ( (fom == FOM_CCANO) || (fom == FOM_CRDANO)
		  || (fom == FOM_RANO) || (fom == FOM_RANORSPLIT) )
		{
			Reflection *refl1_bij = NULL;
			Reflection *refl2_bij = NULL;
			signed int hb, kb, lb;

			if ( is_centric(h, k, l, sym) ) {
				ncen++;
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
				nbij++;
				continue;
			}
		}

		refl1_acc = add_refl(list1_acc, h, k, l);
		copy_data(refl1_acc, refl1);
		set_intensity(refl1_acc, val1);

		refl2_acc = add_refl(list2_acc, h, k, l);
		copy_data(refl2_acc, refl2);
		set_intensity(refl2_acc, val2);

		if ( val1 > max_I ) max_I = val1;
		if ( val1 < min_I ) min_I = val1;

		ncom++;

	}

	gsl_set_error_handler_off();

	if ( nrej > 0 ) {
		STATUS("Discarded %i reflection pairs because either or both"
		       " versions had I/sigma(I) < %f.\n", nrej, sigma_cutoff);
	}

	if ( config_ignorenegs && (nneg > 0) ) {
		STATUS("Discarded %i reflection pairs because either or both"
		       " versions had negative intensities.\n", nneg);
	}

	if ( config_zeronegs && (nneg > 0) ) {
		STATUS("For %i reflection pairs, either or both versions had"
		       " negative intensities which were set to zero.\n", nneg);
	}

	if ( nres > 0 ) {
		STATUS("%i reflection pairs rejected because either or both"
		       " versions were outside the resolution range.\n", nres);
	}

	if ( nbij > 0 ) {
		STATUS("%i reflection pairs rejected because either or both"
		       " versions did not have Bijvoet partners.\n", nres);
	}

	if ( ncen > 0 ) {
		STATUS("%i reflection pairs rejected because they were"
		       " centric.\n", ncen);
	}

	STATUS("%i reflection pairs accepted.\n", ncom);

	resolution_limits(list1_acc, cell, &rmin, &rmax);
	resolution_limits(list2_acc, cell, &rmin, &rmax);
	STATUS("Accepted resolution range: %f to %f nm^-1"
	       " (%.2f to %.2f Angstroms).\n",
	       rmin/1e9, rmax/1e9, 1e10/rmin, 1e10/rmax);

	reflist_free(list1_raw);
	reflist_free(list2_raw);
	reflist_free(list1);
	reflist_free(list2);

	if ( rmin_fix >= 0.0 ) {
		rmin = rmin_fix;
	}
	if ( rmax_fix >= 0.0 ) {
		rmax = rmax_fix;
	}
	if ( (rmin_fix>=0.0) || (rmax_fix>=0.0) ) {
		STATUS("Fixed resolution range: %f to %f nm^-1"
		       " (%.2f to %.2f Angstroms).\n",
		       rmin/1e9, rmax/1e9, 1e10/rmin, 1e10/rmax);
	}
	do_fom(list1_acc, list2_acc, cell, rmin, rmax, fom, config_unity,
	       nshells, shell_file, config_intshells, min_I, max_I, sym);

	free(shell_file);
	reflist_free(list1_acc);
	reflist_free(list2_acc);

	return 0;
}
