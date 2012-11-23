/*
 * compare_hkl.c
 *
 * Compare reflection lists
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
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
	FOM_CCSTAR
};


static enum fom get_fom(const char *s)
{
	if ( strcasecmp(s, "r1i") == 0 ) return FOM_R1I;
	if ( strcasecmp(s, "r1f") == 0 ) return FOM_R1F;
	if ( strcasecmp(s, "r2") == 0 ) return FOM_R2;
	if ( strcasecmp(s, "rsplit") == 0 ) return FOM_RSPLIT;
	if ( strcasecmp(s, "cc") == 0 ) return FOM_CC;
	if ( strcasecmp(s, "ccstar") == 0 ) return FOM_CCSTAR;

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
"  -p, --pdb=<filename>       PDB file to use.\n"
"      --fom=<FoM>            Calculate this figure of merit  Choose from:.\n"
"                              R1I, R1F, R2, Rsplit,\n"
"                              CCF, CCI, CCFstar, CCIstar.\n"
"      --nshells=<n>          Use <n> resolution shells.\n"
"  -u                         Force scale factor to 1.\n"
"      --shell-file=<file>    Write resolution shells to <file>.\n"
"\n"
"You can control which reflections are included in the calculation:\n"
"\n"
"      --ignore-negs          Ignore reflections with negative intensities.\n"
"      --zero-negs            Set negative intensities to zero.\n"
"      --sigma-cutoff=<n>     Discard reflections with I/sigma(I) < n.\n"
"      --rmin=<res>           Set a lower resolution limit (m^-1).\n"
"      --rmax=<res>           Set an upper resolution limit (m^-1).\n"
"\n"
"  -h, --help                 Display this help message.\n"
);
}


struct fom_context
{
	enum fom fom;
	int nshells;

	/* For R-factors */
	double *num;
	double *den;

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

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
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
                       int bin)
{
	double f1, f2;

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

	}
}


static double fom_overall(struct fom_context *fctx)
{
	double overall_num = INFINITY;
	double overall_den = 0.0;
	int i;
	double *overall_vec1;
	double *overall_vec2;
	int overall_n;
	double cc = INFINITY;

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
		overall_num = 0.0;
		overall_den = 0.0;
		for ( i=0; i<fctx->nshells; i++ ) {
			overall_num += fctx->num[i];
			overall_den += fctx->den[i];
		}
		break;

		case FOM_CC :
		case FOM_CCSTAR :
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
		return cc;

		case FOM_CCSTAR :
		return sqrt((2.0*cc)/(1.0+cc));

	}

	ERROR("This point is never reached.\n");
	abort();
}


static double fom_shell(struct fom_context *fctx, int i)
{
	double cc;

	switch ( fctx->fom ) {

		case FOM_R1I :
		case FOM_R1F :
		return fctx->num[i]/fctx->den[i];

		case FOM_R2 :
		return sqrt(fctx->num[i]/fctx->den[i]);

		case FOM_RSPLIT :
		return 2.0*(fctx->num[i]/fctx->den[i]) / sqrt(2.0);

		case FOM_CC :
		return gsl_stats_correlation(fctx->vec1[i], 1, fctx->vec2[i], 1,
		                             fctx->n[i]);

		case FOM_CCSTAR :
		cc = gsl_stats_correlation(fctx->vec1[i], 1, fctx->vec2[i], 1,
		                           fctx->n[i]);
		return sqrt((2.0*cc)/(1.0+cc));

	}

	ERROR("This point is never reached.\n");
	abort();
}


static void do_fom(RefList *list1, RefList *list2, UnitCell *cell,
                   double rmin, double rmax, enum fom fom,
                   int config_unity, int nshells, const char *filename)
{
	int *cts;
	double *rmins;
	double *rmaxs;
	double total_vol, vol_per_shell;
	int i;
	Reflection *refl1;
	RefListIterator *iter;
	FILE *fh;
	double scale;
	struct fom_context *fctx;

	cts = malloc(nshells*sizeof(int));
	rmins = malloc(nshells*sizeof(double));
	rmaxs = malloc(nshells*sizeof(double));
	fctx = init_fom(fom, num_reflections(list1), nshells);

	if ( (fctx==NULL) || (cts==NULL) || (rmins==NULL) || (rmaxs==NULL) )
	{
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

	/* Calculate the resolution bins */
	total_vol = pow(rmax, 3.0) - pow(rmin, 3.0);
	vol_per_shell = total_vol / nshells;
	rmins[0] = rmin;
	for ( i=1; i<nshells; i++ ) {

		double r;

		r = vol_per_shell + pow(rmins[i-1], 3.0);
		r = pow(r, 1.0/3.0);

		/* Shells of constant volume */
		rmaxs[i-1] = r;
		rmins[i] = r;

		/* Shells of constant thickness */
		//rmins[i] = rmins[i-1] + (rmax-rmin)/NBINS;
		//rmaxs[i-1] = rmins[i-1] + (rmax-rmin)/NBINS;

	}
	rmaxs[nshells-1] = rmax;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		double d;
		int bin;
		double i1, i2;
		int j;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		d = 2.0 * resolution(cell, h, k, l);

		bin = -1;
		for ( j=0; j<nshells; j++ ) {
			if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
				bin = j;
				break;
			}
		}

		/* Allow for slight rounding errors */
		if ( (bin == -1) && (d <= rmins[0]) ) bin = 0;
		assert(bin != -1);

		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;

		i1 = get_intensity(refl1);
		i2 = scale * get_intensity(refl2);

		add_to_fom(fctx, i1, i2, bin);
		cts[bin]++;
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

	}

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		return;
	}

	switch ( fom ) {

		case FOM_R1I :
		fprintf(fh, "1/d centre    R1(I)/%%       nref      d / A\n");
		break;

		case FOM_R1F :
		fprintf(fh, "1/d centre    R1(F)/%%       nref      d / A\n");
		break;

		case FOM_R2 :
		fprintf(fh, "1/d centre       R2/%%       nref      d / A\n");
		break;

		case FOM_RSPLIT :
		fprintf(fh, "1/d centre   Rsplit/%%       nref      d / A\n");
		break;

		case FOM_CC :
		fprintf(fh, "1/d centre         CC       nref      d / A\n");
		break;

		case FOM_CCSTAR :
		fprintf(fh, "1/d centre        CC*       nref      d / A\n");
		break;

	}

	for ( i=0; i<nshells; i++ ) {

		double r, cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;

		r = fom_shell(fctx, i);

		switch ( fom ) {

		case FOM_R1I :
		case FOM_R1F :
		case FOM_R2 :
		case FOM_RSPLIT :
		fprintf(fh, "%10.3f %10.2f %10i %10.2f\n",
		        cen*1.0e-9, r*100.0, cts[i], (1.0/cen)*1e10);

		break;

		case FOM_CC :
		case FOM_CCSTAR :
		fprintf(fh, "%10.3f %10.7f %10i %10.2f\n",
		        cen*1.0e-9, r, cts[i], (1.0/cen)*1e10);

		break;

	}

	}

	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *afile = NULL;
	char *bfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	int ncom, nrej, nneg, nres;
	RefList *list1_acc;
	RefList *list2_acc;
	RefList *list1;
	RefList *list2;
	RefList *list1_raw;
	RefList *list2_raw;
	enum fom fom = FOM_R1I;
	char *pdb = NULL;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	double rmin, rmax;
	Reflection *refl1;
	RefListIterator *iter;
	float sigma_cutoff = -INFINITY;
	int config_ignorenegs = 0;
	int config_zeronegs = 0;
	int config_unity = 0;
	int nshells = 10;
	char *shell_file = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},
		{"rmin",               1, NULL,                2},
		{"rmax",               1, NULL,                3},
		{"fom",                1, NULL,                4},
		{"sigma-cutoff",       1, NULL,                5},
		{"nshells",            1, NULL,                6},
		{"shell-file",         1, NULL,                7},
		{"ignore-negs",        0, &config_ignorenegs,  1},
		{"zero-negs",          0, &config_zeronegs,    1},
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

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'p' :
			pdb = strdup(optarg);
			break;

			case 'u' :
			config_unity = 1;
			break;

			case 0 :
			break;

			case 2 :
			if ( sscanf(optarg, "%e", &rmin_fix) != 1 ) {
				ERROR("Invalid value for --rmin\n");
				return 1;
			}
			break;

			case 3 :
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
			break;
		}
	}

	if ( sym_str == NULL ) {
		sym_str = strdup("1");
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	afile = strdup(argv[optind++]);
	bfile = strdup(argv[optind]);

	if ( pdb == NULL ) {
		ERROR("You must provide a PDB file with the unit cell.\n");
		exit(1);
	}

	if ( shell_file == NULL ) shell_file = strdup("shells.dat");

	cell = load_cell_from_pdb(pdb);
	free(pdb);

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

		refl1_acc = add_refl(list1_acc, h, k, l);
		copy_data(refl1_acc, refl1);
		set_intensity(refl1_acc, val1);

		refl2_acc = add_refl(list2_acc, h, k, l);
		copy_data(refl2_acc, refl2);
		set_intensity(refl2_acc, val2);

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
		STATUS("%i reflection pairs rejected because either or both "
		       " versions were outside the resolution range.\n", nres);
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
	       nshells, shell_file);

	free(shell_file);
	reflist_free(list1_acc);
	reflist_free(list2_acc);

	return 0;
}
