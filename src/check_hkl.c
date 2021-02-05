/*
 * check_hkl.c
 *
 * Characterise reflection lists
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_fit.h>
#include <assert.h>

#include <utils.h>
#include <symmetry.h>
#include <reflist.h>
#include <reflist-utils.h>
#include <cell-utils.h>
#include <fom.h>

#include "version.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.hkl>\n\n", s);
	printf(
"Characterise an intensity list.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"  -y, --symmetry=<sym>       The symmetry of the input file.\n"
"  -p, --pdb=<filename>       Unit cell file to use (PDB or CrystFEL format).\n"
"      --rmin=<res>           Low resolution cutoff (1/d in m^-1).\n"
"      --rmax=<res>           High resolution cutoff (1/d in m^-1).\n"
"      --lowres=<n>           Low resolution cutoff in (d in A).\n"
"      --highres=<n>          High resolution cutoff in (d in A).\n"
"      --sigma-cutoff=<n>     Discard reflections with I/sigma(I) < n.\n"
"      --nshells=<n>          Use <n> resolution shells or bins.\n"
"      --wilson               Calculate a Wilson plot\n"
"      --ltest                Perform an L-test (for twinning)\n"
"      --shell-file=<file>    Write results table to <file>.\n"
"      --ignore-negs          Ignore reflections with negative intensities.\n"
"      --zero-negs            Set negative intensities to zero.\n"
"\n");
}


static int add_ltest(RefList *list, double i1, int *bins, int nbins,
                     double step, const SymOpList *sym, double *lt, double *l2t,
                     signed int h1, signed int k1, signed int l1,
                     signed int h2, signed int k2, signed int l2)
{
	Reflection *refl;
	double i2, L;
	int bin;

	if ( SERIAL(h1, k1, l1) > SERIAL(h2, k2, l2) ) return 0;

	refl = find_refl(list, h2, k2, l2);
	if ( refl == NULL ) {
		signed int h, k, l;
		if ( !find_equiv_in_list(list, h2, k2, l2, sym, &h, &k, &l) ) {
			return 0;
		}
		refl = find_refl(list, h, k, l);
	}

	i2 = get_intensity(refl);
	L = (i1-i2) / (i1+i2);
	if ( isnan(L) ) {
		/* This happens with --zero-negs and two negative intensities,
		 * because L=(0-0)/(0+0) */
		return 0;
	}

	bin = fabs(L)/step;
	if ( (bin < 0) || (isnan(L)) ) {
		bin = 0;
	} else if ( bin >= nbins ) {
		bin = nbins-1;
	}
	bins[bin]++;

	*lt += fabs(L);
	*l2t += pow(L, 2.0);

	return 1;
}


static void l_test(RefList *list, UnitCell *cell, const SymOpList *sym,
                   double rmin_fix, double rmax_fix, int nbins,
                   const char *filename)
{
	Reflection *refl;
	RefListIterator *iter;
	int *bins;
	FILE *fh;
	int npairs, i;
	double tot;
	double lt = 0.0;
	double l2t = 0.0;
	const double step = 1.0/nbins;
	int hd, kd, ld;
	const char cen = cell_get_centering(cell);

	bins = malloc(nbins*sizeof(int));
	if ( bins == NULL ) return;

	for ( i=0; i<nbins; i++ ) bins[i] = 0;

	if ( cen == 'P' ) { hd = 1;  kd = 1;  ld = 1; }
	if ( cen == 'R' ) { hd = 1;  kd = 1;  ld = 1; }

	if ( cen == 'A' ) { hd = 1;  kd = 2;  ld = 2; }
	if ( cen == 'B' ) { hd = 2;  kd = 1;  ld = 2; }
	if ( cen == 'C' ) { hd = 2;  kd = 2;  ld = 1; }

	if ( cen == 'I' ) { hd = 2;  kd = 2;  ld = 2; }
	if ( cen == 'F' ) { hd = 2;  kd = 2;  ld = 2; }

	/* Obverse setting */
	if ( cen == 'H' ) { hd = 3;  kd = 3;  ld = 3; }

	npairs = 0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double i1;

		get_indices(refl, &h, &k, &l);
		i1 = get_intensity(refl);

		npairs += add_ltest(list, i1, bins, nbins, step, sym, &lt, &l2t,
		                    h, k, l, h-hd, k, l);
		npairs += add_ltest(list, i1, bins, nbins, step, sym, &lt, &l2t,
		                    h, k, l, h+hd, k, l);
		npairs += add_ltest(list, i1, bins, nbins, step, sym, &lt, &l2t,
		                    h, k, l, h, k-kd, l);
		npairs += add_ltest(list, i1, bins, nbins, step, sym, &lt, &l2t,
		                    h, k, l, h, k+kd, l);
		npairs += add_ltest(list, i1, bins, nbins, step, sym, &lt, &l2t,
		                    h, k, l, h, k, l-ld);
		npairs += add_ltest(list, i1, bins, nbins, step, sym, &lt, &l2t,
		                    h, k, l, h, k, l+ld);
	}
	STATUS("%i pairs\n", npairs);
	STATUS("<|L|> = %.3f (ideal untwinned %.3f, twinned %.3f)\n",
	       lt/npairs, 1.0/2.0, 3.0/8.0);
	STATUS("<L^2> = %.3f (ideal untwinned %.3f, twinned %.3f)\n",
	       l2t/npairs, 1.0/3.0, 1.0/5.0);

	fh = fopen(filename, "w");
	if ( fh == NULL ) return;
	tot = 0.0;
	fprintf(fh, "  |L|  N(|L|) untwinned twinned\n");
	fprintf(fh, "%.3f %7.3f %9.3f %7.3f\n", 0.0, tot, 0.0, 0.0);
	for ( i=0; i<nbins; i++ ) {
		double l = (i+1)*step;
		tot += (double)bins[i]/npairs;
		fprintf(fh, "%.3f %7.3f %9.3f %7.3f\n",
		        l, tot, l, l*(3-l*l)/2.0);
	}
	fclose(fh);

	free(bins);
}


static double get_sfac(char el, double s)
{
	double sfac[9];
	double sf;
	double s2 = pow(s/1e10, 2.0);  /* s^2 in A^-2 */

	switch ( el ) {

		case 'C' :
		sfac[0] = 2.31000;   sfac[1] = 20.8439;
		sfac[2] = 1.02000;   sfac[3] = 10.2075;
		sfac[4] = 1.58860;   sfac[5] = 0.568700;
		sfac[6] = 0.865000;  sfac[7] = 51.6512;
		sfac[8] = 0.215600;
		break;

		case 'N' :
		sfac[0] = 12.2126;  sfac[1] = 0.005700;
		sfac[2] = 3.13220;  sfac[3] = 9.89330;
		sfac[4] = 2.01250;  sfac[5] = 28.9975;
		sfac[6] = 1.16630;  sfac[7] = 0.582600;
		sfac[8] = -11.529;
		break;

		case 'O' :
		sfac[0] = 3.04850;  sfac[1] = 13.2771;
		sfac[2] = 2.28680;  sfac[3] = 5.70110;
		sfac[4] = 1.54630;  sfac[5] = 0.323900;
		sfac[6] = 0.867000; sfac[7] = 322.9098;
		sfac[8] = 0.250800;
		break;

		case 'H' :
		sfac[0] = 0.489918;  sfac[1] = 20.6593;
		sfac[2] = 0.262003;  sfac[3] = 7.74039;
		sfac[4] = 0.196767;  sfac[5] = 49.5519;
		sfac[6] = 0.049879;  sfac[7] = 2.20159;
		sfac[8] = 0.001305;
		break;

		default :
		ERROR("Unrecognised atom '%c'\n", el);
		abort();
	}

	sf  = sfac[0] * exp(-sfac[1]*s2);
	sf += sfac[2] * exp(-sfac[3]*s2);
	sf += sfac[4] * exp(-sfac[5]*s2);
	sf += sfac[6] * exp(-sfac[7]*s2);
	sf += sfac[8];

	return sf;
}


static void wilson_plot(RefList *list, UnitCell *cell, const SymOpList *sym,
                        double rmin_fix, double rmax_fix, int nbins,
                        const char *filename)
{
	double rmin, rmax;
	double s2min, s2max, s2step;
	Reflection *refl;
	RefListIterator *iter;
	FILE *fh;
	double *plot_i, *s2;
	int *plot_n;
	int i, ngen;
	SymOpMask *mask;

	resolution_limits(list, cell, &rmin, &rmax);
	STATUS("1/d goes from %f to %f nm^-1\n", rmin/1e9, rmax/1e9);

	/* Widen the range just a little bit */
	rmin -= 0.001e9;
	rmax += 0.001e9;

	/* Fixed resolution shells if needed */
	if ( rmin_fix > 0.0 ) rmin = rmin_fix;
	if ( rmax_fix > 0.0 ) rmax = rmax_fix;

	s2min = pow(rmin/2.0, 2.0);
	s2max = pow(rmax/2.0, 2.0);
	s2step = (s2max - s2min)/nbins;

	plot_i = malloc(nbins*sizeof(double));
	if ( plot_i == NULL ) return;
	for ( i=0; i<nbins; i++ ) plot_i[i] = 0.0;

	s2 = malloc(nbins*sizeof(double));
	if ( s2 == NULL ) return;

	plot_n = calloc(nbins, sizeof(int));
	if ( plot_n == NULL ) return;

	ngen = num_equivs(sym, NULL);
	mask = new_symopmask(sym);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double s, intensity, E;
		int bin;
		int e;

		get_indices(refl, &h, &k, &l);

		s = resolution(cell, h, k, l);  /* This gives "s" directly */
		intensity = get_intensity(refl);

		bin = (pow(s, 2.0) - s2min)/s2step;

		special_position(sym, mask, h, k, l);
		e = ngen / num_equivs(sym, mask);

		/* Average atoms per residue from Rupp BMC 1st ed p356 */
		E  = 5.00*get_sfac('C', s);
		E += 1.35*get_sfac('N', s);
		E += 1.50*get_sfac('O', s);
		E += 8.00*get_sfac('H', s);

		if ( bin == nbins ) bin = nbins-1;
		assert(bin < nbins);

		plot_i[bin] += intensity / (e*E);
		plot_n[bin]++;

	}

	free_symopmask(mask);

	for ( i=0; i<nbins; i++ ) {
		plot_i[i] = log(plot_i[i] / plot_n[i]);
		s2[i] = s2min + (i+0.5)*s2step;
	}

	fh = fopen(filename, "w");
	if ( fh == NULL ) return;
	fprintf(fh, "  n  s^2/A^-2      d/A  ln <I>/eE  nrefl\n");
	for ( i=0; i<nbins; i++ ) {
		fprintf(fh, "%3i  %8.6f %8.4f %10f %6i\n", i,
		        s2[i]/1e20, 0.5e10/sqrt(s2[i]), plot_i[i], plot_n[i]);
	}
	fclose(fh);

	if ( rmax < 3.125e9 ) {
		ERROR("Resolution too low to estimate B factor\n");
	} else {
		int bs = (pow(3.125e9/2.0, 2.0) - s2min)/s2step + 1;
		double lnk, minus2B, cov00, cov01, cov11, sumsq;
		double B;
		if ( nbins - bs < 3 ) {
			ERROR("Not enough bins to estimate B factor\n");
			ERROR("Resolution of data, or number of bins, is too "
			      "low.\n");
		} else {
			double *s2fit;
			double *lnifit;
			int nbfit = 0;
			s2fit = malloc(nbins*sizeof(double));
			lnifit = malloc(nbins*sizeof(double));
			if ( (s2fit==NULL) || (lnifit==NULL) ) return;
			for ( i=0; i<nbins-bs; i++ ) {
				if ( isnan(plot_i[bs+i]) ) continue;
				s2fit[i] = s2[bs+i];
				lnifit[i] = plot_i[bs+i];
				nbfit++;
			}
			if ( nbfit < 3 ) {
				ERROR("Too many bits had invalid values.\n");
				return;
			}
			gsl_fit_linear(s2fit, 1, lnifit, 1, nbfit,
				       &lnk, &minus2B, &cov00, &cov01, &cov11,
				       &sumsq);
			B = -minus2B/2.0;
			STATUS("ln k = %.2f\n", lnk);
			STATUS("B = %.2f A^2\n", B*1e20);
		}
	}


	free(plot_i);
	free(plot_n);
}


static void plot_shells(RefList *list, UnitCell *cell, const SymOpList *sym,
                        double rmin_fix, double rmax_fix, int nshells,
			const char *shell_file)
{
	double rmin, rmax;
	int i;
	FILE *fh;
	struct fom_shells *shells;
	struct fom_context *nmeas_ctx;
	struct fom_context *red_ctx;
	struct fom_context *snr_ctx;
	struct fom_context *mean_ctx;
	struct fom_context *compl_ctx;

	fh = fopen(shell_file, "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", shell_file);
		return;
	}

	resolution_limits(list, cell, &rmin, &rmax);
	STATUS("1/d goes from %f to %f nm^-1\n", rmin/1e9, rmax/1e9);

	/* Fixed resolution shells if needed */
	if ( rmin_fix > 0.0 ) rmin = rmin_fix;
	if ( rmax_fix > 0.0 ) rmax = rmax_fix;

	shells = fom_make_resolution_shells(rmin, rmax,  nshells);

	STATUS("Overall values within specified resolution range:\n");

	nmeas_ctx = fom_calculate(list, NULL, cell, shells,
	                          FOM_NUM_MEASUREMENTS, 0, sym);
	red_ctx = fom_calculate(list, NULL, cell, shells,
	                        FOM_REDUNDANCY, 0, sym);
	snr_ctx = fom_calculate(list, NULL, cell, shells,
	                        FOM_SNR, 0, sym);
	mean_ctx = fom_calculate(list, NULL, cell, shells,
	                         FOM_MEAN_INTENSITY, 0, sym);
	compl_ctx = fom_calculate(list, NULL, cell, shells,
	                          FOM_COMPLETENESS, 0, sym);

	STATUS("%.0f measurements in total.\n",
	       fom_overall_value(nmeas_ctx));
	STATUS("%li reflections in total.\n",
	       fom_overall_num_reflections(compl_ctx));
	STATUS("%li reflections possible.\n",
	       fom_overall_num_possible(compl_ctx));
	STATUS("Overall <snr> = %f\n", fom_overall_value(snr_ctx));
	STATUS("Overall redundancy = %f measurements/unique reflection\n",
	       fom_overall_value(red_ctx));
	STATUS("Overall completeness = %f %%\n",
	       100.0*fom_overall_value(compl_ctx));

	fprintf(fh, "Center 1/nm  # refs Possible  Compl       "
		    "Meas   Red   SNR     Mean I     d(A)    "
		    "Min 1/nm   Max 1/nm\n");
	for ( i=0; i<nshells; i++ ) {

		long int measured, possible;

		measured = fom_shell_num_reflections(compl_ctx, i);
		possible = fom_shell_num_possible(compl_ctx, i);

		fprintf(fh, "%10.3f %8li %8li %6.2f %10.0f %5.1f"
		            " %5.2f %10.2f %8.2f  %10.3f %10.3f\n",
		        fom_shell_centre(shells, i)*1.0e-9,
		        measured,
		        possible,
		        100.0*fom_shell_value(compl_ctx, i),
		        fom_shell_value(nmeas_ctx, i),
		        fom_shell_value(red_ctx, i),
		        fom_shell_value(snr_ctx, i),
		        fom_shell_value(mean_ctx, i),
		        (1.0/fom_shell_centre(shells, i))*1e10,
			shells->rmins[i]*1.0e-9,
			shells->rmaxs[i]*1.0e-9);

	}

	fclose(fh);

	STATUS("Resolution shell information written to %s.\n", shell_file);
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
	char *file = NULL;
	char *sym_str = NULL;
	char *sym_str_fromfile = NULL;
	SymOpList *sym;
	RefList *raw_list;
	RefList *list;
	char *cellfile = NULL;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	float sigma_cutoff = -INFINITY;
	int nshells = 10;
	int have_nshells = 0;
	char *shell_file = NULL;
	int wilson = 0;
	int ltest = 0;
	int ignorenegs = 0;
	int zeronegs = 0;
	float highres, lowres;
	struct fom_rejections rej;

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                9 },
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},

		{"rmin",               1, NULL,                2},
		{"rmax",               1, NULL,                3},
		{"sigma-cutoff",       1, NULL,                4},
		{"nshells",            1, NULL,                5},
		{"shell-file",         1, NULL,                6},
		{"highres",            1, NULL,                7},
		{"lowres",             1, NULL,                8},

		{"wilson",             0, &wilson,             1},
		{"ltest",              0, &ltest,              1},
		{"ignore-negs",        0, &ignorenegs,         1},
		{"zero-negs",          0, &zeronegs,           1},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hy:p:", longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 9 :
			printf("CrystFEL: %s\n",
			       crystfel_version_string());
			printf("%s\n",
			       crystfel_licence_string());
			return 0;

			case 'y' :
			sym_str = strdup(optarg);
			break;

			case 'p' :
			cellfile = strdup(optarg);
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
			if ( sscanf(optarg, "%f", &sigma_cutoff) != 1 ) {
				ERROR("Invalid value for --sigma-cutoff\n");
				return 1;
			}
			STATUS("WARNING: You are using --sigma-cutoff.  "
			       "Be aware that the figures of merit will not "
			       "reflect the entire data set!\n");
			break;

			case 5 :
			if ( sscanf(optarg, "%i", &nshells) != 1 ) {
				ERROR("Invalid value for --nshells\n");
				return 1;
			}
			have_nshells = 1;
			break;

			case 6 :
			shell_file = strdup(optarg);
			break;

			case 7 :
			check_highres();
			if ( sscanf(optarg, "%e", &highres) != 1 ) {
				ERROR("Invalid value for --highres\n");
				return 1;
			}
			rmax_fix = 1.0 / (highres/1e10);
			break;

			case 8 :
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

	if ( argc != (optind+1) ) {
		ERROR("Please provide exactly one HKL file to check.\n");
		return 1;
	}

	if ( !ltest && (ignorenegs || zeronegs) ) {
		ERROR("WARNING: You are using --zero-negs or --ignore-negs "
		      "even though it's not required.\n");
		ERROR("The figures of merit will not reflect the entire data "
		      "set!\n");
	}

	file = strdup(argv[optind++]);

	if ( cellfile == NULL ) {
		ERROR("You need to provide a unit cell.\n");
		return 1;
	}
	cell = load_cell_from_file(cellfile);
	if ( cell == NULL ) {
		ERROR("Failed to load cell.\n");
		return 1;
	}
	free(cellfile);

	raw_list = read_reflections_2(file, &sym_str_fromfile);
	if ( raw_list == NULL ) {
		ERROR("Couldn't read file '%s'\n", file);
		return 1;
	}
	free(file);

	if ( sym_str == NULL ) {
		if ( sym_str_fromfile != NULL ) {
			STATUS("Using symmetry from reflection file: %s\n",
			       sym_str_fromfile);
			sym_str = sym_str_fromfile;
		} else {
			sym_str = strdup("1");
		}
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	if ( shell_file == NULL ) shell_file = strdup("shells.dat");

	/* Check that the intensities have the correct symmetry */
	if ( check_list_symmetry(raw_list, sym) ) {
		ERROR("The input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(sym));
		if ( cell_get_lattice_type(cell) == L_MONOCLINIC ) {
			ERROR("You may need to specify the unique axis in your "
			      "point group.  The default is unique axis c.\n");
			ERROR("See 'man crystfel' for more details.\n");
		}
		return 1;
	}

	/* Reject some reflections */
	rej = fom_select_reflections(raw_list, &list, cell, sym,
	                             rmin_fix, rmax_fix, sigma_cutoff,
	                             ignorenegs, zeronegs, 0);

	STATUS("Discarded %i reflections (out of %i) with I/sigma(I) < %f\n",
	       rej.low_snr, num_reflections(raw_list), sigma_cutoff);
	reflist_free(raw_list);

	if ( rej.negative_deleted > 0 ) {
		STATUS("Discarded %i reflections because they had negative "
		       "intensities.\n", rej.negative_deleted);
	}

	if ( rej.negative_zeroed > 0 ) {
		STATUS("Set %i negative intensities to zero\n",
		       rej.negative_zeroed);
	}

	if ( rej.outside_resolution_range > 0 ) {
		STATUS("%i reflections rejected because they were outside the "
		       "resolution range.\n", rej.outside_resolution_range);
	}

	if ( rej.nan_inf_value ) {
		STATUS("WARNING: %i reflections had infinite or invalid values"
		       " of I or sigma(I).\n", rej.nan_inf_value);
	}

	if ( wilson ) {
		if ( !have_nshells ) nshells = 50;
		wilson_plot(list, cell, sym, rmin_fix, rmax_fix, nshells,
		            shell_file);
	} else if ( ltest ) {
		if ( !have_nshells ) nshells = 50;
		if ( !ignorenegs && !zeronegs ) {
			ERROR("For the L-test you must specify either"
			       "--ignore-negs or --zero-negs.\n");
			return 1;
		}
		l_test(list, cell, sym, rmin_fix, rmax_fix, nshells,
		       shell_file);
	} else {
		plot_shells(list, cell, sym, rmin_fix, rmax_fix, nshells,
		            shell_file);
	}

	free_symoplist(sym);
	reflist_free(list);
	cell_free(cell);
	free(shell_file);

	return 0;
}
