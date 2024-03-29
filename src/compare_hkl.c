/*
 * compare_hkl.c
 *
 * Compare reflection lists
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>

#include <utils.h>
#include <symmetry.h>
#include <reflist-utils.h>
#include <cell-utils.h>
#include <fom.h>

#include "version.h"


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
"			       CCano, CRDano, Rano, Rano/Rsplit, d1sig,\n"
"                              d2sig\n"
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
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
);
}


static enum fom_type fom_type_from_string(const char *s)
{
	if ( strcasecmp(s, "r1i") == 0 ) return FOM_R1I;
	if ( strcasecmp(s, "r1f") == 0 ) return FOM_R1F;
	if ( strcasecmp(s, "r2") == 0 ) return FOM_R2;
	if ( strcasecmp(s, "rsplit") == 0 ) return FOM_RSPLIT;
	if ( strcasecmp(s, "cc") == 0 ) return FOM_CC;
	if ( strcasecmp(s, "cc1/2") == 0 ) return FOM_CC;
	if ( strcasecmp(s, "cchalf") == 0 ) return FOM_CC;
	if ( strcasecmp(s, "ccstar") == 0 ) return FOM_CCSTAR;
	if ( strcasecmp(s, "cc*") == 0 ) return FOM_CCSTAR;
	if ( strcasecmp(s, "ccano") == 0 ) return FOM_CCANO;
	if ( strcasecmp(s, "crdano") == 0 ) return FOM_CRDANO;
	if ( strcasecmp(s, "rano") == 0 ) return FOM_RANO;
	if ( strcasecmp(s, "rano/rsplit") == 0 ) return FOM_RANORSPLIT;
	if ( strcasecmp(s, "d1sig") == 0 ) return FOM_D1SIG;
	if ( strcasecmp(s, "d2sig") == 0 ) return FOM_D2SIG;

	ERROR("Unknown figure of merit '%s'.\n", s);
	exit(1);
}

static void do_fom(RefList *list1, RefList *list2, UnitCell *cell,
                   double rmin, double rmax, enum fom_type fom,
                   int config_unity, int nshells, const char *filename,
                   SymOpList *sym)
{
	struct fom_shells *shells;
	struct fom_context *fctx;
	FILE *fh;
	int i;
	const char *t1, *t2;

	/* Calculate the bins */
	shells = fom_make_resolution_shells(rmin, rmax, nshells);

	if ( shells == NULL ) {
		ERROR("Failed to set up shells.\n");
		return;
	}

	fctx = fom_calculate(list1, list2, cell, shells, fom,
	                     config_unity, sym);

	switch ( fom ) {

		case FOM_R1I :
		STATUS("Overall R1(I) = %.2f %%\n", 100.0*fom_overall_value(fctx));
		break;

		case FOM_R1F :
		STATUS("Overall R1(F) = %.2f %%\n", 100.0*fom_overall_value(fctx));
		break;

		case FOM_R2 :
		STATUS("Overall R(2) = %.2f %%\n", 100.0*fom_overall_value(fctx));
		break;

		case FOM_RSPLIT :
		STATUS("Overall Rsplit = %.2f %%\n", 100.0*fom_overall_value(fctx));
		break;

		case FOM_CC :
		STATUS("Overall CC = %.7f\n", fom_overall_value(fctx));
		break;

		case FOM_CCSTAR :
		STATUS("Overall CC* = %.7f\n", fom_overall_value(fctx));
		break;

		case FOM_CCANO :
		STATUS("Overall CCano = %.7f\n", fom_overall_value(fctx));
		break;

		case FOM_CRDANO :
		STATUS("Overall CRDano = %.7f\n", fom_overall_value(fctx));
		break;

		case FOM_RANO :
		STATUS("Overall Rano =  %.2f %%\n", 100.0*fom_overall_value(fctx));
		break;

		case FOM_RANORSPLIT :
		STATUS("Overall Rano/Rsplit =  %.7f\n", fom_overall_value(fctx));
		break;

		case FOM_D1SIG :
		STATUS("Fraction of differences less than 1 sigma = %.7f %%\n",
		       100.0*fom_overall_value(fctx));
		break;

		case FOM_D2SIG :
		STATUS("Fraction of differences less than 2 sigma = %.7f %%\n",
		       100.0*fom_overall_value(fctx));
		break;

		default :
		ERROR("Unhandled figure of merit type\n");
		break;

	}

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		return;
	}

	t1 = "  1/d centre";
	t2 = "      d / A   Min 1/nm    Max 1/nm";

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

		case FOM_D1SIG :
		fprintf(fh, "%s D<1sigma/%%     nref%s\n", t1, t2);
		break;

		case FOM_D2SIG :
		fprintf(fh, "%s D<2sigma/%%     nref%s\n", t1, t2);
		break;

		default :
		break;

	}

	for ( i=0; i<nshells; i++ ) {

		double r, cen;

		cen = fom_shell_centre(shells, i);
		r = fom_shell_value(fctx, i);

		switch ( fom ) {

			case FOM_R1I :
			case FOM_R1F :
			case FOM_R2 :
			case FOM_RSPLIT :
			case FOM_RANO :
			fprintf(fh, "%10.3f %10.2f %10i %10.2f "
			        "%10.3f  %10.3f\n",
			        cen*1.0e-9, r*100.0,
			        fom_shell_num_reflections(fctx, i),
			        (1.0/cen)*1e10,
			        shells->rmins[i]*1.0e-9,
			        shells->rmaxs[i]*1.0e-9);
			break;

			case FOM_CC :
			case FOM_CCSTAR :
			case FOM_CCANO :
			case FOM_CRDANO :
			fprintf(fh, "%10.3f %10.7f %10i %10.2f "
			        "%10.3f  %10.3f\n",
			        cen*1.0e-9, r,
			        fom_shell_num_reflections(fctx, i),
			        (1.0/cen)*1e10,
			        shells->rmins[i]*1.0e-9,
			        shells->rmaxs[i]*1.0e-9);
			break;

			case FOM_RANORSPLIT :
			fprintf(fh, "%10.3f    %10.7f %10i %10.2f "
			        "%10.3f  %10.3f\n",
			        cen*1.0e-9, r,
			        fom_shell_num_reflections(fctx, i),
			        (1.0/cen)*1e10,
			        shells->rmins[i]*1.0e-9,
			        shells->rmaxs[i]*1.0e-9);
			break;

			case FOM_D1SIG :
			case FOM_D2SIG :
			fprintf(fh, "%10.3f %10.2f %10i %10.2f "
			        "%10.3f  %10.3f\n",
			        cen*1.0e-9, r*100.0,
			        fom_shell_num_reflections(fctx, i),
			        (1.0/cen)*1e10,
			        shells->rmins[i]*1.0e-9,
			        shells->rmaxs[i]*1.0e-9);
			break;

			default :
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
	char *sym_str_fromfile = NULL;
	char *sym_str_fromfile1 = NULL;
	char *sym_str_fromfile2 = NULL;
	SymOpList *sym;
	RefList *list1_acc;
	RefList *list2_acc;
	RefList *list1;
	RefList *list2;
	RefList *list1_raw;
	RefList *list2_raw;
	enum fom_type fom = FOM_R1I;
	char *cellfile = NULL;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	double rmin, rmax;
	float sigma_cutoff = -INFINITY;
	int config_ignorenegs = 0;
	int config_zeronegs = 0;
	int config_unity = 0;
	int nshells = 10;
	char *shell_file = NULL;
	float highres, lowres;
	int mul_cutoff = 0;
	int anom;
	struct fom_rejections rej;

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
		{"min-measurements",   1, NULL,               11},
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

			case 10 :
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
			fom = fom_type_from_string(optarg);
			break;

			case 5 :
			if ( sscanf(optarg, "%f", &sigma_cutoff) != 1 ) {
				ERROR("Invalid value for --sigma-cutoff\n");
				return 1;
			}
			STATUS("WARNING: You are using --sigma-cutoff.  "
			       "Be aware that the figures of merit will not "
			       "reflect the entire data set!\n");
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

			case 11 :
			if ( sscanf(optarg, "%i", &mul_cutoff) != 1 ) {
				ERROR("Invalid value for --min-measurements\n");
				return 1;
			}
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
			return 1;

			case FOM_R2 :
			case FOM_R1I :
			case FOM_RSPLIT :
			case FOM_CC :
			case FOM_CCSTAR :
			case FOM_CCANO :
			case FOM_CRDANO :
			case FOM_RANO :
			case FOM_RANORSPLIT :
			case FOM_D1SIG :
			case FOM_D2SIG :
			break;

			default :
			ERROR("Unhandled figure of merit!\n");
			return 1;
		}
	}

	if ( (fom != FOM_R1F) && (config_ignorenegs || config_zeronegs) ) {
		ERROR("WARNING: You are using --zero-negs or --ignore-negs "
		      "even though your chosen figure of merit does not "
		      "require it.\n");
		ERROR("The figures of merit will not reflect the entire data "
		      "set!\n");
	}

	afile = strdup(argv[optind++]);
	bfile = strdup(argv[optind]);

	if ( shell_file == NULL ) shell_file = strdup("shells.dat");

	cell = load_cell_from_file(cellfile);
	if ( cellfile == NULL ) {
		ERROR("You must provide a unit cell.\n");
		exit(1);
	}
	if ( cell == NULL ) {
		ERROR("Failed to load cell.\n");
		return 1;
	}
	free(cellfile);

	list1_raw = read_reflections_2(afile, &sym_str_fromfile1);
	if ( list1_raw == NULL ) {
		ERROR("Couldn't read file '%s'\n", afile);
		return 1;
	}

	list2_raw = read_reflections_2(bfile, &sym_str_fromfile2);
	if ( list2_raw == NULL ) {
		ERROR("Couldn't read file '%s'\n", bfile);
		return 1;
	}

	if ( (sym_str_fromfile1 != NULL) && (sym_str_fromfile2 != NULL) ) {
		if ( strcmp(sym_str_fromfile1, sym_str_fromfile2) != 0 ) {
			ERROR("The symmetries of the two list do not match:\n");
			ERROR(" %s: %s\n", afile, sym_str_fromfile1);
			ERROR(" %s: %s\n", bfile, sym_str_fromfile2);
			return 1;
		}
		sym_str_fromfile = sym_str_fromfile1;
		free(sym_str_fromfile2);
	}

	if ( sym_str == NULL ) {
		if ( sym_str_fromfile != NULL ) {
			STATUS("Using symmetry from reflection files: %s\n",
			       sym_str_fromfile);
			sym_str = sym_str_fromfile;
		} else {
			sym_str = strdup("1");
		}
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
			case FOM_D1SIG :
			case FOM_D2SIG :
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

			default :
			ERROR("Unhandled figure of merit type!\n");
			return 1;
		}
	}


	/* Check that the intensities have the correct symmetry */
	if ( check_list_symmetry(list1_raw, sym) ) {
		ERROR("The first input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(sym));
		if ( cell_get_lattice_type(cell) == L_MONOCLINIC ) {
			ERROR("You may need to specify the unique axis in your "
			      "point group.  The default is unique axis c.\n");
			ERROR("See 'man crystfel' for more details.\n");
		}
		return 1;
	}
	if ( check_list_symmetry(list2_raw, sym) ) {
		ERROR("The second input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(sym));
		if ( cell_get_lattice_type(cell) == L_MONOCLINIC ) {
			ERROR("You may need to specify the unique axis in your "
			      "point group.  The default is unique axis c.\n");
			ERROR("See 'man crystfel' for more details.\n");
		}
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
	reflist_free(list1_raw);
	reflist_free(list2_raw);

	anom = ( (fom == FOM_CCANO) || (fom == FOM_CRDANO)
	        || (fom == FOM_RANO) || (fom == FOM_RANORSPLIT) );
	rej = fom_select_reflection_pairs(list1, list2, &list1_acc, &list2_acc,
	                                  cell, sym,
	                                  anom, rmin_fix, rmax_fix, sigma_cutoff,
	                                  config_ignorenegs, config_zeronegs,
	                                  mul_cutoff);
	reflist_free(list1);
	reflist_free(list2);

	gsl_set_error_handler_off();

	if ( rej.low_snr > 0 ) {
		STATUS("Discarded %i reflection pairs because either or both"
		       " versions had I/sigma(I) < %f.\n",
		       rej.low_snr, sigma_cutoff);
	}

	if ( rej.negative_deleted > 0 ) {
		STATUS("Discarded %i reflection pairs because either or both"
		       " versions had negative intensities.\n",
		       rej.negative_deleted);
	}

	if ( rej.negative_zeroed > 0 ) {
		STATUS("For %i reflection pairs, either or both versions had"
		       " negative intensities which were set to zero.\n",
		       rej.negative_zeroed);
	}

	if ( rej.few_measurements > 0 ) {
		STATUS("%i reflection pairs rejected because either or both"
		       " versions had too few measurements.\n",
		       rej.few_measurements);
	}

	if ( rej.outside_resolution_range > 0 ) {
		STATUS("%i reflection pairs rejected because either or both"
		       " versions were outside the resolution range.\n",
		       rej.outside_resolution_range);
	}

	if ( rej.no_bijvoet > 0 ) {
		STATUS("%i reflection pairs rejected because either or both"
		       " versions did not have Bijvoet partners.\n",
		       rej.no_bijvoet);
	}

	if ( rej.centric > 0 ) {
		STATUS("%i reflection pairs rejected because they were"
		       " centric.\n", rej.centric);
	}

	STATUS("%i reflection pairs accepted.\n", rej.common);

	resolution_limits(list1_acc, cell, &rmin, &rmax);
	resolution_limits(list2_acc, cell, &rmin, &rmax);
	STATUS("Accepted resolution range: %f to %f nm^-1"
	       " (%.2f to %.2f Angstroms).\n",
	       rmin/1e9, rmax/1e9, 1e10/rmin, 1e10/rmax);

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
	       nshells, shell_file, sym);

	free(shell_file);
	reflist_free(list1_acc);
	reflist_free(list2_acc);

	return 0;
}
