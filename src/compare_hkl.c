/*
 * compare_hkl.c
 *
 * Compare reflection lists
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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

#include "utils.h"
#include "statistics.h"
#include "symmetry.h"
#include "reflist-utils.h"


/* Number of bins for plot of resolution shells */
#define NBINS (10)


enum r_shell
{
	R_SHELL_NONE,
	R_SHELL_R1I,
	R_SHELL_R1F,
	R_SHELL_RSPLIT,
};


static enum r_shell get_r_shell(const char *s)
{
	if ( strcmp(s, "r1i") == 0 ) return R_SHELL_R1I;
	if ( strcmp(s, "r1f") == 0 ) return R_SHELL_R1F;
	if ( strcmp(s, "rsplit") == 0 ) return R_SHELL_RSPLIT;

	ERROR("Unknown R-factor '%s' - try '--shells=rsplit', or --help for"
	      " more possibilities.\n", s);
	exit(1);
}


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file1.hkl> <file2.hkl>\n\n", s);
	printf(
"Compare intensity lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -o, --ratio=<filename>     Specify output filename for ratios.\n"
"  -y, --symmetry=<sym>       The symmetry of both the input files.\n"
"  -p, --pdb=<filename>       PDB file to use (default: molecule.pdb).\n"
"      --shells=<FoM>         Plot this figure of merit in resolution shells.\n"
"                              Choose from: 'Rsplit', 'R1f' and 'R1i'.\n"
"      --rmin=<res>           Fix lower resolution limit for --shells (m^-1).\n"
"      --rmax=<res>           Fix upper resolution limit for --shells (m^-1).\n"
"\n");
}


static void plot_shells(RefList *list1, RefList *list2, double scale,
                        UnitCell *cell, double rmin_fix, double rmax_fix,
                        enum r_shell config_shells)
{
	double num[NBINS];
	int cts[NBINS];
	unsigned int measurements[NBINS];
	unsigned int measured[NBINS];
	double total_vol, vol_per_shell;
	double rmins[NBINS];
	double rmaxs[NBINS];
	double snr[NBINS];
	double rmin, rmax;
	int i;
	Reflection *refl1;
	RefListIterator *iter;
	FILE *fh;
	double den;
	int ctot, nout;

	if ( cell == NULL ) {
		ERROR("Need the unit cell to plot resolution shells.\n");
		return;
	}

	for ( i=0; i<NBINS; i++ ) {
		num[i] = 0.0;
		cts[i] = 0;
		measured[i] = 0;
		measurements[i] = 0;
		snr[i] = 0;
	}

	/* Find resolution limits */
	resolution_limits(list1, cell, &rmin, &rmax);
	STATUS("1/d goes from %f to %f nm^-1\n", rmin/1e9, rmax/1e9);

	/* Widen the range just a little bit */
	rmin -= 0.001e9;
	rmax += 0.001e9;

	/* Fixed resolution shells if needed */
	if ( rmin_fix > 0.0 ) rmin = rmin_fix;
	if ( rmax_fix > 0.0 ) rmax = rmax_fix;

	/* Calculate the resolution bins */
	total_vol = pow(rmax, 3.0) - pow(rmin, 3.0);
	vol_per_shell = total_vol / NBINS;
	rmins[0] = rmin;
	for ( i=1; i<NBINS; i++ ) {

		double r;

		r = vol_per_shell + pow(rmins[i-1], 3.0);
		r = pow(r, 1.0/3.0);

		/* Shells of constant volume */
		rmaxs[i-1] = r;
		rmins[i] = r;

		/* Shells of constant thickness */
		//rmins[i] = rmins[i-1] + (rmax-rmin)/NBINS;
		//rmaxs[i-1] = rmins[i-1] + (rmax-rmin)/NBINS;

		STATUS("Shell %i: %f to %f\n", i-1,
		       rmins[i-1]/1e9, rmaxs[i-1]/1e9);

	}
	rmaxs[NBINS-1] = rmax;
	STATUS("Shell %i: %f to %f\n", NBINS-1,
	       rmins[NBINS-1]/1e9, rmaxs[NBINS-1]/1e9);

	den = 0.0;  ctot = 0;  nout = 0;
	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		double d;
		int bin;
		double i1, i2, f1, f2;
		int j;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		d = 2.0 * resolution(cell, h, k, l);

		bin = -1;
		for ( j=0; j<NBINS; j++ ) {
			if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
				bin = j;
				break;
			}
		}

		/* Outside resolution range? */
		if ( bin == -1 ) {
			nout++;
			continue;
		}

		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		f1 = i1 > 0.0 ? sqrt(i1) : 0.0;
		f2 = i2 > 0.0 ? sqrt(i2) : 0.0;

		switch ( config_shells ) {

			case R_SHELL_RSPLIT :
			num[bin] += fabs(i1 - scale*i2);
			den += i1 + scale*i2;
			break;

			case R_SHELL_R1I :
			num[bin] += fabs(i1 - scale*i2);
			den += i1;
			break;

			case R_SHELL_R1F :
			num[bin] += fabs(f1 - scale*f2);
			den += f1;
			break;

			default : break;

		}

		ctot++;
		cts[bin]++;
	}

	if ( nout ) {
		STATUS("Warning; %i reflections outside resolution range.\n",
		       nout);
	}

	fh = fopen("shells.dat", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open 'shells.dat'\n");
		return;
	}

	switch ( config_shells ) {

		case R_SHELL_RSPLIT :
		fprintf(fh, "1/d centre   Rsplit / %%\n");
		break;

		case R_SHELL_R1I :
		fprintf(fh, "1/d centre   R1(I) / %%\n");
		break;

		case R_SHELL_R1F :
		fprintf(fh, "1/d centre   R1(F) ignoring -ves / %%\n");
		break;

		default :
		fprintf(fh, "1/d centre   0.0\n");
		break;

	}

	for ( i=0; i<NBINS; i++ ) {

		double r, cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;

		switch ( config_shells ) {

			case R_SHELL_RSPLIT :
			r = (2.0*(num[i]/den)*((double)ctot/cts[i]))/sqrt(2.0);
			break;

			case R_SHELL_R1I :
			case R_SHELL_R1F :
			r = (num[i]/den) * ((double)ctot/cts[i]);
			break;

			default :
			r = 0.0;
			break;

		}

		fprintf(fh, "%10.3f %10.2f %10i\n",
		        cen*1.0e-9, r*100.0, cts[i]);

	}

	fclose(fh);

	STATUS("Resolution shell information written to shells.dat.\n");
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *ratiofile = NULL;
	char *afile = NULL;
	char *bfile = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	double scale, scale_r2, scale_rdig, R1, R2, R1i, Rdiff, pearson;
	double scale_rintint, scale_r1i, scale_r1, scale_r1fi;
	int ncom;
	RefList *list1;
	RefList *list2;
	RefList *list1_raw;
	RefList *list2_raw;
	RefList *ratio;
	enum r_shell config_shells = R_SHELL_NONE;
	char *pdb = NULL;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	Reflection *refl1;
	RefListIterator *iter;
	int config_unity = 0;
	double scale_for_shells;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"ratio" ,             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"shells",             1, NULL,               4},
		{"pdb",                1, NULL,               'p'},
		{"rmin",               1, NULL,               2},
		{"rmax",               1, NULL,               3},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:y:p:u",
	                        longopts, NULL)) != -1)
	{

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'o' :
			ratiofile = strdup(optarg);
			break;

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
			config_shells = get_r_shell(optarg);
			break;

		default :
			return 1;
		}

	}

	if ( argc != (optind+2) ) {
		ERROR("Please provide exactly two HKL files to compare.\n");
		return 1;
	}

	if ( sym_str == NULL ) {
		sym_str = strdup("1");
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	afile = strdup(argv[optind++]);
	bfile = strdup(argv[optind]);

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

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

	list1 = asymmetric_indices(list1_raw, sym);
	list2 = asymmetric_indices(list2_raw, sym);

	/* Find common reflections and calculate ratio */
	ratio = reflist_new();
	ncom = 0;
	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		double val1, val2;
		Reflection *refl2;
		Reflection *tr;

		get_indices(refl1, &h, &k, &l);

		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;

		ncom++;

		val1 = get_intensity(refl1);
		val2 = get_intensity(refl2);

		/* Add divided version to 'output' list */
		tr = add_refl(ratio, h, k, l);
		set_int(tr, val1/val2);
		set_redundancy(tr, 1);
	}

	if ( ratiofile != NULL ) {
		write_reflist(ratiofile, ratio, cell);
	}
	reflist_free(ratio);

	gsl_set_error_handler_off();

	STATUS("%i,%i reflections: %i in common\n",
	       num_reflections(list1), num_reflections(list2), ncom);

	R1 = stat_r1_ignore(list1, list2, &scale_r1fi, config_unity);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e)"
	       " (ignoring negative intensities)\n",
	       R1*100.0, scale_r1fi);

	R1 = stat_r1_zero(list1, list2, &scale_r1, config_unity);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e)"
	       " (zeroing negative intensities)\n",
	       R1*100.0, scale_r1);

	R2 = stat_r2(list1, list2, &scale_r2, config_unity);
	STATUS("R2(I) = %5.4f %% (scale=%5.2e)\n", R2*100.0, scale_r2);

	R1i = stat_r1_i(list1, list2, &scale_r1i, config_unity);
	STATUS("R1(I) = %5.4f %% (scale=%5.2e)\n", R1i*100.0, scale_r1i);

	Rdiff = stat_rdiff_ignore(list1, list2, &scale_rdig,
	                          config_unity);
	STATUS("Rint(F) = %5.4f %% (scale=%5.2e)"
	       " (ignoring negative intensities)\n",
	       Rdiff*100.0, scale_rdig);

	Rdiff = stat_rdiff_zero(list1, list2, &scale, config_unity);
	STATUS("Rint(F) = %5.4f %% (scale=%5.2e)"
	       " (zeroing negative intensities)\n",
	       Rdiff*100.0, scale);

	Rdiff = stat_rdiff_intensity(list1, list2, &scale_rintint,
	                             config_unity);
	STATUS("Rint(I) = %5.4f %% (scale=%5.2e)\n",
	       Rdiff*100.0, scale_rintint);

	pearson = stat_pearson_i(list1, list2);
	STATUS("Pearson r(I) = %5.4f\n", pearson);

	pearson = stat_pearson_f_ignore(list1, list2);
	STATUS("Pearson r(F) = %5.4f (ignoring negative intensities)\n",
	       pearson);

	pearson = stat_pearson_f_zero(list1, list2);
	STATUS("Pearson r(F) = %5.4f (zeroing negative intensities)\n",
	       pearson);

	switch ( config_shells ) {
		case R_SHELL_R1I : scale_for_shells = scale_r1i;  break;
		case R_SHELL_R1F : scale_for_shells = scale_r1;  break;
		case R_SHELL_RSPLIT : scale_for_shells = scale_rintint;  break;
		default : scale_for_shells = 0.0;
	}

	if ( config_shells != R_SHELL_NONE ) {
		plot_shells(list1, list2, scale_for_shells,
		            cell, rmin_fix, rmax_fix, config_shells);
	}

	return 0;
}
