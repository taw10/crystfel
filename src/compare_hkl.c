/*
 * compare_hkl.c
 *
 * Compare reflection lists
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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

#include "utils.h"
#include "statistics.h"
#include "symmetry.h"
#include "reflist-utils.h"


/* Number of bins for plot of resolution shells */
#define NBINS (10)


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
"      --shells               Plot the figures of merit by resolution.\n"
"      --rmin=<res>           Fix lower resolution limit for --shells (m^-1).\n"
"      --rmax=<res>           Fix upper resolution limit for --shells (m^-1).\n"
"\n");
}


static void plot_shells(RefList *list1, RefList *list2, double scale,
                        UnitCell *cell, const char *sym,
                        double rmin_fix, double rmax_fix)
{
	double num[NBINS];
	int cts[NBINS];
	int possible[NBINS];
	unsigned int *counted;
	unsigned int measurements[NBINS];
	unsigned int measured[NBINS];
	double total_vol, vol_per_shell;
	double rmins[NBINS];
	double rmaxs[NBINS];
	double snr[NBINS];
	double den;
	double rmin, rmax;
	signed int h, k, l;
	int i;
	int ctot;
	FILE *fh;
	int nout = 0;
	Reflection *refl1;
	RefListIterator *iter;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	int hmax, kmax, lmax;

	if ( cell == NULL ) {
		ERROR("Need the unit cell to plot resolution shells.\n");
		return;
	}

	fh = fopen("shells.dat", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open 'shells.dat'\n");
		return;
	}

	for ( i=0; i<NBINS; i++ ) {
		num[i] = 0.0;
		cts[i] = 0;
		possible[i] = 0;
		measured[i] = 0;
		measurements[i] = 0;
		snr[i] = 0;
	}

	/* Iterate over all common reflections and calculate min and max
	 * resolution */
	rmin = +INFINITY;  rmax = 0.0;
	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) ) {

		signed int h, k, l;
		double d;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		d = resolution(cell, h, k, l) * 2.0;
		if ( d > rmax ) rmax = d;
		if ( d < rmin ) rmin = d;

	}

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

	/* Count the number of reflections possible in each shell */
	counted = new_list_count();
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	hmax = rmax / modulus(asx, asy, asz);
	kmax = rmax / modulus(bsx, bsy, bsz);
	lmax = rmax / modulus(csx, csy, csz);

	for ( h=-hmax; h<hmax; h++ ) {
	for ( k=-kmax; k<kmax; k++ ) {
	for ( l=-lmax; l<lmax; l++ ) {

		double d;
		signed int hs, ks, ls;
		int bin;

		get_asymm(h, k, l, &hs, &ks, &ls, sym);
		if ( lookup_count(counted, hs, ks, ls) ) continue;
		set_count(counted, hs, ks, ls, 1);

		d = resolution(cell, h, k, l) * 2.0;

		bin = -1;
		for ( i=0; i<NBINS; i++ ) {
			if ( (d>rmins[i]) && (d<=rmaxs[i]) ) {
				bin = i;
				break;
			}
		}
		if ( bin == -1 ) continue;

		possible[bin]++;

	}
	}
	}
	free(counted);

	den = 0.0;
	ctot = 0;
	nout = 0;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) ) {

		signed int h, k, l;
		double d;
		int bin;
		double i1, i2;
		int j;
		Reflection *refl2;

		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;  /* No common reflection */

		d = resolution(cell, h, k, l) * 2.0;

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

		i1 = get_intensity(refl1);
		i2 = get_intensity(refl2);
		i2 *= scale;

		num[bin] += fabs(i1 - i2);
		den += i1;
		ctot++;
		cts[bin]++;

	}

	if ( nout ) {
		STATUS("Warning; %i reflections outside resolution range.\n",
		       nout);
	}

	for ( i=0; i<NBINS; i++ ) {

		double r, cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;
		r = (num[i]/den)*((double)ctot/cts[i]);
		fprintf(fh, "%f %f\n", cen*1.0e-9, r*100.0);

	}

	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *ratiofile = NULL;
	char *afile = NULL;
	char *bfile = NULL;
	char *sym = NULL;
	double scale, scale_r2, scale_rdig, R1, R2, R1i, Rdiff, pearson;
	double scale_rintint, scale_r1i, scale_r1, scale_r1fi;
	int ncom;
	RefList *list1;
	RefList *list2;
	RefList *list2_transformed;
	RefList *ratio;
	RefList *deleteme;
	int config_shells = 0;
	char *pdb = NULL;
	int rej1 = 0;
	int rej2 = 0;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	Reflection *refl1;
	RefListIterator *iter;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"ratio" ,             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"shells",             0, &config_shells,     1},
		{"pdb",                1, NULL,               'p'},
		{"rmin",               1, NULL,               2},
		{"rmax",               1, NULL,               3},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:y:p:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'o' :
			ratiofile = strdup(optarg);
			break;

		case 'y' :
			sym = strdup(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
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

		default :
			return 1;
		}

	}

	if ( argc != (optind+2) ) {
		ERROR("Please provide exactly two HKL files to compare.\n");
		return 1;
	}

	if ( sym == NULL ) {
		sym = strdup("1");
	}

	afile = strdup(argv[optind++]);
	bfile = strdup(argv[optind]);

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

	cell = load_cell_from_pdb(pdb);
	free(pdb);

	list1 = read_reflections(afile);
	if ( list1 == NULL ) {
		ERROR("Couldn't read file '%s'\n", afile);
		return 1;
	}

	list2 = read_reflections(bfile);
	if ( list2 == NULL ) {
		ERROR("Couldn't read file '%s'\n", bfile);
		return 1;
	}

	/* Check that the intensities have the correct symmetry */
	if ( check_list_symmetry(list1, sym) ) {
		ERROR("The first input reflection list does not appear to"
		      " have symmetry %s\n", sym);
		return 1;
	}
	if ( check_list_symmetry(list2, sym) ) {
		ERROR("The second input reflection list does not appear to"
		      " have symmetry %s\n", sym);
		return 1;
	}

	/* Find common reflections (taking symmetry into account) */
	list2_transformed = reflist_new();
	ratio = reflist_new();
	deleteme = reflist_new();
	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) ) {

		signed int h, k, l;
		signed int he, ke, le;
		double val1, val2;
		double sig1, sig2;
		int ig = 0;
		double d;
		Reflection *refl2;
		Reflection *tr;

		get_indices(refl1, &h, &k, &l);

		if ( !find_equiv_in_list(list2, h, k, l, sym, &he, &ke, &le) ) {
			/* No common reflection */
			add_refl(deleteme, h, k, l);
			continue;
		}

		refl2 = find_refl(list2, he, ke, le);
		assert(refl2 != NULL);

		val1 = get_intensity(refl1);
		val2 = get_intensity(refl2);
		sig1 = get_esd_intensity(refl1);
		sig2 = get_esd_intensity(refl2);

		if ( val1 < 3.0 * sig1 ) {
			rej1++;
			ig = 1;
		}
		if ( val2 < 3.0 * sig2 ) {
			rej2++;
			ig = 1;
		}

		d = 0.5/resolution(cell, h, k, l);
		if ( d > 55.0e-10 ) ig = 1;
		//if ( d < 15.0e-10 ) ig = 1;

		//if ( ig ) continue;

		/* Add the old data from 'refl2' to a new list with the same
		 * indices as its equivalent in 'list1' */
		tr = add_refl(list2_transformed, h, k, l);
		copy_data(tr, refl2);

		/* Add divided version to 'output' list */
		tr = add_refl(ratio, h, k, l);
		set_int(tr, val1/val2);

	}

	STATUS("%i reflections in '%s' had I < 3.0*sigma(I)\n", rej1, afile);
	STATUS("%i reflections in '%s' had I < 3.0*sigma(I)\n", rej2, bfile);

	ncom = num_reflections(list2_transformed);
	STATUS("%i,%i reflections: %i in common\n",
	       num_reflections(list1), num_reflections(list2), ncom);

	/* Trim reflections from 'list1' which had no equivalents in 'list2' */
	for ( refl1 = first_refl(deleteme, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) ) {

		signed int h, k, l;
		Reflection *del;

		get_indices(refl1, &h, &k, &l);
		del = find_refl(list1, h, k, l);
		assert(del != NULL);

		delete_refl(del);

	}
	reflist_free(deleteme);
	reflist_free(list2);

	R1 = stat_r1_ignore(list1, list2_transformed, &scale_r1fi);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e)"
	       " (ignoring negative intensities)\n",
	       R1*100.0, scale_r1fi);

	R1 = stat_r1_zero(list1, list2_transformed, &scale_r1);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e)"
	       " (zeroing negative intensities)\n",
	       R1*100.0, scale_r1);

	R2 = stat_r2(list1, list2_transformed, &scale_r2);
	STATUS("R2(I) = %5.4f %% (scale=%5.2e)\n", R2*100.0, scale_r2);

	R1i = stat_r1_i(list1, list2_transformed, &scale_r1i);
	STATUS("R1(I) = %5.4f %% (scale=%5.2e)\n", R1i*100.0, scale_r1i);

	Rdiff = stat_rdiff_ignore(list1, list2_transformed, &scale_rdig);
	STATUS("Rint(F) = %5.4f %% (scale=%5.2e)"
	       " (ignoring negative intensities)\n",
	       Rdiff*100.0, scale_rdig);

	Rdiff = stat_rdiff_zero(list1, list2_transformed, &scale);
	STATUS("Rint(F) = %5.4f %% (scale=%5.2e)"
	       " (zeroing negative intensities)\n",
	       Rdiff*100.0, scale);

	Rdiff = stat_rdiff_intensity(list1, list2_transformed, &scale_rintint);
	STATUS("Rint(I) = %5.4f %% (scale=%5.2e)\n",
	       Rdiff*100.0, scale_rintint);

	pearson = stat_pearson_i(list1, list2_transformed);
	STATUS("Pearson r(I) = %5.4f\n", pearson);

	pearson = stat_pearson_f_ignore(list1, list2_transformed);
	STATUS("Pearson r(F) = %5.4f (ignoring negative intensities)\n",
	       pearson);

	pearson = stat_pearson_f_zero(list1, list2_transformed);
	STATUS("Pearson r(F) = %5.4f (zeroing negative intensities)\n",
	       pearson);

	if ( config_shells ) {
		plot_shells(list1, list2_transformed, scale_r1fi,
		            cell, sym, rmin_fix, rmax_fix);
	}

	if ( ratiofile != NULL ) {
		write_reflist(ratiofile, ratio, cell);
	}

	return 0;
}
