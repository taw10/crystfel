/*
 * check_hkl.c
 *
 * Characterise reflection lists
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

#include "utils.h"
#include "statistics.h"
#include "symmetry.h"
#include "reflist.h"
#include "reflist-utils.h"


/* Number of bins for plot of resolution shells */
#define NBINS (10)


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file.hkl>\n\n", s);
	printf(
"Characterise an intensity list.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -y, --symmetry=<sym>       The symmetry of both the input files.\n"
"  -p, --pdb=<filename>       PDB file to use (default: molecule.pdb).\n"
"      --rmin=<res>           Fix lower resolution limit for --shells (m^-1).\n"
"      --rmax=<res>           Fix upper resolution limit for --shells (m^-1).\n"
"\n");
}


static void plot_shells(RefList *list, UnitCell *cell, const char *sym,
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
	double mean[NBINS];
	double var[NBINS];
	double rmin, rmax;
	signed int h, k, l;
	int i;
	FILE *fh;
	double snr_total = 0;
	int nmeas = 0;
	int nmeastot = 0;
	int nout = 0;
	Reflection *refl;
	RefListIterator *iter;

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
		var[i] = 0;
		mean[i] = 0;
	}

	/* Iterate over all common reflections and calculate min and max
	 * resolution */
	rmin = +INFINITY;  rmax = 0.0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double d;

		get_indices(refl, &h, &k, &l);

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
	for ( h=-150; h<=+150; h++ ) {
	for ( k=-150; k<=+150; k++ ) {
	for ( l=-150; l<=+150; l++ ) {

		double d;
		signed int hs, ks, ls;
		int bin;

		d = resolution(cell, h, k, l) * 2.0;

		bin = -1;
		for ( i=0; i<NBINS; i++ ) {
			if ( (d>rmins[i]) && (d<=rmaxs[i]) ) {
				bin = i;
				break;
			}
		}
		if ( bin == -1 ) continue;

		get_asymm(h, k, l, &hs, &ks, &ls, sym);
		if ( lookup_count(counted, hs, ks, ls) ) continue;
		set_count(counted, hs, ks, ls, 1);

		possible[bin]++;

	}
	}
	}
	free(counted);

	/* Calculate means */
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double d;
		int bin;
		int j;

		get_indices(refl, &h, &k, &l);

		d = resolution(cell, h, k, l) * 2.0;

		bin = -1;
		for ( j=0; j<NBINS; j++ ) {
			if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
				bin = j;
				break;
			}
		}
		if ( bin == -1 ) {
			nout++;
			continue;
		}

		measured[bin]++;
		mean[bin] += get_intensity(refl);

	}

	for ( i=0; i<NBINS; i++ ) {
		mean[i] /= (double)measured[i];
	}

	/* Characterise the data set */
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double d;
		int bin;
		int j;
		double val, esd;

		get_indices(refl, &h, &k, &l);

		d = resolution(cell, h, k, l) * 2.0;
		val = get_intensity(refl);
		esd = get_esd_intensity(refl);

		bin = -1;
		for ( j=0; j<NBINS; j++ ) {
			if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
				bin = j;
				break;
			}
		}
		if ( bin == -1 ) {
			nout++;
			continue;
		}

		if ( !isfinite(val/esd) ) continue;

		/* measured[bin] was done earlier */
		measurements[bin] += get_redundancy(refl);
		snr[bin] += val / esd;
		snr_total += val / esd;
		nmeas++;
		nmeastot += get_redundancy(refl);

		var[bin] += pow(val-mean[bin], 2.0);

	}
	STATUS("overall <snr> = %f\n", snr_total/(double)nmeas);
	STATUS("%i measurements in total.\n", nmeastot);
	STATUS("%i reflections in total.\n", nmeas);

	if ( nout ) {
		STATUS("Warning; %i reflections outside resolution range.\n",
		       nout);
	}

	fprintf(fh, "1/d centre   # refs Possible  Compl       "
		    "Meas   Red   SNR    Std dev       Mean\n");
	for ( i=0; i<NBINS; i++ ) {

		double cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;
		fprintf(fh, "%10.3f %8i %8i %6.2f %10i %5.1f %5.2f %10.2f %10.2f\n",
		        cen*1.0e-9,
		        measured[i],
		        possible[i],
		        100.0*(float)measured[i]/possible[i],
		        measurements[i],
		        (float)measurements[i]/measured[i],
		        snr[i]/(double)measured[i],
		        sqrt(var[i]/measured[i]),
		        mean[i]);

	}

	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *file = NULL;
	char *sym = NULL;
	RefList *raw_list;
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	char *pdb = NULL;
	int rej = 0;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},
		{"rmin",               1, NULL,               2},
		{"rmax",               1, NULL,               3},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hy:p:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

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

	if ( argc != (optind+1) ) {
		ERROR("Please provide exactly one HKL file to check.\n");
		return 1;
	}

	if ( sym == NULL ) {
		sym = strdup("1");
	}

	file = strdup(argv[optind++]);

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}
	cell = load_cell_from_pdb(pdb);
	free(pdb);

	raw_list = read_reflections(file);
	if ( raw_list == NULL ) {
		ERROR("Couldn't read file '%s'\n", file);
		return 1;
	}

	/* Check that the intensities have the correct symmetry */
	if ( check_list_symmetry(raw_list, sym) ) {
		ERROR("The input reflection list does not appear to"
		      " have symmetry %s\n", sym);
		return 1;
	}

	/* Reject some reflections */
	list = reflist_new();
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double val, sig;
		int ig = 0;
		double d;
		Reflection *new;

		get_indices(refl, &h, &k, &l);

		val = get_intensity(refl);
		sig = get_esd_intensity(refl);

		if ( val < 3.0 * sig ) {
			rej++;
			ig = 1;
		}

		d = 0.5/resolution(cell, h, k, l);
		if ( d > 55.0e-10 ) ig = 1;
		//if ( d < 15.0e-10 ) ig = 1;

		//if ( ig ) continue;

		new = add_refl(list, h, k, l);
		copy_data(new, refl);

	}

	plot_shells(list, cell, sym, rmin_fix, rmax_fix);

	return 0;
}
