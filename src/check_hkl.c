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
#include "sfac.h"
#include "reflections.h"
#include "statistics.h"
#include "symmetry.h"


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


static void plot_shells(const double *ref, ReflItemList *items, UnitCell *cell,
                        const char *sym, unsigned int *counts,
                        const double *sigma, double rmin_fix, double rmax_fix)
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

	rmin = +INFINITY;
	rmax = 0.0;
	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double d;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

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

	/* Characterise the data set */
	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double d;
		int bin;
		int j;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

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
		measurements[bin] += lookup_count(counts, h, k, l);
		snr[bin] += (lookup_intensity(ref1, h, k, l) /
		                              lookup_intensity(sigma, h, k, l));
		snr_total += (lookup_intensity(ref1, h, k, l) /
		                              lookup_intensity(sigma, h, k, l));
		nmeas++;
		nmeastot += lookup_count(counts, h, k, l);

	}
	STATUS("overall <snr> = %f\n", snr_total/(double)nmeas);
	STATUS("%i measurements in total.\n", nmeastot);
	STATUS("%i reflections in total.\n", nmeas);

	if ( nout ) {
		STATUS("Warning; %i reflections outside resolution range.\n",
		       nout);
	}

	for ( i=0; i<NBINS; i++ ) {

		double r, cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;
		r = (num[i]/den)*((double)ctot/cts[i]);
		fprintf(fh, "%f %i %i %5.2f %i %f %f\n", cen*1.0e-9, measured[i],
		        possible[i], 100.0*(float)measured[i]/possible[i],
		        measurements[i], (float)measurements[i]/measured[i],
		        (snr[i]/(double)measured[i]));

	}

	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	double *ref;
	UnitCell *cell;
	char *file = NULL;
	char *sym = NULL;
	int i;
	ReflItemList *items;
	ReflItemList *good_items;
	char *pdb = NULL;
	double *esd;
	int rej = 0;
	unsigned int *cts;
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
		ERROR("Please provide exactly one HKL files to check.\n");
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

	ref = new_list_intensity();
	esd = new_list_sigma();
	cts = new_list_count();
	items = read_reflections(file, ref, NULL, cts, esd);
	if ( items == NULL ) {
		ERROR("Couldn't open file '%s'\n", file);
		return 1;
	}

	/* Reject reflections */
	good_items = new_items();
	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double val, sig;
		int ig = 0;
		double d;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		val = lookup_intensity(ref, h, k, l);
		sig = lookup_sigma(esd, h, k, l);

		if ( val < 3.0 * sig ) {
			rej++;
			ig = 1;
		}

		d = 0.5/resolution(cell, h, k, l);
		if ( d > 55.0e-10 ) ig = 1;
		//if ( d < 15.0e-10 ) ig = 1;

		//if ( ig ) continue;

		add_item(good_items, h, k, l);

	}

	plot_shells(ref, items, cell, sym, cts, esd, rmin_fix, rmax_fix);

	return 0;
}
