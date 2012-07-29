/*
 * check_hkl.c
 *
 * Characterise reflection lists
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
"  -y, --symmetry=<sym>       The symmetry of the input file.\n"
"  -p, --pdb=<filename>       PDB file to use.\n"
"      --rmin=<res>           Fix lower resolution limit for resolution shells. (m^-1).\n"
"      --rmax=<res>           Fix upper resolution limit for resolution shells. (m^-1).\n"
"      --sigma-cutoff=<n>     Discard reflections with I/sigma(I) < n.\n"
"\n");
}


static void plot_shells(RefList *list, UnitCell *cell, const SymOpList *sym,
                        double rmin_fix, double rmax_fix)
{
	double num[NBINS];
	int cts[NBINS];
	int possible[NBINS];
	unsigned int measurements[NBINS];
	unsigned int measured[NBINS];
	unsigned int snr_measured[NBINS];
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
	int nrefl = 0;
	int nmeastot = 0;
	int nout = 0;
	int nsilly = 0;
	Reflection *refl;
	RefListIterator *iter;
	RefList *counted;
	int hmax, kmax, lmax;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

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
		snr_measured[i] = 0;
		measurements[i] = 0;
		snr[i] = 0;
		var[i] = 0;
		mean[i] = 0;
	}

	resolution_limits(list, cell, &rmin, &rmax);
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
	counted = reflist_new();
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

		d = 2.0 * resolution(cell, h, k, l);

		bin = -1;
		for ( i=0; i<NBINS; i++ ) {
			if ( (d>rmins[i]) && (d<=rmaxs[i]) ) {
				bin = i;
				break;
			}
		}
		if ( bin == -1 ) continue;

		get_asymm(sym, h, k, l, &hs, &ks, &ls);
		if ( find_refl(counted, hs, ks, ls) != NULL ) continue;
		add_refl(counted, hs, ks, ls);

		possible[bin]++;

	}
	}
	}
	reflist_free(counted);

	/* Calculate means */
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double d, val, esd;
		int bin;
		int j;

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
		if ( bin == -1 ) continue;

		measured[bin]++;
		mean[bin] += get_intensity(refl);

		if ( !isfinite(val/esd) ) nsilly++;

	}

	for ( i=0; i<NBINS; i++ ) {
		mean[i] /= (double)measured[i];
	}

	/* Characterise the data set */
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
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

		/* measured[bin] was done earlier */
		measurements[bin] += get_redundancy(refl);

		if ( isfinite(val/esd) ) {
			snr[bin] += val / esd;
			snr_total += val / esd;
			snr_measured[bin]++;
		}

		nrefl++;
		nmeastot += get_redundancy(refl);

		var[bin] += pow(val-mean[bin], 2.0);

	}
	STATUS("overall <snr> = %f\n", snr_total/(double)nrefl);
	STATUS("%i measurements in total.\n", nmeastot);
	STATUS("%i reflections in total.\n", nrefl);

	if ( nout ) {
		STATUS("Warning; %i reflections outside resolution range.\n",
		       nout);
	}

	if ( nsilly ) {
		STATUS("Warning; %i reflections had infinite or invalid values"
		       " of I/sigma(I).\n", nsilly);
	}


	fprintf(fh, "1/d centre   # refs Possible  Compl       "
		    "Meas   Red   SNR    Std dev       Mean     d(A)\n");
	for ( i=0; i<NBINS; i++ ) {

		double cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;
		fprintf(fh, "%10.3f %8i %8i %6.2f %10i %5.1f"
		            " %5.2f %10.2f %10.2f %8.2f\n",
		        cen*1.0e-9,
		        measured[i],
		        possible[i],
		        100.0*(double)measured[i]/possible[i],
		        measurements[i],
		        (double)measurements[i]/measured[i],
		        snr[i]/(double)snr_measured[i],
		        sqrt(var[i]/measured[i]),
		        mean[i], (1.0/cen)*1e10);

	}

	fclose(fh);

	STATUS("Resolution shell information written to shells.dat.\n");
}


int main(int argc, char *argv[])
{
	int c;
	UnitCell *cell;
	char *file = NULL;
	char *sym_str = NULL;
	SymOpList *sym;
	RefList *raw_list;
	RefList *list;
	Reflection *refl;
	RefListIterator *iter;
	char *pdb = NULL;
	int rej = 0;
	float rmin_fix = -1.0;
	float rmax_fix = -1.0;
	float sigma_cutoff = -INFINITY;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},
		{"rmin",               1, NULL,                2},
		{"rmax",               1, NULL,                3},
		{"sigma-cutoff",       1, NULL,                4},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hy:p:", longopts, NULL)) != -1) {

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
			if ( sscanf(optarg, "%f", &sigma_cutoff) != 1 ) {
				ERROR("Invalid value for --sigma-cutoff\n");
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

	if ( sym_str == NULL ) {
		sym_str = strdup("1");
	}
	sym = get_pointgroup(sym_str);
	free(sym_str);

	file = strdup(argv[optind++]);

	if ( pdb == NULL ) {
		ERROR("You need to provide a PDB file containing"
		       " the unit cell.\n");
		return 1;
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
		      " have symmetry %s\n", symmetry_name(sym));
		return 1;
	}

	/* Reject some reflections */
	list = reflist_new();
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

		if ( val < sigma_cutoff * sig ) {
			rej++;
			ig = 1;
		}

		if ( ig ) continue;

		new = add_refl(list, h, k, l);
		copy_data(new, refl);

	}
	STATUS("Discarded %i reflections (out of %i) with I/sigma(I) < %f\n",
	       rej, num_reflections(raw_list), sigma_cutoff);

	plot_shells(list, cell, sym, rmin_fix, rmax_fix);

	return 0;
}
