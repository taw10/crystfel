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

#include "utils.h"
#include "sfac.h"
#include "reflections.h"
#include "statistics.h"
#include "symmetry.h"


/* Number of bins for plot of resolution shells */
#define NBINS (10)


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file1.hkl> <file2.hkl>\n\n", s);
	printf(
"Compare intensity lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -o, --output=<filename>    Specify output filename for correction factor.\n"
"  -y, --symmetry=<sym>       The symmetry of both the input files.\n"
"  -p, --pdb=<filename>       PDB file to use (default: molecule.pdb).\n"
"      --shells               Plot the figures of merit by resolution.\n"
"\n");
}


static void plot_shells(const double *ref1, const double *ref2,
                        ReflItemList *items, double scale, UnitCell *cell,
                        const char *sym, ReflItemList *characterise,
                        unsigned int *char_counts, const double *sigma)
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
	double snr_total = 0;
	int nmeas = 0;
	int nmeastot = 0;

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

	STATUS("%f -> %f\n", rmin/1e9, rmax/1e9);

	/* Increase the max just a little bit */
	rmax += 0.001e9;

	/* FIXME: Fixed resolution shells */
	rmin = 0.120e9;
	rmax = 1.172e9;

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

#if 0
	/* FIXME: Fixed resolution shells */
	rmins[0] = 0.121065;  rmaxs[0] = 0.552486;
	rmins[1] = 0.552486;  rmaxs[1] = 0.690186;
	rmins[2] = 0.690186;  rmaxs[2] = 0.787787;
	rmins[3] = 0.787787;  rmaxs[3] = 0.865813;
	rmins[4] = 0.865813;  rmaxs[4] = 0.931853;
	rmins[5] = 0.931853;  rmaxs[5] = 0.989663;
	rmins[6] = 0.989663;  rmaxs[6] = 1.041409;
	rmins[7] = 1.041409;  rmaxs[7] = 1.088467;
	rmins[8] = 1.088467;  rmaxs[8] = 1.131775;
	rmins[9] = 1.131775;  rmaxs[9] = 1.172000;
	for ( i=0; i<NBINS; i++ ) {
		rmins[i] *= 1e9;
		rmaxs[i] *= 1e9;
	}
#endif

	/* Count the number of reflections possible in each shell */
	counted = new_list_count();
	for ( h=-50; h<=+50; h++ ) {
	for ( k=-50; k<=+50; k++ ) {
	for ( l=-50; l<=+50; l++ ) {

		double d;
		signed int hs, ks, ls;
		int bin;

		/* FIXME: Reflection condition */
		if ( (h==0) && (k==0) && (l%2) ) continue;

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

	/* Characterise the first data set (only) */
	for ( i=0; i<num_items(characterise); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double d;
		int bin;
		int j;

		it = get_item(characterise, i);
		h = it->h;  k = it->k;  l = it->l;

		/* FIXME: Reflection condition */
		if ( (h==0) && (k==0) && (l%2) ) continue;

		d = resolution(cell, h, k, l) * 2.0;

		bin = -1;
		for ( j=0; j<NBINS; j++ ) {
			if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
				bin = j;
				break;
			}
		}
		if ( bin == -1 ) {
			ERROR("Warnung! %i %i %i %f\n", h, k, l, d/1e9);
			continue;
		}

		measured[bin]++;
		measurements[bin] += lookup_count(char_counts, h, k, l);
		snr[bin] += (lookup_intensity(ref1, h, k, l) /
		                              lookup_intensity(sigma, h, k, l));
		snr_total += (lookup_intensity(ref1, h, k, l) /
		                              lookup_intensity(sigma, h, k, l));
		nmeas++;
		nmeastot += lookup_count(char_counts, h, k, l);

	}
	STATUS("overall <snr> = %f\n", snr_total/(double)nmeas);
	STATUS("%i measurements in total.\n", nmeastot);
	STATUS("%i reflections in total.\n", nmeas);

	den = 0.0;
	ctot = 0;
	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double d;
		int bin;
		double i1, i2, f1, f2;
		int j;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		/* FIXME: Reflection condition */
		if ( (h==0) && (k==0) && (l%2) ) continue;

		d = resolution(cell, h, k, l) * 2.0;

		bin = -1;
		for ( j=0; j<NBINS; j++ ) {
			if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
				bin = j;
				break;
			}
		}
		if ( bin == -1 ) {
			ERROR("Warnung! %i %i %i %f\n", h, k, l, d/1e9);
			abort();
			continue;
		}

		i1 = lookup_intensity(ref1, h, k, l);
		//if ( i1 < 0.0 ) continue;
		//f1 = sqrt(i1);
		i2 = lookup_intensity(ref2, h, k, l);
		//if ( i2 < 0.0 ) continue;
		//f2 = sqrt(i2);
		i2 *= scale;

		num[bin] += fabs(i1 - i2);
		den += i1;// + i2) / 2.0;
		ctot++;
		cts[bin]++;

	}

	for ( i=0; i<NBINS; i++ ) {

		double r, cen;
		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;
		r = (num[i]/den)*((double)ctot/cts[i]);
		fprintf(fh, "%f %f %i %i %i %f %f\n", cen*1.0e-9, r*100.0,
		                            measured[i],
		                            possible[i], measurements[i],
		                            (float)measurements[i]/measured[i],
		                            (snr[i]/(double)measured[i]));

	}

	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	double *ref1;
	double *ref2;
	double *ref2_transformed;
	double *out;
	UnitCell *cell;
	char *outfile = NULL;
	char *afile = NULL;
	char *bfile = NULL;
	char *sym = NULL;
	double scale, scale_r2, scale_rdig, R1, R2, R1i, Rdiff, pearson;
	double scale_rintint, scale_r1i, scale_r1, scale_r1fi;
	int i, ncom;
	ReflItemList *i1, *i2, *icommon;
	int config_shells = 0;
	char *pdb = NULL;
	double *esd1;
	double *esd2;
	int rej1 = 0;
	int rej2 = 0;
	unsigned int *cts1;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"shells",             0, &config_shells,     1},
		{"pdb",                1, NULL,               'p'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:y:p:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'o' :
			outfile = strdup(optarg);
			break;

		case 'y' :
			sym = strdup(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
			break;

		case 0 :
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

	ref1 = new_list_intensity();
	esd1 = new_list_sigma();
	cts1 = new_list_count();
	i1 = read_reflections(afile, ref1, NULL, cts1, esd1);
	if ( i1 == NULL ) {
		ERROR("Couldn't open file '%s'\n", afile);
		return 1;
	}
	ref2 = new_list_intensity();
	esd2 = new_list_sigma();
	i2 = read_reflections(bfile, ref2, NULL, NULL, esd2);
	if ( i2 == NULL ) {
		ERROR("Couldn't open file '%s'\n", bfile);
		return 1;
	}

	/* List for output scale factor map */
	out = new_list_intensity();

	/* Find common reflections (taking symmetry into account) */
	icommon = new_items();
	ref2_transformed = new_list_intensity();
	for ( i=0; i<num_items(i1); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		signed int he, ke, le;
		double val1, val2;
		double sig1, sig2;
		int ig = 0;
		double d;

		it = get_item(i1, i);
		h = it->h;  k = it->k;  l = it->l;

		if ( !find_unique_equiv(i2, h, k, l, sym, &he, &ke, &le) ) {
			//STATUS("%i %i %i not matched (%f nm).\n", h, k, l,
			//       1.0/(2.0*resolution(cell, h, k, l)/1e9));
			continue;
		}

		val1 = lookup_intensity(ref1, h, k, l);
		val2 = lookup_intensity(ref2, he, ke, le);
		sig1 = lookup_sigma(esd1, h, k, l);
		sig2 = lookup_sigma(esd2, he, ke, le);

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

		set_intensity(ref2_transformed, h, k, l, val2);
		set_intensity(out, h, k, l, val1/val2);
		add_item(icommon, h, k, l);

	}
	ncom = num_items(icommon);
	STATUS("%i reflections with I < 3.0*sigma(I) rejected from '%s'\n",
	       rej1, afile);
	STATUS("%i reflections with I < 3.0*sigma(I) rejected from '%s'\n",
	       rej2, bfile);

	STATUS("%i,%i reflections: %i in common\n",
	       num_items(i1), num_items(i2), ncom);

	R1 = stat_r1_ignore(ref1, ref2_transformed, icommon, &scale_r1fi);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e) (ignoring negative intensities)\n",
	       R1*100.0, scale_r1fi);

	R1 = stat_r1_zero(ref1, ref2_transformed, icommon, &scale_r1);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e) (zeroing negative intensities)\n",
	       R1*100.0, scale_r1);

	R2 = stat_r2(ref1, ref2_transformed, icommon, &scale_r2);
	STATUS("R2(I) = %5.4f %% (scale=%5.2e)\n", R2*100.0, scale_r2);

	R1i = stat_r1_i(ref1, ref2_transformed, icommon, &scale_r1i);
	STATUS("R1(I) = %5.4f %% (scale=%5.2e)\n", R1i*100.0, scale_r1i);

	Rdiff = stat_rdiff_ignore(ref1, ref2_transformed, icommon, &scale_rdig);
	STATUS("Rint(F) = %5.4f %% (scale=%5.2e) (ignoring negative intensities)\n",
	       Rdiff*100.0, scale_rdig);

	Rdiff = stat_rdiff_zero(ref1, ref2_transformed, icommon, &scale);
	STATUS("Rint(F) = %5.4f %% (scale=%5.2e) (zeroing negative intensities)\n",
	       Rdiff*100.0, scale);

	Rdiff = stat_rdiff_intensity(ref1, ref2_transformed, icommon,
	                             &scale_rintint);
	STATUS("Rint(I) = %5.4f %% (scale=%5.2e)\n",
	       Rdiff*100.0, scale_rintint);

	pearson = stat_pearson_i(ref1, ref2_transformed, icommon);
	STATUS("Pearson r(I) = %5.4f\n", pearson);

	pearson = stat_pearson_f_ignore(ref1, ref2_transformed, icommon);
	STATUS("Pearson r(F) = %5.4f (ignoring negative intensities)\n",
	       pearson);

	pearson = stat_pearson_f_zero(ref1, ref2_transformed, icommon);
	STATUS("Pearson r(F) = %5.4f (zeroing negative intensities)\n",
	       pearson);

	if ( config_shells ) {
		plot_shells(ref1, ref2_transformed, icommon, scale_r1fi,
		            cell, sym, i1, cts1, esd1);
	}

	if ( outfile != NULL ) {
		write_reflections(outfile, icommon, out, NULL, NULL, cell);
	}

	return 0;
}
