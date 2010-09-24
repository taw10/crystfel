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
                        ReflItemList *items, double scale, UnitCell *cell)
{
	double num[NBINS];
	double den[NBINS];
	double rmin, rmax, rstep;
	int i;
	FILE *fh;

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
		den[i] = 0.0;
	}

	rmin = +INFINITY;
	rmax = 0.0;
	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double res;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		res = 2.0*resolution(cell, h, k, l);
		if ( res > rmax ) rmax = res;
		if ( res < rmin ) rmin = res;

	}
	rstep = (rmax-rmin) / NBINS;

	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		double res;
		int bin;
		double i1, i2;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		res = 2.0*resolution(cell, h, k, l);

		bin = (res-rmin)/rstep;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = scale * lookup_intensity(ref2, h, k, l);

		num[bin] += pow(i1 - i2, 2.0);
		den[bin] += pow(i1, 2.0);

	}

	for ( i=0; i<NBINS; i++ ) {

		double r2, cen;
		cen = rmin + rstep*i + rstep/2.0;
		r2 = sqrt(num[i]/den[i]);
		fprintf(fh, "%f %f\n", cen/1.0e9, r2*100.0);

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
	double scale, scale_r2, R1, R2, R1i, Rdiff, pearson;
	int i, ncom;
	ReflItemList *i1, *i2, *icommon;
	int config_shells = 0;
	char *pdb = NULL;
	double *esd1;
	double *esd2;
	int rej1 = 0;
	int rej2 = 0;

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
	i1 = read_reflections(afile, ref1, NULL, NULL, esd1);
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
		if ( ig ) continue;

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

	R1 = stat_r1_ignore(ref1, ref2_transformed, icommon, &scale);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e) (ignoring negative intensities)\n",
	       R1*100.0, scale);

	R1 = stat_r1_zero(ref1, ref2_transformed, icommon, &scale);
	STATUS("R1(F) = %5.4f %% (scale=%5.2e) (zeroing negative intensities)\n",
	       R1*100.0, scale);

	R2 = stat_r2(ref1, ref2_transformed, icommon, &scale_r2);
	STATUS("R2(I) = %5.4f %% (scale=%5.2e)\n", R2*100.0, scale_r2);

	R1i = stat_r1_i(ref1, ref2_transformed, icommon, &scale);
	STATUS("R1(I) = %5.4f %% (scale=%5.2e)\n", R1i*100.0, scale);

	Rdiff = stat_rdiff_ignore(ref1, ref2_transformed, icommon, &scale);
	STATUS("Rdiff(F) = %5.4f %% (scale=%5.2e) (ignoring negative intensities)\n",
	       Rdiff*100.0, scale);

	Rdiff = stat_rdiff_zero(ref1, ref2_transformed, icommon, &scale);
	STATUS("Rdiff(F) = %5.4f %% (scale=%5.2e) (zeroing negative intensities)\n",
	       Rdiff*100.0, scale);

	pearson = stat_pearson_i(ref1, ref2_transformed, icommon);
	STATUS("Pearson r(I) = %5.4f\n", pearson);

	pearson = stat_pearson_f_ignore(ref1, ref2_transformed, icommon);
	STATUS("Pearson r(F) = %5.4f (ignoring negative intensities)\n",
	       pearson);

	pearson = stat_pearson_f_zero(ref1, ref2_transformed, icommon);
	STATUS("Pearson r(F) = %5.4f (zeroing negative intensities)\n",
	       pearson);

	if ( config_shells ) {
		plot_shells(ref1, ref2_transformed, icommon, scale_r2, cell);
	}

	if ( outfile != NULL ) {
		write_reflections(outfile, icommon, out, NULL, NULL, cell);
	}

	return 0;
}
