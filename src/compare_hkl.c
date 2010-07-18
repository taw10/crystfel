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


static void show_help(const char *s)
{
	printf("Syntax: %s [options] <file1.hkl> <file2.hkl>\n\n", s);
	printf(
"Compare intensity lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -o, --output=<filename>    Specify output filename for correction factor.\n"
"  -y, --symmetry=<sym>       The symmetry of both the input files.\n"
"\n");
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
	double scale, R2, Rmerge, pearson;
	int i, ncom;
	ReflItemList *i1, *i2, *icommon;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:y:", longopts, NULL)) != -1) {

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

	afile = strdup(argv[optind++]);
	bfile = strdup(argv[optind]);

	cell = load_cell_from_pdb("molecule.pdb");
	ref1 = new_list_intensity();
	i1 = read_reflections(afile, ref1, NULL, NULL);
	if ( i1 == NULL ) {
		ERROR("Couldn't open file '%s'\n", afile);
		return 1;
	}
	ref2 = new_list_intensity();
	i2 = read_reflections(bfile, ref2, NULL, NULL);
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

		it = get_item(i1, i);
		h = it->h;  k = it->k;  l = it->l;

		if ( !find_unique_equiv(i2, h, k, l, sym, &he, &ke, &le) ) {
			continue;
		}

		val1 = lookup_intensity(ref1, h, k, l);
		val2 = lookup_intensity(ref2, he, ke, le);
		set_intensity(ref2_transformed, h, k, l, val2);
		set_intensity(out, h, k, l, val1/val2);
		add_item(icommon, h, k, l);

	}
	ncom = num_items(icommon);

	STATUS("%i,%i reflections: %i in common\n",
	       num_items(i1), num_items(i2), ncom);
	R2 = stat_r2(ref1, ref2_transformed, icommon, &scale);
	STATUS("R2 = %5.4f %% (scale=%5.2e)\n", R2*100.0, scale);
	Rmerge = stat_rmerge(ref1, ref2_transformed, icommon, &scale);
	STATUS("Rmerge = %5.4f %% (scale=%5.2e)\n", Rmerge*100.0, scale);
	pearson = stat_pearson(ref1, ref2_transformed, icommon);
	STATUS("Pearson r = %5.4f\n", pearson);

	if ( outfile != NULL ) {
		write_reflections(outfile, icommon, out, NULL, NULL, cell);
	}

	return 0;
}
