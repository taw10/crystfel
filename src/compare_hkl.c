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
	double *out;
	UnitCell *cell;
	char *outfile = NULL;
	char *afile = NULL;
	char *bfile = NULL;
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
	while ((c = getopt_long(argc, argv, "ho:", longopts, NULL)) != -1) {

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
	if ( ref1 == NULL ) {
		ERROR("Couldn't open file '%s'\n", afile);
		return 1;
	}
	ref2 = new_list_intensity();
	i2 = read_reflections(bfile, ref2, NULL, NULL);
	if ( ref2 == NULL ) {
		ERROR("Couldn't open file '%s'\n", bfile);
		return 1;
	}

	/* Find common reflections */
	icommon = intersection_items(i1, i2);
	ncom = num_items(icommon);

	/* List for output scale factor map */
	out = new_list_intensity();

	for ( i=0; i<ncom; i++ ) {

		double i1, i2;
		struct refl_item *it;
		signed int h, k, l;

		it = get_item(icommon, i);
		h = it->h;  k = it->k;  l = it->l;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = lookup_intensity(ref2, h, k, l);

		set_intensity(out, h, k, l, i1/i2);

	}

	STATUS("%i,%i reflections: %i in common\n",
	       num_items(i1), num_items(i2), ncom);
	R2 = stat_r2(ref1, ref2, icommon, &scale);
	STATUS("R2 = %5.4f %% (scale=%5.2e)\n", R2*100.0, scale);
	Rmerge = stat_rmerge(ref1, ref2, icommon, &scale);
	STATUS("Rmerge = %5.4f %% (scale=%5.2e)\n", Rmerge*100.0, scale);
	pearson = stat_pearson(ref1, ref2, icommon);
	STATUS("Pearson r = %5.4f\n", pearson);

	if ( outfile != NULL ) {
		write_reflections(outfile, icommon, out, NULL, NULL, cell);
	}

	return 0;
}
