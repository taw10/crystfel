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


static void show_help(const char *s)
{
	printf("Syntax: %s [options] -a <file1.hkl> -b <file2.hkl>\n\n", s);
	printf(
"Compare intensity lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"  -o, --output=<filename>    Specify output filename for correction factor.\n"
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
	signed int h, k, l;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:a:b:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'o' : {
			outfile = strdup(optarg);
			break;
		}

		case 'a' : {
			afile = strdup(optarg);
			break;
		}

		case 'b' : {
			bfile = strdup(optarg);
			break;
		}

		case 0 : {
			break;
		}

		default : {
			return 1;
		}
		}

	}

	if ( outfile == NULL ) {
		ERROR("You must specify the output filename with -o\n");
		return 1;
	}

	cell = load_cell_from_pdb("molecule.pdb");
	ref1 = read_reflections(afile, NULL);
	if ( ref1 == NULL ) {
		ERROR("Couldn't open file '%s'\n", afile);
		return 1;
	}
	ref2 = read_reflections(bfile, NULL);
	if ( ref2 == NULL ) {
		ERROR("Couldn't open file '%s'\n", bfile);
		return 1;
	}
	out = new_list_intensity();

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {

		double i1, i2;

		i1 = lookup_intensity(ref1, h, k, l);
		i2 = lookup_intensity(ref2, h, k, l);

		if ( (i1 != 0.0) && (i2 != 0.0) ) {
			set_intensity(out, h, k, l, i1/i2);
		}

	}
	}
	}

	write_reflections(outfile, NULL, out, 1, cell);

	return 0;
}
