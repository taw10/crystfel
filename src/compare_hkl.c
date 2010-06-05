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
	double scale, R;
	unsigned int *c1;
	unsigned int *c2;
	int i;
	int nc1, nc2, ncom;

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

	cell = load_cell_from_pdb("molecule.pdb");
	c1 = new_list_count();
	ref1 = read_reflections(afile, c1, NULL);
	if ( ref1 == NULL ) {
		ERROR("Couldn't open file '%s'\n", afile);
		return 1;
	}
	c2 = new_list_count();
	ref2 = read_reflections(bfile, c2, NULL);
	if ( ref2 == NULL ) {
		ERROR("Couldn't open file '%s'\n", bfile);
		return 1;
	}
	out = new_list_intensity();

	set_count(c1, 0, 0, 0, 0);
	set_count(c2, 0, 0, 0, 0);

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

	nc1 = 0;
	nc2 = 0;
	ncom = 0;
	for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
		nc1 += c1[i];
		nc2 += c2[i];
		ncom += c1[i] && c2[i];
	}
	STATUS("%i,%i reflections: %i in common\n", nc1, nc2, ncom);
	R = stat_r2(ref1, c1, ref2, c2, &scale);
	STATUS("R2 = %5.4f %% (scale=%5.2f)\n", R*100.0, scale);

	if ( outfile != NULL ) {
		write_reflections(outfile, NULL, out, NULL, 1, cell, 1);
	}

	return 0;
}
