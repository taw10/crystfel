/*
 * get_hkl.c
 *
 * Small program to write out a list of h,k,l,I values given a structure
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
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Write idealised intensity lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -t, --template=<filename>  Only include reflections mentioned in file.\n"
"      --poisson              Simulate Poisson samples.\n"
"      --twin                 Generate twinned data.\n"
"  -o  --output=<filename>    Output filename (default: stdout).\n"
"      --zone-axis            Generate hk0 intensities only (and add\n"
"                              Synth2D-style header.\n"
);
}


static int template_reflections(const char *filename, unsigned int *counts)
{
	char *rval;
	FILE *fh;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		return 1;
	}

	do {

		char line[1024];
		int r;
		signed int h, k, l;

		rval = fgets(line, 1023, fh);

		r = sscanf(line, "%i %i %i", &h, &k, &l);
		if ( r != 3 ) continue;

		set_count(counts, h, k, l, 1);

	} while ( rval != NULL );

	fclose(fh);

	return 0;
}


/* Apply Poisson noise to all reflections */
static void noisify_reflections(double *ref)
{
	signed int h, k, l;

	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {

		double val;
		int c;

		val = lookup_intensity(ref, h, k, l);
		c = poisson_noise(val);
		set_intensity(ref, h, k, l, c);

	}
	}
	progress_bar(h+INDMAX, 2*INDMAX, "Simulating noise");
	}
}


int main(int argc, char *argv[])
{
	int c;
	double *ideal_ref;
	struct molecule *mol;
	char *template = NULL;
	int config_noisify = 0;
	int config_twin = 0;
	int config_za = 0;
	char *output = NULL;
	unsigned int *counts;
	signed int h, k, l;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"template",           1, NULL,               't'},
		{"poisson",            0, &config_noisify,     1},
		{"output",             1, NULL,               'o'},
		{"twin",               0, &config_twin,        1},
		{"zone-axis",          0, &config_za,          1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ht:o:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 't' : {
			template = strdup(optarg);
			break;
		}

		case 'o' : {
			output = strdup(optarg);
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

	mol = load_molecule();
	ideal_ref = get_reflections(mol, eV_to_J(1790.0), 1/(0.6e-9));

	counts = new_list_count();

	if ( template != NULL ) {

		if ( template_reflections(template, counts) != 0 ) {
			ERROR("Failed to template reflections.\n");
			return 1;
		}

	} else {

		for ( h=-INDMAX; h<=INDMAX; h++ ) {
		for ( k=-INDMAX; k<=INDMAX; k++ ) {
		for ( l=-INDMAX; l<=INDMAX; l++ ) {
			set_count(counts, h, k, l, 1);
		}
		}
		}

	}

	if ( config_noisify ) noisify_reflections(ideal_ref);

	if ( config_twin ) {

		for ( h=-INDMAX; h<=INDMAX; h++ ) {
		for ( k=-INDMAX; k<=INDMAX; k++ ) {
		for ( l=-INDMAX; l<=INDMAX; l++ ) {

			double a, b, c, d;
			double t;

			if ( abs(h+k) > INDMAX ) {
				set_intensity(ideal_ref, h, k, l, 0.0);
				continue;
			}

			a = lookup_intensity(ideal_ref, h, k, l);
			b = lookup_intensity(ideal_ref, k, h, -l);
			c = lookup_intensity(ideal_ref, -(h+k), k, -l);
			d = lookup_intensity(ideal_ref, -(h+k), h, l);

			t = (a+b+c+d)/4.0;

			set_intensity(ideal_ref, h, k, l, t);
			set_intensity(ideal_ref, k, h, -l, t);
			set_intensity(ideal_ref, -(h+k), h, l, t);
			set_intensity(ideal_ref, -(h+k), k, -l, t);

		}
		}
		}

	}

	write_reflections(output, counts, ideal_ref, config_za, mol->cell);

	return 0;
}
