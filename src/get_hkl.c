/*
 * get_hkl.c
 *
 * Small program to write out a list of h,k,l,I values given a structure
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
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
"  -h, --help              Display this help message.\n");
}


int main(int argc, char *argv[])
{
	int c;
	double *ref;
	struct molecule *mol;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:r", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
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
	get_reflections_cached(mol, eV_to_J(2.0e3));
	ref = ideal_intensities(mol->reflections);
	write_reflections("results/ideal-reflections.hkl", NULL, ref, 1,
	                  mol->cell);

	return 0;
}
