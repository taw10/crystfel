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
#include "symmetry.h"


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
"  -y, --symmetry=<sym>       The symmetry of the input file (-i).\n"
"  -w, --twin=<sym>           Generate twinned data according to the given\n"
"                              point group.\n"
"  -o, --output=<filename>    Output filename (default: stdout).\n"
"  -i, --intensities=<file>   Read intensities from file instead of\n"
"                              calculating them from scratch.  You might use\n"
"                              this if you need to apply noise or twinning.\n"
"  -p, --pdb=<file>           PDB file from which to get the structure.\n"
);
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


static void scold_user_about_symmetry(signed int h, signed int k, signed int l,
                                      signed int he, signed int ke,
                                      signed int le)
{
	ERROR("Merohedrally equivalent reflection (%i %i %i) found for "
	      "%i %i %i.\n", he, ke, le, h, k, l);
	ERROR("This indicates that you lied to me about the symmetry of the "
	      "input reflections.  ");
	ERROR("I won't be able to give you a meaningful result in this "
	      "situation, so I'm going to give up right now.  ");
	ERROR("Please reconsider your previous processing of the data, and "
	      "perhaps try again with a lower symmetry for the '-y' option.\n");
}


static int find_unique_equiv(ReflItemList *items, signed int h, signed int k,
                             signed int l, const char *mero, signed int *hu,
                             signed int *ku, signed int *lu)
{
	int i;

	for ( i=0; i<num_equivs(h, k, l, mero); i++ ) {

		signed int he, ke, le;
		get_equiv(h, k, l, &he, &ke, &le, mero, i);
		if ( find_item(items, he, ke, le) ) {
			*hu = he;  *ku = ke;  *lu = le;
			return 1;
		}

	}

	return 0;
}


static ReflItemList *twin_reflections(double *ref, ReflItemList *items,
                                      const char *holo, const char *mero)
{
	int i;
	ReflItemList *new;

	new = new_items();

	for ( i=0; i<num_items(items); i++ ) {

		double mean;
		struct refl_item *it;
		signed int h, k, l;
		int n, j;
		int skip;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		/* None of the equivalent reflections should exist in the
		 * input dataset.  That would indicate that the user lied about
		 * the input symmetry.
		 *
		 * Start from j=1 to ignore the reflection itself.
		 */
		for ( j=1; j<num_equivs(h, k, l, mero); j++ ) {

			signed int he, ke, le;
			get_equiv(h, k, l, &he, &ke, &le, mero, j);
			if ( !find_item(items, he, ke, le) ) continue;

			scold_user_about_symmetry(h, k, l, he, ke, le);
			abort();

		}
		/* It doesn't matter if the reflection wasn't actually the one
		 * we define as being in the asymmetric unit cell, as long as
		 * things aren't confused by there being more than one of it.
		 */

		n = num_equivs(h, k, l, holo);

		mean = 0.0;
		skip = 0;
		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;
			signed int hu, ku, lu;

			get_equiv(h, k, l, &he, &ke, &le, holo, j);

			/* Do we have this reflection?
			 * We might not have the particular (merohedral)
			 * equivalent which belongs to our definition of the
			 * asymmetric unit cell, so check them all.
			 *
			 * We checked earlier that there's only one of these
			 * for each reflection.
			 */
			if ( !find_unique_equiv(items, he, ke, le, mero,
			                        &hu, &ku, &lu) ) {
				/* Don't have this reflection, so bail out */
				ERROR("Twinning %i %i %i requires the %i %i %i "
				      "reflection (or an equivalent in %s), "
				      "which I don't have. %i %i %i won't "
				      "appear in the output\n",
				      h, k, l, he, ke, le, mero, h, k, l);
				skip = 1;
				break;
			}

			mean += lookup_intensity(ref, hu, ku, lu);

		}

		if ( !skip ) {

			mean /= (double)n;

			set_intensity(ref, h, k, l, mean);
			add_item(new, h, k, l);

		}

	}

	return new;
}


int main(int argc, char *argv[])
{
	int c;
	double *ideal_ref;
	double *phases;
	struct molecule *mol;
	char *template = NULL;
	int config_noisify = 0;
	char *holo = NULL;
	char *mero = NULL;
	char *output = NULL;
	char *input = NULL;
	char *filename = NULL;
	ReflItemList *input_items;
	ReflItemList *write_items;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"template",           1, NULL,               't'},
		{"poisson",            0, &config_noisify,     1},
		{"output",             1, NULL,               'o'},
		{"twin",               1, NULL,               'w'},
		{"symmetry",           1, NULL,               'y'},
		{"intensities",        1, NULL,               'i'},
		{"pdb",                1, NULL,               'p'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ht:o:i:p:w:y:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 't' :
			template = strdup(optarg);
			break;

		case 'o' :
			output = strdup(optarg);
			break;

		case 'i' :
			input = strdup(optarg);
			break;

		case 'p' :
			filename = strdup(optarg);
			break;

		case 'y' :
			mero = strdup(optarg);
			break;

		case 'w' :
			holo = strdup(optarg);
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( filename == NULL ) {
		filename = strdup("molecule.pdb");
	}

	mol = load_molecule(filename);
	phases = new_list_phase();
	if ( input == NULL ) {
		input_items = new_items();
		ideal_ref = get_reflections(mol, eV_to_J(1790.0), 1/(0.05e-9),
		                            phases, input_items);
	} else {
		ideal_ref = new_list_intensity();
		phases = new_list_phase();
		input_items = read_reflections(input, ideal_ref, phases, NULL);
		free(input);
	}

	if ( config_noisify ) noisify_reflections(ideal_ref);

	if ( holo != NULL ) {
		ReflItemList *new;
		STATUS("Twinning from %s into %s\n", mero, holo);
		new = twin_reflections(ideal_ref, input_items, holo, mero);
		delete_items(input_items);
		input_items = new;
	}

	if ( template ) {
		/* Write out only reflections which are in the template
		 * (and which we have in the input) */
		ReflItemList *template_items;
		template_items = read_reflections(template, NULL, NULL, NULL);
		write_items = intersection_items(input_items, template_items);
		delete_items(template_items);
	} else {
		/* Write out all reflections */
		write_items = new_items();
		/* (quick way of copying a list) */
		union_items(write_items, input_items);
	}

	write_reflections(output, write_items, ideal_ref, phases, NULL,
	                  mol->cell);

	delete_items(input_items);
	delete_items(write_items);

	return 0;
}
