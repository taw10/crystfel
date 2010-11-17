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
#include "beam-parameters.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Create reflections lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -t, --template=<filename>  Only include reflections mentioned in file.\n"
"      --poisson              Simulate Poisson samples.\n"
"      --noise                Add 10%% random noise.\n"
"  -y, --symmetry=<sym>       The symmetry of the input file (-i).\n"
"  -w, --twin=<sym>           Generate twinned data according to the given\n"
"                              point group.\n"
"  -e, --expand=<sym>         Expand reflections to this point group.\n"
"  -o, --output=<filename>    Output filename (default: stdout).\n"
"  -i, --intensities=<file>   Read intensities from file instead of\n"
"                              calculating them from scratch.  You might use\n"
"                              this if you need to apply noise or twinning.\n"
"  -p, --pdb=<file>           PDB file from which to get the structure.\n"
"      --no-phases            Do not try to use phases in the input file.\n"
"      --multiplicity         Multiply intensities by the number of\n"
"                              equivalent reflections.\n"
"  -b, --beam=<file>         Get beam parameters from file (used for sigmas).\n"
);
}


/* Apply Poisson noise to all reflections */
static void poisson_reflections(double *ref, ReflItemList *items)
{
	int i;
	const int n = num_items(items);

	for ( i=0; i<n; i++ ) {

		struct refl_item *it;
		double val;
		int c;

		it = get_item(items, i);

		val = lookup_intensity(ref, it->h, it->k, it->l);
		c = poisson_noise(val);
		set_intensity(ref, it->h, it->k, it->l, c);

		progress_bar(i, n-1, "Simulating noise");

	}
}


/* Apply 10% uniform noise to all reflections */
static void noise_reflections(double *ref, ReflItemList *items)
{
	int i;
	const int n = num_items(items);

	for ( i=0; i<n; i++ ) {

		struct refl_item *it;
		double val;
		double r;

		it = get_item(items, i);

		val = lookup_intensity(ref, it->h, it->k, it->l);

		r = (double)random()/RAND_MAX;
		val += 0.1 * val * r;

		set_intensity(ref, it->h, it->k, it->l, val);

		progress_bar(i, n-1, "Simulating noise");

	}
}


static ReflItemList *twin_reflections(double *ref, ReflItemList *items,
                                      const char *holo, const char *mero)
{
	int i;
	ReflItemList *new;

	new = new_items();

	if ( num_general_equivs(holo) < num_general_equivs(mero) ) {
		ERROR("%s is not a subgroup of %s!\n", mero, holo);
		return NULL;
	}

	for ( i=0; i<num_items(items); i++ ) {

		double mean;
		struct refl_item *it;
		signed int h, k, l;
		int n, j;
		int skip;

		it = get_item(items, i);

		/* There is a many-to-one correspondence between reflections
		 * in the merohedral and holohedral groups.  Do the calculation
		 * only once for each reflection in the holohedral group, which
		 * contains fewer reflections.
		 */
		get_asymm(it->h, it->k, it->l, &h, &k, &l, holo);
		if ( find_item(new, h, k, l) ) continue;

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


static ReflItemList *expand_reflections(double *ref, ReflItemList *items,
                                        const char *target, const char *initial)
{
	int i;
	ReflItemList *new;

	new = new_items();

	if ( num_general_equivs(target) > num_general_equivs(initial) ) {
		ERROR("%s is not a subgroup of %s!\n", initial, target);
		return NULL;
	}

	for ( i=0; i<num_items(items); i++ ) {

		struct refl_item *it;
		signed int h, k, l;
		signed int hd, kd, ld;
		int n, j;
		double intensity;

		it = get_item(items, i);
		h = it->h;  k = it->k;  l = it->l;

		/* Actually we don't really care what the equivalent is,
		 * we just want to be sure that there is nly be one version of
		 * this reflection. */
		find_unique_equiv(items, h, k, l, initial, &hd, &kd, &ld);

		/* Now find out how many reflections need to be filled in */
		n = num_equivs(h, k, l, initial);
		intensity = lookup_intensity(ref, h, k, l);

		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;

			/* Get the equivalent */
			get_equiv(h, k, l, &he, &ke, &le, initial, j);

			/* Put it into the asymmetric unit for the target */
			get_asymm(he, ke, le, &he, &ke, &le, target);

			/* Make sure the intensity is in the right place */
			set_intensity(ref, he, ke, le, intensity);

			/* Add the reflection, but only once */
			if ( !find_item(new, he, ke, le) ) {
				add_item(new, he, ke, le);
			}

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
	int config_noise = 0;
	int config_poisson = 0;
	int config_nophase = 0;
	int config_multi = 0;
	char *holo = NULL;
	char *mero = NULL;
	char *expand = NULL;
	char *output = NULL;
	char *input = NULL;
	char *filename = NULL;
	ReflItemList *input_items;
	ReflItemList *write_items;
	UnitCell *cell = NULL;
	double adu_per_photon;
	struct beam_params *beam = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"template",           1, NULL,               't'},
		{"poisson",            0, &config_poisson,     1},
		{"noise",              0, &config_noise,       1},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"twin",               1, NULL,               'w'},
		{"expand",             1, NULL,               'e'},
		{"intensities",        1, NULL,               'i'},
		{"pdb",                1, NULL,               'p'},
		{"no-phases",          0, &config_nophase,     1},
		{"multiplicity",       0, &config_multi,       1},
		{"beam",               1, NULL,               'b'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ht:o:i:p:w:y:e:",
	                        longopts, NULL)) != -1) {

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

		case 'e' :
			expand = strdup(optarg);
			break;

		case 'b' :
			beam = get_beam_parameters(optarg);
			if ( beam == NULL ) {
				ERROR("Failed to load beam parameters"
				      " from '%s'\n", optarg);
				return 1;
			}
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

	if ( (holo != NULL) && (expand != NULL) ) {
		ERROR("You cannot 'twin' and 'expand' at the same time.\n");
		ERROR("Decide which one you want to do first.\n");
		exit(1);
	}

	mol = load_molecule(filename);
	cell = load_cell_from_pdb(filename);
	if ( !config_nophase ) {
		phases = new_list_phase();
	} else {
		phases = NULL;
	}
	if ( input == NULL ) {
		input_items = new_items();
		ideal_ref = get_reflections(mol, eV_to_J(1790.0), 1/(0.05e-9),
		                            phases, input_items);
	} else {
		ideal_ref = new_list_intensity();
		input_items = read_reflections(input, ideal_ref, phases,
		                               NULL, NULL);
		free(input);
	}

	if ( config_poisson ) poisson_reflections(ideal_ref, input_items);
	if ( config_noise ) noise_reflections(ideal_ref, input_items);

	if ( holo != NULL ) {
		ReflItemList *new;
		STATUS("Twinning from %s into %s\n", mero, holo);
		new = twin_reflections(ideal_ref, input_items, holo, mero);
		delete_items(input_items);
		input_items = new;
	}

	if ( expand != NULL ) {
		ReflItemList *new;
		STATUS("Expanding from %s into %s\n", mero, expand);
		new = expand_reflections(ideal_ref, input_items, expand, mero);
		delete_items(input_items);
		input_items = new;
	}

	if ( config_multi ) {

		int i;

		for ( i=0; i<num_items(input_items); i++ ) {

			struct refl_item *it;
			double inty;

			it = get_item(input_items, i);
			inty = lookup_intensity(ideal_ref, it->h, it->k, it->l);
			inty *= num_equivs(it->h, it->k, it->l, mero);
			set_intensity(ideal_ref, it->h, it->k, it->l, inty);
			STATUS("%i %i %i %i\n", it->h, it->k, it->l,
			       num_equivs(it->h, it->k, it->l, mero));

		}
	}

	if ( template ) {
		/* Write out only reflections which are in the template
		 * (and which we have in the input) */
		ReflItemList *template_items;
		template_items = read_reflections(template,
		                                  NULL, NULL, NULL, NULL);
		write_items = intersection_items(input_items, template_items);
		delete_items(template_items);
	} else {
		/* Write out all reflections */
		write_items = new_items();
		/* (quick way of copying a list) */
		union_items(write_items, input_items);
	}

	if ( beam == NULL ) {
		adu_per_photon = 167.0;
		STATUS("No beam parameters file provided (use -b), "
		       "so I'm assuming 167.0 ADU per photon.\n");
	} else {
		adu_per_photon = beam->adu_per_photon;
	}

	write_reflections(output, write_items, ideal_ref, NULL, phases,
	                  NULL, cell);

	delete_items(input_items);
	delete_items(write_items);

	return 0;
}
