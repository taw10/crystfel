/*
 * get_hkl.c
 *
 * Small program to manipulate reflection lists
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
#include "reflist-utils.h"
#include "symmetry.h"
#include "beam-parameters.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Manipulate reflection lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -i, --input=<file>         Read reflections from <file>.\n"
"  -y, --symmetry=<sym>       The symmetry of the input reflection list.\n"
"\n"
"You can add noise to the reflections with either of:\n"
"      --poisson              Simulate Poisson samples.\n"
"      --noise                Add 10%% random noise.\n"
"\n"
"To calculate Poisson samples accurately, you must also give:\n"
"  -b, --beam=<file>          Get beam parameters from file.\n"
"\n"
"You can artificially 'twin' the reflections, or expand them out.  You can also"
" do both, in which case the 'twinning' will be done first:\n"
"  -w, --twin=<sym>           Generate twinned data according to the given\n"
"                              point group.\n"
"  -e, --expand=<sym>         Expand reflections to this point group.\n"
"\n"
"You can restrict which reflections are written out:\n"
"  -t, --template=<filename>  Only include reflections mentioned in file.\n"
"\n"
"You might sometimes need to do this:\n"
"      --multiplicity         Multiply intensities by the number of\n"
"                              equivalent reflections.\n"
"\n"
"Don't forget to specify the output filename:\n"
"  -o, --output=<filename>    Output filename (default: stdout).\n"
"  -p, --pdb=<filename>       Use unit cell parameters from this PDB file to\n"
"                              generate resolution values in the output file.\n"
);
}


/* Apply Poisson noise to all reflections */
static void poisson_reflections(RefList *list, double adu_per_photon)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double val, c;

		val = get_intensity(refl);

		c = adu_per_photon * poisson_noise(val/adu_per_photon);
		set_int(refl, c);

	}
}


/* Apply 10% uniform noise to all reflections */
static void noise_reflections(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double val, r;

		val = get_intensity(refl);

		r = (double)random()/RAND_MAX;
		val += 0.1 * val * r;

		set_int(refl, val);

	}
}


static RefList *template_reflections(RefList *list, RefList *template)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;

	out = reflist_new();

	for ( refl = first_refl(template, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		Reflection *new;
		Reflection *old;

		get_indices(refl, &h, &k, &l);

		old = find_refl(list, h, k, l);
		if ( old == NULL ) continue;

		new = add_refl(out, h, k, l);
		copy_data(new, old);

	}

	return out;
}


static RefList *twin_reflections(RefList *in,
                                 const char *holo, const char *mero)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;

	if ( num_general_equivs(holo) < num_general_equivs(mero) ) {
		ERROR("%s is not a subgroup of %s!\n", mero, holo);
		return NULL;
	}

	out = reflist_new();

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double total, sigma;
		signed int h, k, l;
		int n, j;
		int skip;

		get_indices(refl, &h, &k, &l);

		/* There is a many-to-one correspondence between reflections
		 * in the merohedral and holohedral groups.  Do the calculation
		 * only once for each reflection in the holohedral group, which
		 * contains fewer reflections.
		 */
		get_asymm(h, k, l, &h, &k, &l, holo);
		if ( find_refl(out, h, k, l) != NULL ) continue;

		total = 0.0;
		sigma = 0.0;
		skip = 0;
		n = num_equivs(h, k, l, holo);
		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;
			signed int hu, ku, lu;

			get_equiv(h, k, l, &he, &ke, &le, holo, j);

			/* Do we have this reflection?
			 * We might not have the particular (merohedral)
			 * equivalent which belongs to our definition of the
			 * asymmetric unit cell, so check them all.
			 */
			if ( !find_equiv_in_list(in, he, ke, le, mero,
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

			total += get_intensity(refl);
			sigma += pow(get_esd_intensity(refl), 2.0);

		}

		if ( !skip ) {

			Reflection *new = add_refl(out, h, k, l);
			set_int(new, total);
			set_esd_intensity(new, sqrt(sigma));

		}

	}

	return out;
}


static RefList *expand_reflections(RefList *in,
                                   const char *target, const char *initial)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;

	if ( num_general_equivs(target) > num_general_equivs(initial) ) {
		ERROR("%s is not a subgroup of %s!\n", initial, target);
		return NULL;
	}

	out = reflist_new();

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		int n, j;
		double intensity;

		get_indices(refl, &h, &k, &l);
		intensity = get_intensity(refl);

		n = num_equivs(h, k, l, initial);

		/* For each equivalent in the higher symmetry group */
		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;
			Reflection *new;

			/* Get the equivalent */
			get_equiv(h, k, l, &he, &ke, &le, initial, j);

			/* Put it into the asymmetric unit for the target */
			get_asymm(he, ke, le, &he, &ke, &le, target);

			/* Make sure the intensity is in the right place */
			new = add_refl(out, he, ke, le);
			copy_data(new, refl);

		}

	}

	return out;
}


int main(int argc, char *argv[])
{
	int c;
	int config_noise = 0;
	int config_poisson = 0;
	int config_multi = 0;
	char *holo = NULL;
	char *mero = NULL;
	char *expand = NULL;
	char *input_file = NULL;
	char *template = NULL;
	char *output = NULL;
	char *beamfile = NULL;
	struct beam_params *beam = NULL;
	RefList *input;
	UnitCell *cell = NULL;

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
		{"multiplicity",       0, &config_multi,       1},
		{"beam",               1, NULL,               'b'},
		{"pdb",                1, NULL,               'p'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ht:o:i:w:y:e:b:",
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
			input_file = strdup(optarg);
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
			beamfile = strdup(optarg);
			break;

		case 'p' :
			cell = load_cell_from_pdb(optarg);
			if ( cell == NULL ) {
				ERROR("Failed to get cell from '%s'\n", optarg);
				return 1;
			}
			break;

		case 0 :
			break;

		default :
			return 1;
		}

	}

	if ( (holo != NULL) && (expand != NULL) ) {
		ERROR("You cannot 'twin' and 'expand' at the same time.\n");
		ERROR("Decide which one you want to do first.\n");
		return 1;
	}

	if ( beamfile != NULL ) {
		beam = get_beam_parameters(beamfile);
		if ( beam == NULL ) {
			ERROR("Failed to load beam parameters from '%s'\n",
			      beamfile);
			return 1;
		}
	}

	if ( cell == NULL ) {
		ERROR("You need to give a PDB file with the unit cell.\n");
		return 1;
	}

	input = read_reflections(input_file);
	free(input_file);
	if ( check_list_symmetry(input, mero) ) {
		ERROR("The input reflection list does not appear to"
		      " have symmetry %s\n", mero);
		return 1;
	}

	if ( config_poisson ) {
		if ( beam != NULL ) {
			poisson_reflections(input, beam->adu_per_photon);
		} else {
			ERROR("You must give a beam parameters file in order"
			      " to calculate Poisson noise.\n");
			return 1;
		}
	}

	if ( config_noise ) noise_reflections(input);

	if ( holo != NULL ) {

		RefList *new;
		STATUS("Twinning from %s into %s\n", mero, holo);
		new = twin_reflections(input, holo, mero);

		/* Replace old with new */
		reflist_free(input);
		input = new;

		/* The symmetry of the list has changed */
		free(mero);
		mero = holo;

	}

	if ( expand != NULL ) {

		RefList *new;
		STATUS("Expanding from %s into %s\n", mero, expand);
		new = expand_reflections(input, expand, mero);

		/* Replace old with new */
		reflist_free(input);
		input = new;

	}

	if ( config_multi ) {

		Reflection *refl;
		RefListIterator *iter;

		for ( refl = first_refl(input, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			double inty;
			signed int h, k, l;

			get_indices(refl, &h, &k, &l);
			inty = get_intensity(refl);

			inty *= (double)num_equivs(h, k, l, mero);
			set_int(refl, inty);

		}
	}

	if ( template ) {

		RefList *t = read_reflections(template);
		RefList *new = template_reflections(input, t);
		reflist_free(input);
		input = new;

	}

	write_reflist(output, input, cell);

	reflist_free(input);

	return 0;
}
