/*
 * get_hkl.c
 *
 * Small program to manipulate reflection lists
 *
 * Copyright © 2013-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2015 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
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

#include "version.h"
#include "utils.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "cell.h"
#include "cell-utils.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Manipulate reflection lists.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"\n"
"  -i, --input=<file>         Read reflections from <file>.\n"
"  -y, --symmetry=<sym>       The symmetry of the input reflection list.\n"
"  -p, --pdb=<file>           PDB file with cell parameters (needed when\n"
"                              using a resolution cutoff)\n"
"\n"
"You can add noise to the reflections with either of:\n"
"      --poisson              Simulate Poisson samples.\n"
"      --noise                Add 10%% random noise.\n"
"\n"
"To calculate Poisson samples accurately, you must also give:\n"
"      --adu-per-photon=<n>   Number of ADU per photon.\n"
"\n"
"You can artificially 'twin' the reflections, or expand them out.\n"
"  -w, --twin=<sym>           Generate twinned data according to the given\n"
"                              point group.\n"
"  -e, --expand=<sym>         Expand reflections to this point group.\n"
"      --no-need-all-parts    Output a twinned reflection even if not all\n"
"                              the necessary equivalents were present.\n"
"\n"
"You can reindex the reflections according to an operation, e.g. k,h,-l:\n"
"      --reindex=<op>         Reindex according to <op>.\n"
"\n"
"Use this option with care, and only if you understand why it might sometimes\n"
" be necessary:\n"
"  --trim-centrics            Remove reflections which are duplicated in the\n"
"                              point group specified with the '-y' option.\n"
"\n"
"You can restrict which reflections are written out:\n"
"  -t, --template=<filename>  Only include reflections mentioned in file.\n"
"      --cutoff-angstroms=<n> Only include reflections with d < n Angstroms.\n"
"\n"
"You might sometimes need to do this:\n"
"      --multiplicity         Multiply intensities by the number of\n"
"                              equivalent reflections.\n"
"\n"
"Don't forget to specify the output filename:\n"
"  -o, --output=<filename>    Output filename (default: stdout).\n"
);
}


/* Apply Poisson noise to all reflections */
static void poisson_reflections(RefList *list, double adu_per_photon)
{
	Reflection *refl;
	RefListIterator *iter;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double val, c;

		val = get_intensity(refl);

		c = adu_per_photon * poisson_noise(rng, val/adu_per_photon);
		set_intensity(refl, c);

	}

	gsl_rng_free(rng);
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

		set_intensity(refl, val);

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


static RefList *twin_reflections(RefList *in, int need_all_parts,
                                 const SymOpList *holo, const SymOpList *mero)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;
	SymOpMask *m;
	int n;

	out = reflist_new();

	/* No need to free and reallocate this for every reflection */
	m = new_symopmask(holo);

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double total, sigma;
		int multi;
		signed int h, k, l;
		int j;
		int skip;

		/* Figure out where to put the twinned version, and check it's
		 * not there already. */
		get_indices(refl, &h, &k, &l);
		get_asymm(holo, h, k, l, &h, &k, &l);
		if ( find_refl(out, h, k, l) != NULL ) continue;

		special_position(holo, m, h, k, l);
		n = num_equivs(holo, m);

		total = 0.0;
		sigma = 0.0;
		multi = 0;
		skip = 0;

		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;
			signed int hu, ku, lu;
			Reflection *part;
			int r;

			get_equiv(holo, m, j, h, k, l, &he, &ke, &le);
			get_asymm(mero, he, ke, le, &he, &ke, &le);

			/* Do we have this reflection?
			 * We might not have the particular (merohedral)
			 * equivalent which belongs to our definition of the
			 * asymmetric unit cell, so check them all.
			 */
			r = find_equiv_in_list(in, he, ke, le, mero,
			                       &hu, &ku, &lu);

			if ( need_all_parts && !r ) {

				ERROR("Twinning %i %i %i requires the %i %i %i "
				      "reflection (or an equivalent in %s), "
				      "which I don't have.\n",
				      h, k, l, he, ke, le, symmetry_name(mero));

				skip = 1;
				break;

			}

			if ( r ) {

				double i, sigi;
				int mult;

				part = find_refl(in, hu, ku, lu);

				i = get_intensity(part);
				sigi = get_esd_intensity(part);
				mult = get_redundancy(part);

				total += mult*i;
				sigma += pow(sigi*mult, 2.0);
				multi += mult;

				set_intensity(part, 0.0);
				set_esd_intensity(part, 0.0);
				set_redundancy(part, 0);
			}

		}

		if ( !skip ) {

			Reflection *new = add_refl(out, h, k, l);
			set_intensity(new, total/multi);
			set_esd_intensity(new, sqrt(sigma)/multi);
			set_redundancy(new, multi);

		}

	}

	return out;
}


static RefList *expand_reflections(RefList *in, const SymOpList *initial,
                                                const SymOpList *target)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;
	SymOpMask *m;
	int phase_warning = 0;

	if ( !is_subgroup(initial, target) ) {
		ERROR("%s is not a subgroup of %s!\n", symmetry_name(target),
		                                       symmetry_name(initial));
		return NULL;
	}

	out = reflist_new();
	m = new_symopmask(initial);

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		int n, j;

		get_indices(refl, &h, &k, &l);

		special_position(initial, m, h, k, l);
		n = num_equivs(initial, m);

		/* For each equivalent in the higher symmetry group */
		for ( j=0; j<n; j++ ) {

			signed int he, ke, le;
			Reflection *copy;
			int have_phase;
			double ph;

			/* Get the equivalent */
			get_equiv(initial, m, j, h, k, l, &he, &ke, &le);

			/* Put it into the asymmetric unit for the target */
			get_asymm(target, he, ke, le, &he, &ke, &le);

			if ( find_refl(out, he, ke, le) != NULL ) continue;

			/* Make sure the intensity is in the right place */
			copy = add_refl(out, he, ke, le);
			copy_data(copy, refl);

			ph = get_phase(refl, &have_phase);
			if ( have_phase ) {
				set_phase(copy, ph);
				if ( !phase_warning ) {
					ERROR("WARNING: get_hkl can't expand "
					      "phase values correctly when the "
					      "structure contains glides or "
					      "screw axes.\n");
					phase_warning = 1;
				}
			}

		}

	}

	free_symopmask(m);

	return out;
}


static RefList *trim_centrics(RefList *in, const SymOpList *sym)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;
	long long int nref = 0;
	long long int ntrim = 0;

	out = reflist_new();

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		signed int ha, ka, la;
		Reflection *new;

		get_indices(refl, &h, &k, &l);

		/* Put it into the asymmetric unit */
		get_asymm(sym, h, k, l, &ha, &ka, &la);
		nref++;

		new = find_refl(out, ha, ka, la);
		if ( new != NULL ) {
			ntrim++;
			continue;
		}

		/* Add new reflection under asymmetric (unique) indices */
		new = add_refl(out, ha, ka, la);
		copy_data(new, refl);
	}

	STATUS("Trimmed %lli out of %lli reflections.\n", ntrim, nref);

	return out;
}


int main(int argc, char *argv[])
{
	int c;
	int config_noise = 0;
	int config_poisson = 0;
	int config_multi = 0;
	int config_trimc = 0;
	int config_nap = 1;
	char *holo_str = NULL;
	char *mero_str = NULL;
	char *expand_str = NULL;
	SymOpList *holo = NULL;
	SymOpList *mero = NULL;
	SymOpList *expand = NULL;
	char *input_file = NULL;
	char *template = NULL;
	char *output = NULL;
	RefList *input;
	double adu_per_photon = 0.0;
	int have_adu_per_photon = 0;
	int have_cutoff_iso = 0;
	int have_cutoff_aniso = 0;
	char *cutoff_str = NULL;
	double cutiso = 0.0;
	float cutn1, cutn2, cutn3;
	char *cellfile = NULL;
	char *reindex_str = NULL;
	SymOpList *reindex = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                5 },
		{"template",           1, NULL,               't'},
		{"poisson",            0, &config_poisson,     1},
		{"noise",              0, &config_noise,       1},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"twin",               1, NULL,               'w'},
		{"expand",             1, NULL,               'e'},
		{"intensities",        1, NULL,               'i'},
		{"pdb",                1, NULL,               'p'},
		{"multiplicity",       0, &config_multi,       1},
		{"trim-centrics",      0, &config_trimc,       1},
		{"no-need-all-parts",  0, &config_nap,         0},
		{"adu-per-photon",     1, NULL,                2},
		{"cutoff-angstroms",   1, NULL,                3},
		{"reindex",            1, NULL,                4},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ht:o:i:w:y:e:p:",
	                        longopts, NULL)) != -1) {

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 5 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
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
			mero_str = strdup(optarg);
			break;

			case 'w' :
			holo_str = strdup(optarg);
			break;

			case 'e' :
			expand_str = strdup(optarg);
			break;

			case 'p' :
			cellfile = strdup(optarg);
			break;

			case 2 :
			adu_per_photon = strtof(optarg, NULL);
			have_adu_per_photon = 1;
			break;

			case 3 :
			cutoff_str = strdup(optarg);
			break;

			case 4 :
			reindex_str = strdup(optarg);
			break;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( cutoff_str != NULL ) {

		int r;

		r = sscanf(cutoff_str, "%f,%f,%f", &cutn1, &cutn2, &cutn3);
		if ( r == 3 ) {

			have_cutoff_aniso = 1;

			/* Convert Angstroms -> m */
			cutn1 /= 1e10;  cutn2 /= 1e10;  cutn3 /= 1e10;

		} else {

			char *rval;

			errno = 0;
			cutiso = strtod(cutoff_str, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --cutoff-angstroms.\n");
				return 1;
			}

			have_cutoff_iso = 1;

		}
		free(cutoff_str);

	}

	if ( (holo_str != NULL) && (expand_str != NULL) ) {
		ERROR("You cannot 'twin' and 'expand' at the same time.\n");
		ERROR("Decide which one you want to do first.\n");
		return 1;
	}

	if ( holo_str != NULL ) {
		holo = get_pointgroup(holo_str);
		free(holo_str);
	} else {
		holo = NULL;
	}
	if ( mero_str != NULL ) {
		mero = get_pointgroup(mero_str);
		free(mero_str);
	} else {
		mero = NULL;
	}
	if ( expand_str != NULL ) {
		expand = get_pointgroup(expand_str);
		free(expand_str);
	} else {
		expand = NULL;
	}

	if ( reindex_str != NULL ) {
		reindex = parse_symmetry_operations(reindex_str);
		if ( reindex == NULL ) return 1;
		set_symmetry_name(reindex, "Reindex");
	}

	if ( (expand != NULL) || (holo != NULL) || config_trimc
	  || config_multi ) {
		if ( mero == NULL ) {
			ERROR("You must specify the point group with -y.\n");
		}
	}

	input = read_reflections(input_file);
	if ( input == NULL ) {
		ERROR("Problem reading input file %s\n", input_file);
		return 1;
	}
	free(input_file);

	STATUS("%i reflections in input.\n", num_reflections(input));

	if ( (mero != NULL) && !config_trimc
	  && check_list_symmetry(input, mero) )
	{
		ERROR("The input reflection list does not appear to"
		      " have symmetry %s\n", symmetry_name(mero));
		ERROR("If your unit cell is monoclinic, you may need to specify"
		      " the unique axis for your point group.  The default is"
		      " unique axis c.\n");
		ERROR("See 'man crystfel' for more details.\n");
		return 1;
	}

	if ( config_poisson ) {
		if ( have_adu_per_photon ) {
			poisson_reflections(input, adu_per_photon);
		} else {
			ERROR("You must give the number of ADU per photon to "
			      "use --poisson.\n");
			return 1;
		}
	}

	if ( config_noise ) noise_reflections(input);

	if ( holo != NULL ) {

		RefList *new;
		STATUS("Twinning from %s into %s\n", symmetry_name(mero),
		                                     symmetry_name(holo));
		new = twin_reflections(input, config_nap, holo, mero);

		/* Replace old with new */
		reflist_free(input);
		input = new;

		/* The symmetry of the list has changed */
		free(mero);
		mero = holo;

	}

	if ( expand != NULL ) {

		RefList *new;
		STATUS("Expanding from %s into %s\n", symmetry_name(mero),
		                                      symmetry_name(expand));
		new = expand_reflections(input, mero, expand);

		/* Replace old with new */
		reflist_free(input);
		input = new;

	}

	if ( config_trimc ) {

		RefList *new;

		/* Can't do this if point group is invalid */
		if ( mero == NULL ) {
			ERROR("Need point group to trim centrics.\n");
			return 1;
		}

		STATUS("Trimming duplicate reflections in %s\n",
		       symmetry_name(mero));
		new = trim_centrics(input, mero);
		reflist_free(input);
		input = new;
		STATUS("%i output reflections\n", num_reflections(input));

	}

	if ( config_multi ) {

		Reflection *refl;
		RefListIterator *iter;
		SymOpMask *m;

		m = new_symopmask(mero);

		for ( refl = first_refl(input, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			double inty;
			signed int h, k, l;

			get_indices(refl, &h, &k, &l);
			inty = get_intensity(refl);

			special_position(mero, m, h, k, l);
			inty *= (double)num_equivs(mero, m);
			set_intensity(refl, inty);

		}

		free_symopmask(m);
	}

	if ( template ) {

		RefList *t = read_reflections(template);
		RefList *new = template_reflections(input, t);
		reflist_free(input);
		input = new;

	}

	if ( have_cutoff_iso ) {

		RefList *n;
		Reflection *refl;
		RefListIterator *iter;
		UnitCell *cell;

		if ( cellfile == NULL ) {
			ERROR("You must provide a unit cell when using "
			      "--cutoff-angstroms.\n");
			return 1;
		}

		cell = load_cell_from_file(cellfile);
		if ( cell == NULL ) {
			ERROR("Failed to load cell from '%s'\n", cellfile);
			return 1;
		}
		free(cellfile);

		n = reflist_new();

		for ( refl = first_refl(input, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			double res;
			get_indices(refl, &h, &k, &l);
			res = 2.0 * resolution(cell, h, k, l);
			if ( res < 1e10 / cutiso ) {
				Reflection *a;
				a = add_refl(n, h, k, l);
				copy_data(a, refl);
			}
		}

		cell_free(cell);
		reflist_free(input);
		input = n;

	}

	if ( have_cutoff_aniso ) {

		RefList *n;
		Reflection *refl;
		RefListIterator *iter;
		UnitCell *cell;
		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		double as, bs, cs;

		if ( cellfile == NULL ) {
			ERROR("You must provide a unit cell when using "
			      "--cutoff-angstroms.\n");
			return 1;
		}

		cell = load_cell_from_file(cellfile);
		if ( cell == NULL ) {
			ERROR("Failed to load cell from '%s'\n", cellfile);
			return 1;
		}
		free(cellfile);

		cell_get_reciprocal(cell, &asx, &asy, &asz,
				          &bsx, &bsy, &bsz,
		                          &csx, &csy, &csz);
		as = modulus(asx, asy, asz);
		bs = modulus(bsx, bsy, bsz);
		cs = modulus(csx, csy, csz);

		n = reflist_new();

		for ( refl = first_refl(input, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			double sum;

			get_indices(refl, &h, &k, &l);

			sum  = pow(h*as*cutn1, 2.0);
			sum += pow(k*bs*cutn2, 2.0);
			sum += pow(l*cs*cutn3, 2.0);

			if ( sum < 1.0 ) {
				Reflection *a;
				a = add_refl(n, h, k, l);
				copy_data(a, refl);
			}
		}

		cell_free(cell);
		reflist_free(input);
		input = n;

	}

	if ( reindex != NULL ) {

		RefList *n;
		Reflection *refl;
		RefListIterator *iter;

		n = reflist_new();

		for ( refl = first_refl(input, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			Reflection *rn;

			get_indices(refl, &h, &k, &l);
			get_equiv(reindex, NULL, 0, h, k, l, &h, &k, &l);
			rn = add_refl(n, h, k, l);
			copy_data(rn, refl);

		}

		reflist_free(input);
		input = n;

	}

	write_reflist(output, input);

	reflist_free(input);

	return 0;
}
