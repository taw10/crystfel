/*
 * process_hkl.c
 *
 * Assemble and process FEL Bragg intensities
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
#include "statistics.h"
#include "sfac.h"
#include "reflections.h"
#include "likelihood.h"
#include "symmetry.h"


/* Number of divisions for R vs |q| graphs */
#define RVQDV (20)


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Assemble and process FEL Bragg intensities.\n"
"\n"
"  -h, --help                Display this help message.\n"
"  -i, --input=<filename>    Specify input filename (\"-\" for stdin).\n"
"  -o, --output=<filename>   Specify output filename for merged intensities\n"
"                             (don't specify for no output).\n"
"  -p, --pdb=<filename>      PDB file to use (default: molecule.pdb).\n"
"\n"
"      --max-only            Take the integrated intensity to be equal to the\n"
"                             maximum intensity measured for that reflection.\n"
"                             The default is to use the mean value from all\n"
"                             measurements.\n"
"      --sum                 Sum (rather than average) the intensities for the\n"
"                             final output list.  This is useful for comparing\n"
"                             results to radially summed powder patterns, but\n"
"                             will break R-factor analysis.\n"
"      --stop-after=<n>      Stop after processing n patterns.  Zero means\n"
"                             keep going until the end of the input, and is\n"
"                             the default.\n"
"\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model.\n"
"  -y, --symmetry=<sym>      Merge according to point group <sym>.\n"
);
}


/* Note "holo" needn't actually be a holohedral point group, if you want to try
 * something strange like resolving from a low-symmetry group into an even
 * lower symmetry one.
 */
static ReflItemList *get_twin_possibilities(const char *holo, const char *mero)
{
	ReflItemList *test_items;
	ReflItemList *twins;
	int np;

	np = num_general_equivs(holo) / num_general_equivs(mero);

	test_items = new_items();

	/* Some arbitrarily chosen reflections which can't be special
	 * reflections in any point group, i.e. lots of odd numbers,
	 * prime numbers and so on.  There's probably an analytical
	 * way of working these out, but this will do. */
	add_item(test_items, 1, 2, 3);
	add_item(test_items, 3, 7, 13);
	add_item(test_items, 5, 2, 1);

	twins = get_twins(test_items, holo, mero);
	delete_items(test_items);

	if ( num_items(twins) != np ) {
		ERROR("Whoops! Couldn't find all the twinning possiblities.\n");
		abort();
	}

	return twins;
}


static int resolve_twin(ReflItemList *items, ReflItemList *twins,
                        const double *patt, const double *model,
                        const unsigned int *model_counts,
                        const char *holo, const char *mero)
{
	int n, i;
	double best_fom = 0.0;
	int best_op = 0;

	n = num_items(twins);

	for ( i=0; i<n; i++ ) {

		int j;
		int op;
		double *trial_ints = new_list_intensity();
		unsigned int *trial_counts = new_list_count();
		double fom;

		op = get_item(twins, i)->op;

		for ( j=0; j<num_items(items); j++ ) {

			signed int h, k, l;
			struct refl_item *r = get_item(items, j);

			get_general_equiv(r->h, r->k, r->l, &h, &k, &l,
			                  holo, op);
			get_asymm(h, k, l, &h, &k, &l, mero);

			set_intensity(trial_ints, h, k, l,
			              lookup_intensity(patt, r->h, r->k, r->l));
			set_count(trial_counts, h, k, l, 1);

		}

		fom = stat_pearson(trial_ints, trial_counts,
		                   model, model_counts);

		free(trial_ints);
		free(trial_counts);

		//printf(" %f", fom);
		if ( fom > best_fom ) {
			best_fom = fom;
			best_op = op;
		}

	}
	//printf("\n");

	return best_op;
}


static void merge_pattern(double *model, const double *new,
                          unsigned int *model_counts,
                          ReflItemList *items, int mo, int sum,
                          ReflItemList *twins,
                          const char *holo, const char *mero)
{
	int i;
	int twin;

	if ( twins != NULL ) {
		twin = resolve_twin(items, twins, new, model,
		                     model_counts, holo, mero);
	} else {
		twin = 0;
	}

	for ( i=0; i<num_items(items); i++ ) {

		double intensity;
		signed int hs, ks, ls;
		signed int h, k, l;
		struct refl_item *item;

		item = get_item(items, i);

		hs = item->h;
		ks = item->k;
		ls = item->l;

		get_general_equiv(hs, ks, ls, &h, &k, &l, holo, twin);
		get_asymm(h, k, l, &h, &k, &l, mero);

		intensity = lookup_intensity(new, h, k, l);

		if ( !mo ) {
			integrate_intensity(model, h, k, l, intensity);
			if ( sum ) {
				set_count(model_counts, h, k, l, 1);
			} else {
				integrate_count(model_counts, h, k, l, 1);
			}
		} else {
			if ( intensity > lookup_intensity(model, h, k, l) ) {
				set_intensity(model, h, k, l, intensity);
			}
			set_count(model_counts, h, k, l, 1);
		}

	}
}


static void merge_all(FILE *fh, double *model, unsigned int *model_counts,
                      int config_maxonly, int config_scale, int config_sum,
                      int config_stopafter,
                      ReflItemList *twins, const char *holo, const char *mero,
                      int n_total_patterns)
{
	char *rval;
	float f0;
	int n_nof0 = 0;
	int f0_valid = 0;
	int n_patterns = 0;
	double *new_pattern = new_list_intensity();
	ReflItemList *items = new_items();

	do {

		char line[1024];
		signed int h, k, l;
		float intensity;
		int r;

		rval = fgets(line, 1023, fh);
		if ( (strncmp(line, "Reflections from indexing", 25) == 0)
		    || (strncmp(line, "New pattern", 11) == 0) ) {

			/* Start of first pattern? */
			if ( n_patterns == 0 ) {
				n_patterns++;
				continue;
			}

			/* Assume a default I0 if we don't have one by now */
			if ( config_scale && !f0_valid ) {
				n_nof0++;
				f0 = 1.0;
			}

			/* Scale if requested */
			if ( config_scale ) {
				scale_intensities(model, new_pattern,
				                  model_counts, items, f0,
				                  f0_valid);
			}

			/* Start of second or later pattern */
			merge_pattern(model, new_pattern, model_counts,
			              items, config_maxonly, config_sum, twins,
			              holo, mero);

			if ( n_patterns == config_stopafter ) break;

			n_patterns++;
			clear_items(items);

			progress_bar(n_patterns, n_total_patterns, "Merging");

			f0_valid = 0;

		}

		if ( strncmp(line, "f0 = ", 5) == 0 ) {
			r = sscanf(line, "f0 = %f", &f0);
			if ( r != 1 ) {
				f0 = 1.0;
				f0_valid = 0;
				continue;
			}
			f0_valid = 1;
		}

		r = sscanf(line, "%i %i %i %f", &h, &k, &l, &intensity);
		if ( r != 4 ) continue;

		if ( (h==0) && (k==0) && (l==0) ) continue;

		if ( find_item(items, h, k, l) != 0 ) {
			ERROR("More than one measurement for %i %i %i in"
			      " pattern number %i\n", h, k, l, n_patterns);
		}
		set_intensity(new_pattern, h, k, l, intensity);
		add_item(items, h, k, l);

	} while ( rval != NULL );

	delete_items(items);
	free(new_pattern);

	STATUS("%i patterns had no f0 valid value.\n", n_nof0);
}


static int count_patterns(FILE *fh)
{
	char *rval;

	int n_total_patterns = 0;
	do {
		char line[1024];

		rval = fgets(line, 1023, fh);
		if ( (strncmp(line, "Reflections from indexing", 25) == 0)
		    || (strncmp(line, "New pattern", 11) == 0) ) {
		    n_total_patterns++;
		}
	} while ( rval != NULL );

	return n_total_patterns;
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	FILE *fh;
	double *model;
	unsigned int *model_counts;
	UnitCell *cell;
	int config_maxonly = 0;
	int config_stopafter = 0;
	int config_sum = 0;
	int config_scale = 0;
	unsigned int n_total_patterns;
	char *sym = NULL;
	char *pdb = NULL;
	ReflItemList *twins;
	int i;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"max-only",           0, &config_maxonly,     1},
		{"output-every",       1, NULL,               'e'},
		{"stop-after",         1, NULL,               's'},
		{"sum",                0, &config_sum,         1},
		{"scale",              0, &config_scale,       1},
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:ro:p:y:",
	                        longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			filename = strdup(optarg);
			break;

		case 'o' :
			output = strdup(optarg);
			break;

		case 's' :
			config_stopafter = atoi(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
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

	if ( filename == NULL ) {
		ERROR("Please specify filename using the -i option\n");
		return 1;
	}

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

	if ( sym == NULL ) {
		sym = strdup("1");
	}

	model = new_list_intensity();
	model_counts = new_list_count();

	cell = load_cell_from_pdb(pdb);
	free(pdb);

	/* Show useful symmetry information */
	const char *holo = get_holohedral(sym);
	int np = num_general_equivs(holo) / num_general_equivs(sym);
	if ( np > 1 ) {

		STATUS("Resolving point group %s into %s (%i possibilities)\n",
		       holo, sym, np);
		/* Get the list of twin/Bijvoet possibilities */
		twins = get_twin_possibilities(holo, sym);
		STATUS("Twin/inversion operation indices (from %s) are:", holo);
		for ( i=0; i<num_items(twins); i++ ) {
			STATUS(" %i", get_item(twins, i)->op);
		}
		STATUS("\n");

	} else {
		STATUS("No twin/inversion resolution necessary.\n");
		twins = NULL;
	}

	/* Open the data stream */
	if ( strcmp(filename, "-") == 0 ) {
		fh = stdin;
	} else {
		fh = fopen(filename, "r");
	}
	free(filename);
	if ( fh == NULL ) {
		ERROR("Failed to open input file\n");
		return 1;
	}

	/* Count the number of patterns in the file */
	n_total_patterns = count_patterns(fh);
	STATUS("There are %i patterns to process\n", n_total_patterns);
	rewind(fh);

	merge_all(fh, model, model_counts,
	          config_maxonly, config_scale, config_sum, config_stopafter,
                  twins, holo, sym, n_total_patterns);
	rewind(fh);

	fclose(fh);

	if ( output != NULL ) {
		write_reflections(output, model_counts, model, NULL,
		                  0, cell, 1);
	}

	free(sym);
	free(model);
	free(model_counts);
	free(output);
	free(cell);

	return 0;
}
