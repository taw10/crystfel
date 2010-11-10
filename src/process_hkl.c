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
#include "symmetry.h"
#include "stream.h"
#include "beam-parameters.h"


/* Number of divisions for intensity histograms */
#define NBINS (50)


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
"  -b, --beam=<file>         Get beam parameters from file (used for sigmas).\n"
"\n"
"      --max-only            Take the integrated intensity to be equal to the\n"
"                             maximum intensity measured for that reflection.\n"
"                             The default is to use the mean value from all\n"
"                             measurements.\n"
"      --sum                 Sum (rather than average) the intensities for the\n"
"                             final output list.  This is useful for comparing\n"
"                             results to radially summed powder patterns, but\n"
"                             will break R-factor analysis.\n"
"      --start-after=<n>     Skip n patterns at the start of the stream.\n"
"      --stop-after=<n>      Stop after processing n patterns.  Zero means\n"
"                             keep going until the end of the input, and is\n"
"                             the default.\n"
"  -g, --histogram=<h,k,l>   Calculate the histogram of measurements for this\n"
"                             reflection.\n"
"      --rmerge              Calculate and report Rmerge and Rpim\n"
"\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model.\n"
"  -y, --symmetry=<sym>      Merge according to point group <sym>.\n"
"      --reference=<file>    Compare against intensities from <file> when\n"
"                             scaling or resolving ambiguities.\n"
"                             The symmetry of the reference list must be the\n"
"                             same as that given with '-y'.\n"
"      --outstream=<file>    Write an annotated version of the input stream\n"
"                             to <file>.\n"
);
}


static void plot_histogram(double *vals, int n)
{
	int i;
	double max = -INFINITY;
	double min = +INFINITY;
	double step;
	int histo[NBINS];
	FILE *fh;

	fh = fopen("histogram.dat", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open 'histogram.dat'\n");
		return;
	}

	for ( i=0; i<n; i++ ) {
		if ( vals[i] > max ) max = vals[i];
		if ( vals[i] < min ) min = vals[i];
	}
	STATUS("%f %f\n", min, max);

	for ( i=0; i<NBINS; i++ ) {
		histo[i] = 0;
	}

	step = (max-min)/NBINS;

	for ( i=0; i<n; i++ ) {
		int bin;
		bin = (vals[i]-min)/step;
		histo[bin]++;
	}

	for ( i=0; i<NBINS; i++ ) {
		fprintf(fh, "%f %i\n", min+step*i, histo[i]);
	}

	fclose(fh);
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

	/* Idiot check.  Wouldn't be necessary if I could prove that the above
	 * set of arbitrarily chosen reflections were always general. */
	if ( num_items(twins) != np ) {
		ERROR("Whoops! Couldn't find all the twinning possiblities.\n");
		abort();
	}

	return twins;
}


static int resolve_twin(const double *model, ReflItemList *observed,
                        const double *patt, ReflItemList *items,
                        ReflItemList *twins, const char *holo, const char *mero)
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
		ReflItemList *intersection;

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

		intersection = intersection_items(observed, items);
		fom = stat_pearson_i(trial_ints, model, intersection);
		delete_items(intersection);

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


static void merge_pattern(double *model, ReflItemList *observed,
                          const double *new,  ReflItemList *items,
                          unsigned int *model_counts,  int mo,
                          ReflItemList *twins,
                          const char *holo, const char *mero,
                          ReflItemList *reference, const double *reference_i,
                          double *hist_vals,
                          signed int hist_h, signed int hist_k,
                          signed int hist_l, int *hist_n, double *devs,
                          double *tots, double *means, FILE *outfh)
{
	int i;
	int twin;
	ReflItemList *sym_items = new_items();

	if ( twins != NULL ) {
		if ( reference != NULL ) {
			twin = resolve_twin(reference_i, reference, new, items,
			                    twins, holo, mero);
		} else {
			twin = resolve_twin(model, observed, new, items,
			                    twins, holo, mero);
		}
	} else {
		twin = 0;
	}

	if ( outfh != NULL ) {
		fprintf(outfh, "twin=%i\n", twin);
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

		/* Transform into correct side of the twin law.
		 * "twin" is always zero if no de-twinning is performed. */
		get_general_equiv(hs, ks, ls, &h, &k, &l, holo, twin);

		/* Put into the asymmetric unit for the target group */
		get_asymm(h, k, l, &h, &k, &l, mero);

		/* Read the intensity from the original location
		 * (i.e. before screwing around with symmetry) */
		intensity = lookup_intensity(new, hs, ks, ls);

		/* User asked for max only? */
		if ( !mo ) {
			integrate_intensity(model, h, k, l, intensity);
		} else {
			if ( intensity > lookup_intensity(model, h, k, l) ) {
				set_intensity(model, h, k, l, intensity);
			}
		}

		if ( tots != NULL ) {
			double m;
			m = lookup_intensity(means, h, k, l);
			integrate_intensity(tots, h, k, l, intensity);
			integrate_intensity(devs, h, k, l, fabs(intensity-m));
		}

		if ( !find_item(sym_items, h, k, l) ) {
			add_item(sym_items, h, k, l);
		}

		/* Increase count count */
		integrate_count(model_counts, h, k, l, 1);

		if ( hist_vals != NULL ) {
			int p = *hist_n;
			if ( (h==hist_h) && (k==hist_k) && (l==hist_l) ) {
				hist_vals[p] = intensity;
				*hist_n = p+1;
			}
		}

	}

	/* Dump the reflections in this pattern into the overall list */
	union_items(observed, sym_items);

	delete_items(sym_items);
}


enum {
	SCALE_NONE,
	SCALE_CONSTINT,
	SCALE_INTPERBRAGG,
	SCALE_TWOPASS,
};


static void scale_intensities(const double *model, ReflItemList *model_items,
                              double *new_pattern, ReflItemList *new_items,
                              double f0, int f0_valid, const char *sym)
{
	double s;
	double top = 0.0;
	double bot = 0.0;
	unsigned int i;
	const int scaling = SCALE_INTPERBRAGG;

	for ( i=0; i<num_items(new_items); i++ ) {

		double i1, i2;
		struct refl_item *it;
		signed int hu, ku, lu;

		/* Get the next item in the list of new reflections */
		it = get_item(new_items, i);

		switch ( scaling ) {
		case SCALE_TWOPASS :

			/* Find the (only) partner in the model */
			find_unique_equiv(model_items, it->h, it->k, it->l, sym,
				          &hu, &ku, &lu);

			i1 = lookup_intensity(model, hu, ku, lu);
			i2 = lookup_intensity(new_pattern, it->h, it->k, it->l);

			/* Calculate LSQ estimate of scaling factor */
			top += i1 * i2;
			bot += i2 * i2;

			break;

		case SCALE_CONSTINT :

			/* Sum up the intensity in the pattern */
			i2 = lookup_intensity(new_pattern, it->h, it->k, it->l);
			top += i2;

			break;

		case SCALE_INTPERBRAGG :

			/* Sum up the intensity in the pattern */
			i2 = lookup_intensity(new_pattern, it->h, it->k, it->l);
			top += i2;
			bot += 1.0;

			break;

		}

	}

	switch ( scaling ) {
	case SCALE_TWOPASS :
		s = top / bot;
		break;
	case SCALE_CONSTINT :
		s = 1000.0 / top;
		break;
	case SCALE_INTPERBRAGG :
		s = 1000.0 / (top/bot);
		break;
	}
	//if ( f0_valid ) printf("%f %f\n", s, f0);

	/* Multiply the new pattern up by "s" */
	for ( i=0; i<LIST_SIZE; i++ ) {
		new_pattern[i] *= s;
	}
}


static void merge_all(FILE *fh, double **pmodel, ReflItemList **pobserved,
                      unsigned int **pcounts,
                      int config_maxonly, int config_scale, int config_sum,
                      int config_startafter, int config_stopafter,
                      ReflItemList *twins, const char *holo, const char *mero,
                      int n_total_patterns,
                      ReflItemList *reference, double *reference_i,
                      double *hist_vals,
                      signed int hist_h, signed int hist_k, signed int hist_l,
                      int *hist_i, double *devs, double *tots, double *means,
                      FILE *outfh)
{
	char *rval;
	float f0;
	int n_nof0 = 0;
	int f0_valid = 0;
	int n_patterns = 0;
	double *new_pattern = new_list_intensity();
	ReflItemList *items = new_items();
	ReflItemList *observed = new_items();
	double *model = new_list_intensity();
	unsigned int *counts = new_list_count();
	int i;

	if ( config_startafter != 0 ) {

		do {

			char line[1024];

			rval = fgets(line, 1023, fh);
			if ( (strncmp(line, "Reflections from indexing", 25)
			                                                   == 0)
			    || (strncmp(line, "New pattern", 11) == 0) ) {
				n_patterns++;
			}

			if ( n_patterns == config_startafter ) break;

		} while ( rval != NULL );

	}

	do {

		char line[1024];
		signed int h, k, l;
		float intensity;
		int r;

		rval = fgets(line, 1023, fh);
		if ( ((strncmp(line, "Reflections from indexing", 25) == 0)
		  || (rval == NULL)) && (num_items(items)>0) ) {

			/* Assume a default I0 if we don't have one by now */
			if ( config_scale && !f0_valid ) {
				n_nof0++;
				f0 = 1.0;
			}

			/* Scale if requested */
			if ( config_scale ) {
				if ( reference == NULL ) {
					scale_intensities(model, observed,
					                  new_pattern, items,
					                  f0, f0_valid, mero);
				} else {
					scale_intensities(reference_i,
					                  reference,
					                  new_pattern, items,
					                  f0, f0_valid, mero);
				}
			}

			/* Start of second or later pattern */
			merge_pattern(model, observed, new_pattern, items,
			              counts, config_maxonly,
			              twins, holo, mero,
			              reference, reference_i,
			              hist_vals, hist_h, hist_k, hist_l,
			              hist_i, devs, tots, means, outfh);

			n_patterns++;
			if ( n_patterns == config_stopafter ) break;
			progress_bar(n_patterns, n_total_patterns, "Merging");

			/* Reset for the next pattern */
			clear_items(items);

			f0_valid = 0;

		}

		if ( outfh != NULL ) {
			fprintf(outfh, "%s", line);
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

		/* Not interested in the central beam */
		if ( (h==0) && (k==0) && (l==0) ) continue;

		/* The same raw indices (before mapping into the asymmetric
		 * unit should not turn up twice in one pattern. */
		if ( find_item(items, h, k, l) != 0 ) {
			ERROR("More than one measurement for %i %i %i in"
			      " pattern number %i\n", h, k, l, n_patterns);
		}
		set_intensity(new_pattern, h, k, l, intensity);

		/* NB: This list contains raw indices, before working out
		 * where they belong in the asymmetric unit. */
		add_item(items, h, k, l);

	} while ( rval != NULL );

	delete_items(items);
	free(new_pattern);

	/* Calculate mean intensity if necessary */
	if ( !config_sum && !config_maxonly ) {
		for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
			if ( counts[i] > 0 ) {
				model[i] /= (double)counts[i];
			}
		}
	}

	*pmodel = model;
	*pcounts = counts;
	*pobserved = observed;

	STATUS("%i patterns had no f0 valid value.\n", n_nof0);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	FILE *fh;
	double *model;
	unsigned int *counts;
	UnitCell *cell = NULL;
	int config_maxonly = 0;
	int config_startafter = 0;
	int config_stopafter = 0;
	int config_sum = 0;
	int config_scale = 0;
	int config_rmerge = 0;
	unsigned int n_total_patterns;
	char *sym = NULL;
	char *pdb = NULL;
	ReflItemList *twins;
	ReflItemList *observed;
	int i;
	const char *holo = NULL;
	char *histo = NULL;
	signed int hist_h, hist_k, hist_l;
	double *hist_vals = NULL;
	int hist_i;
	char *outstream = NULL;
	char *reference = NULL;
	ReflItemList *reference_items;
	double *reference_i;
	FILE *outfh = NULL;
	struct beam_params *beam = NULL;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"max-only",           0, &config_maxonly,     1},
		{"output-every",       1, NULL,               'e'},
		{"stop-after",         1, NULL,               's'},
		{"start-after",        1, NULL,               'f'},
		{"sum",                0, &config_sum,         1},
		{"scale",              0, &config_scale,       1},
		{"symmetry",           1, NULL,               'y'},
		{"pdb",                1, NULL,               'p'},
		{"histogram",          1, NULL,               'g'},
		{"rmerge",             0, &config_rmerge,      1},
		{"outstream",          1, NULL,               'a'},
		{"reference",          1, NULL,               'r'},
		{"beam",               1, NULL,               'b'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:ro:p:y:g:f:a:r:b:",
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

		case 'f' :
			config_startafter = atoi(optarg);
			break;

		case 'p' :
			pdb = strdup(optarg);
			break;

		case 'y' :
			sym = strdup(optarg);
			break;

		case 'g' :
			histo = strdup(optarg);
			break;

		case 'r' :
			reference = strdup(optarg);
			break;

		case 'a' :
			outstream = strdup(optarg);
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

	if ( config_sum && config_rmerge ) {
		ERROR("Options --sum and --rmerge do not make sense"
		      " together.\n");
		return 1;
	}

	if ( filename == NULL ) {
		ERROR("Please specify filename using the -i option\n");
		return 1;
	}

	if ( pdb == NULL ) {
		pdb = strdup("molecule.pdb");
	}

	cell = load_cell_from_pdb(pdb);
	free(pdb);

	/* Show useful symmetry information */
	if ( sym != NULL ) {
		holo = get_holohedral(sym);
		int np = num_general_equivs(holo) / num_general_equivs(sym);
		if ( np > 1 ) {

			STATUS("Resolving point group %s into %s "
			       "(%i possibilities)\n",
			       holo, sym, np);
			/* Get the list of twin/Bijvoet possibilities */
			twins = get_twin_possibilities(holo, sym);
			STATUS("Twin/inversion operation indices from %s are:",
			       holo);
			for ( i=0; i<num_items(twins); i++ ) {
				STATUS(" %i", get_item(twins, i)->op);
			}
			STATUS("\n");

		} else {
			STATUS("No twin/inversion resolution necessary.\n");
			twins = NULL;
		}
	} else {
		STATUS("Not performing any twin/inversion resolution.\n");
		twins = NULL;
		sym = strdup("1");
		holo = strdup("1");
	}

	if ( histo != NULL ) {
		int r;
		r = sscanf(histo, "%i,%i,%i", &hist_h, &hist_k, &hist_l);
		if ( r != 3 ) {
			ERROR("Invalid indices for '--histogram'\n");
			return 1;
		}
		hist_vals = malloc(10*1024*sizeof(double));
		free(histo);
		STATUS("Histogramming %i %i %i -> ", hist_h, hist_k, hist_l);
		/* Put into the asymmetric cell for the target group */
		get_asymm(hist_h, hist_k, hist_l,
		          &hist_h, &hist_k, &hist_l, sym);
		STATUS("%i %i %i\n", hist_h, hist_k, hist_l);
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

	/* Read the reference reflections */
	if ( reference != NULL ) {
		reference_i = new_list_intensity();
		reference_items = read_reflections(reference, reference_i,
		                                   NULL, NULL, NULL);
		if ( reference_items == NULL ) {
			ERROR("Couldn't read '%s'\n", reference);
			return 1;
		}
	} else {
		reference_items = NULL;
		reference_i = NULL;
	}

	if ( outstream != NULL ) {
		outfh = fopen(outstream, "w");
		if ( outfh == NULL ) {
			ERROR("Couldn't open '%s'\n", outstream);
			return 1;
		}
	}

	/* Count the number of patterns in the file */
	n_total_patterns = count_patterns(fh);
	STATUS("There are %i patterns to process\n", n_total_patterns);
	rewind(fh);

	hist_i = 0;
	merge_all(fh, &model, &observed, &counts,
	          config_maxonly, config_scale, config_sum,
	          config_startafter, config_stopafter,
                  twins, holo, sym, n_total_patterns,
                  reference_items, reference_i,
                  hist_vals, hist_h, hist_k, hist_l, &hist_i, NULL, NULL, NULL,
                  outfh);
	rewind(fh);

	if ( hist_vals != NULL ) {
		STATUS("%i %i %i was seen %i times.\n", hist_h, hist_k, hist_l,
		                                        hist_i);
		plot_histogram(hist_vals, hist_i);
	}

	if ( output != NULL ) {

		double adu_per_photon;

		if ( beam == NULL ) {
			adu_per_photon = 167.0;
			STATUS("No beam parameters file provided (use -b), "
			       "so I'm assuming 167.0 ADU per photon.\n");
		} else {
			adu_per_photon = beam->adu_per_photon;
		}

		write_reflections(output, observed, model, NULL, counts, cell,
		                  adu_per_photon);
	}

	if ( config_rmerge ) {

		double *devs = new_list_intensity();
		double *tots = new_list_intensity();
		double total_dev = 0.0;
		double total_tot = 0.0;
		double total_pdev = 0.0;
		double total_rdev = 0.0;
		int i;

		STATUS("Extra pass to calculate figures of merit...\n");

		rewind(fh);
		merge_all(fh, &model, &observed, &counts,
		          config_maxonly, config_scale, 0,
		          config_startafter, config_stopafter, twins, holo, sym,
		          n_total_patterns, reference_items, reference_i,
		          NULL, 0, 0, 0, NULL, devs, tots, model, NULL);

		for ( i=0; i<num_items(observed); i++ ) {

			struct refl_item *it;
			signed int h, k, l;
			double dev;
			unsigned int count;

			it = get_item(observed, i);
			h = it->h;  k = it->k,  l = it->l;

			dev = lookup_intensity(devs, h, k, l);
			count = lookup_count(counts, h, k, l);

			if ( count < 2 ) continue;

			total_dev += dev;
			total_pdev += sqrt(1.0/(count-1.0))*dev;
			total_rdev += sqrt(count/(count-1.0))*dev;
			total_tot += lookup_intensity(tots, h, k, l);

		}

		STATUS("Rmerge = %f%%\n", 100.0*total_dev/total_tot);
		STATUS("  Rpim = %f%%\n", 100.0*total_pdev/total_tot);
		STATUS(" Rmeas = %f%%\n", 100.0*total_rdev/total_tot);

	}

	fclose(fh);
	free(sym);
	free(model);
	free(counts);
	free(output);
	if ( cell != NULL ) cell_free(cell);

	return 0;
}
