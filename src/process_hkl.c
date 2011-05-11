/*
 * process_hkl.c
 *
 * Assemble and process FEL Bragg intensities
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 * (c) 2011      Andrew Martin <andrew.martin@desy.de>
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
#include "reflist-utils.h"
#include "symmetry.h"
#include "stream.h"
#include "beam-parameters.h"
#include "reflist.h"
#include "image.h"


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
"                             Default: processed.hkl).\n"
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
	min--;  max++;

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


static void merge_pattern(RefList *model, RefList *new, int max_only,
                          const char *sym,
                          double *hist_vals, signed int hist_h,
                          signed int hist_k, signed int hist_l, int *hist_n,
                          int pass)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(new, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double intensity;
		signed int h, k, l;
		Reflection *model_version;
		double model_int;

		get_indices(refl, &h, &k, &l);

		/* Put into the asymmetric unit for the target group */
		get_asymm(h, k, l, &h, &k, &l, sym);

		model_version = find_refl(model, h, k, l);
		if ( model_version == NULL ) {
			model_version = add_refl(model, h, k, l);
		}

		/* Read the intensity from the original location
		 * (i.e. before screwing around with symmetry) */
		intensity = get_intensity(refl);

		/* Get the current model intensity */
		model_int = get_intensity(model_version);

		if ( pass == 1 ) {

			/* User asked for max only? */
			if ( !max_only ) {
				set_int(model_version, model_int + intensity);
			} else {
				if ( intensity>get_intensity(model_version) ) {
					set_int(model_version, intensity);
				}
			}


			/* Increase redundancy */
			int cur_redundancy = get_redundancy(model_version);
			set_redundancy(model_version, cur_redundancy+1);

		} else if ( pass == 2 ) {

			double dev = get_sum_squared_dev(model_version);

			/* Other ways of estimating the sigma are possible,
			 * choose from:
			 *    dev += pow(get_esd_intensity(refl), 2.0);
			 *    dev += pow(intensity, 2.0);
			 * But alter the other part of the calculation below
			 * as well. */
			dev += pow(intensity - model_int, 2.0);

			set_sum_squared_dev(model_version, dev);

			if ( hist_vals != NULL ) {
				int p = *hist_n;
				if ( (h==hist_h) && (k==hist_k)
				  && (l==hist_l) )
				{
					hist_vals[p] = intensity;
					*hist_n = p+1;
				}

			}

		}

	}


}


enum {
	SCALE_NONE,
	SCALE_CONSTINT,
	SCALE_INTPERBRAGG,
	SCALE_TWOPASS,
};


static void scale_intensities(RefList *model, RefList *new, const char *sym)
{
	double s;
	double top = 0.0;
	double bot = 0.0;
	const int scaling = SCALE_INTPERBRAGG;
	Reflection *refl;
	RefListIterator *iter;
	Reflection *model_version;

	for ( refl = first_refl(new, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double i1, i2;
		signed int hu, ku, lu;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		switch ( scaling ) {
		case SCALE_TWOPASS :

			model_version = find_refl(model, h, k, l);
			if ( model_version == NULL ) continue;

			get_asymm(h, k, l, &hu, &ku, &lu, sym);

			i1 = get_intensity(model_version);
			i2 = get_intensity(refl);

			/* Calculate LSQ estimate of scaling factor */
			top += i1 * i2;
			bot += i2 * i2;

			break;

		case SCALE_CONSTINT :

			/* Sum up the intensity in the pattern */
			i2 = get_intensity(refl);
			top += i2;

			break;

		case SCALE_INTPERBRAGG :

			/* Sum up the intensity in the pattern */
			i2 = get_intensity(refl);
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

	/* Multiply the new pattern up by "s" */
	for ( refl = first_refl(new, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double intensity = get_intensity(refl);
		set_int(refl, intensity*s);

	}
}


static void merge_all(FILE *fh, RefList *model,
                      int config_maxonly, int config_scale, int config_sum,
                      int config_startafter, int config_stopafter,
                      const char *sym,
                      int n_total_patterns,
                      double *hist_vals, signed int hist_h,
                      signed int hist_k, signed int hist_l,
                      int *hist_i, int pass)
{
	int rval;
	int n_patterns = 0;
	int n_used = 0;
	Reflection *refl;
	RefListIterator *iter;

	if ( skip_some_files(fh, config_startafter) ) {
		ERROR("Failed to skip first %i files.\n", config_startafter);
		return;
	}

	do {

		struct image image;

		image.det = NULL;

		/* Get data from next chunk */
		rval = read_chunk(fh, &image);
		if ( rval ) break;

		n_patterns++;

		if ( (image.reflections != NULL) && (image.indexed_cell) ) {

			/* Scale if requested */
			if ( config_scale ) {
				scale_intensities(model, image.reflections,
				                  sym);
			}

			merge_pattern(model, image.reflections, config_maxonly,
			              sym, hist_vals, hist_h, hist_k, hist_l,
			              hist_i, pass);

			n_used++;

		}

		free(image.filename);
		reflist_free(image.reflections);
		image_feature_list_free(image.features);
		cell_free(image.indexed_cell);

		progress_bar(n_patterns, n_total_patterns-config_startafter,
		             "Merging");

	} while ( rval == 0 );

	/* Divide by counts to get mean intensity if necessary */
	if ( (pass == 1) && !config_sum && !config_maxonly ) {

		Reflection *refl;
		RefListIterator *iter;

		for ( refl = first_refl(model, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			double intensity = get_intensity(refl);
			int red = get_redundancy(refl);

			set_int(refl, intensity / (double)red);

		}

	}

	/* Calculate ESDs */
	if ( pass == 2 ) {
		for ( refl = first_refl(model, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) ) {

			double sum_squared_dev = get_sum_squared_dev(refl);
			int red = get_redundancy(refl);
			int h, k, l;
			double esd;
			get_indices(refl,&h,&k,&l);

			/* Other ways of estimating the sigma are possible,
			 * such as:
			 *
			 *    double intensity = get_intensity(refl);
			 *    esd = sqrt( (sum_squared_dev/(double)red)
			 *              - pow(intensity,2.0) ) );
			 *
			 * But alter the other part of the calculation above
			 * as well. */
			esd = sqrt(sum_squared_dev)/(double)red;

			set_esd_intensity(refl, esd);

		}
	}

	if ( pass == 1 ) {
		STATUS("%i of the patterns could be used.\n", n_used);
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	FILE *fh;
	RefList *model;
	UnitCell *cell = NULL;
	int config_maxonly = 0;
	int config_startafter = 0;
	int config_stopafter = 0;
	int config_sum = 0;
	int config_scale = 0;
	unsigned int n_total_patterns;
	char *sym = NULL;
	char *pdb = NULL;
	char *histo = NULL;
	signed int hist_h, hist_k, hist_l;
	double *hist_vals = NULL;
	int hist_i;
	struct beam_params *beam = NULL;
	int space_for_hist = 0;

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
		{"beam",               1, NULL,               'b'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:o:p:y:g:f:b:",
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
		ERROR("Please specify filename using the -i option\n");
		return 1;
	}

	if ( output == NULL ) {
		output = strdup("processed.hkl");
	}

	if ( pdb != NULL ) {
		cell = load_cell_from_pdb(pdb);
		if ( cell == NULL ) {
			ERROR("Failed to load cell from '%s'\n", pdb);
			return 1;
		}
		free(pdb);
	} else {
		cell = NULL;
	}

	if ( sym == NULL ) sym = strdup("1");

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
	if ( n_total_patterns == 0 ) {
		ERROR("No patterns to process.\n");
		return 1;
	}
	STATUS("There are %i patterns to process\n", n_total_patterns);
	rewind(fh);

	model = reflist_new();

	if ( histo != NULL ) {
		int r;
		r = sscanf(histo, "%i,%i,%i", &hist_h, &hist_k, &hist_l);
		if ( r != 3 ) {
			ERROR("Invalid indices for '--histogram'\n");
			return 1;
		}
		space_for_hist = n_total_patterns * num_general_equivs(sym);
		hist_vals = malloc(space_for_hist * sizeof(double));
		free(histo);
		STATUS("Histogramming %i %i %i -> ", hist_h, hist_k, hist_l);
		/* Put into the asymmetric cell for the target group */
		get_asymm(hist_h, hist_k, hist_l,
		          &hist_h, &hist_k, &hist_l, sym);
		STATUS("%i %i %i\n", hist_h, hist_k, hist_l);
	}

	hist_i = 0;
	merge_all(fh, model, config_maxonly, config_scale, config_sum,
	          config_startafter, config_stopafter,
                  sym, n_total_patterns,
                  hist_vals, hist_h, hist_k, hist_l, &hist_i, 1);
	if ( ferror(fh) ) {
		ERROR("Stream read error.\n");
		return 1;
	}
	rewind(fh);
	if ( space_for_hist && (hist_i >= space_for_hist) ) {
		ERROR("Histogram array was too small!\n");
	}

	if ( hist_vals != NULL ) {
		STATUS("%i %i %i was seen %i times.\n", hist_h, hist_k, hist_l,
		                                        hist_i);
		plot_histogram(hist_vals, hist_i);
	}

	STATUS("Extra pass to calculate ESDs...\n");
	rewind(fh);
	merge_all(fh, model, config_maxonly, config_scale, 0,
	          config_startafter, config_stopafter, sym, n_total_patterns,
	          NULL, 0, 0, 0, NULL, 2);
	if ( ferror(fh) ) {
		ERROR("Stream read error.\n");
		return 1;
	}

	write_reflist(output, model, cell);

	fclose(fh);
	free(sym);
	reflist_free(model);
	free(output);
	if ( cell != NULL ) cell_free(cell);

	return 0;
}
