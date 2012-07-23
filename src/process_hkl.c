/*
 * process_hkl.c
 *
 * Assemble and process FEL Bragg intensities
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *   2011      Andrew Martin <andrew.martin@desy.de>
 *   2012      Lorenzo Galli <lorenzo.galli@desy.de>
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

#include "utils.h"
#include "statistics.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "stream.h"
#include "reflist.h"
#include "image.h"


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
"  -y, --symmetry=<sym>      Merge according to point group <sym>.\n"
"\n"
"      --start-after=<n>     Skip n patterns at the start of the stream.\n"
"      --stop-after=<n>      Stop after processing n patterns.\n"
"  -g, --histogram=<h,k,l>   Calculate the histogram of measurements for this\n"
"                             reflection.\n"
"  -z, --hist-parameters     Set the range for the histogram and the number of\n"
"          =<min,max,nbins>   bins. \n"
"\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model.\n"
"      --reference=<file>    Compare against intensities from <file> when\n"
"                             scaling. \n"
"      --no-polarisation     Disable polarisation correction.\n"
);
}


static void plot_histogram(double *vals, int n, float hist_min, float hist_max,
                           int nbins)
{
	int i;
	double max = -INFINITY;
	double min = +INFINITY;
	double step;
	int histo[nbins];
	FILE *fh;

	fh = fopen("histogram.dat", "w");
	if ( fh == NULL ) {
		ERROR("Couldn't open 'histogram.dat'\n");
		return;
	}

	if ( hist_min == hist_max ) {
		for ( i=0; i<n; i++ ) {
			if ( vals[i] > max ) max = vals[i];
			if ( vals[i] < min ) min = vals[i];
		}
	} else {
		min = hist_min;
		max = hist_max;
	}
	STATUS("min max nbins: %f %f %i\n", min, max, nbins);
	min--;  max++;

	for ( i=0; i<nbins; i++ ) {
		histo[i] = 0;
	}

	step = (max-min)/nbins;

	for ( i=0; i<n; i++ ) {
		int bin;
		if ( (vals[i] > min) && (vals[i] < max) ) {
			bin = (vals[i]-min)/step;
			histo[bin]++;
		}
	}

	for ( i=0; i<nbins; i++ ) {
		fprintf(fh, "%f %i\n", min+step*i, histo[i]);
	}

	fclose(fh);
}


static double scale_intensities(RefList *reference, RefList *new,
                              const SymOpList *sym)
{
	double s;
	double top = 0.0;
	double bot = 0.0;
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(new, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		double i1, i2;
		signed int hu, ku, lu;
		signed int h, k, l;
		Reflection *reference_version;

		get_indices(refl, &h, &k, &l);
		get_asymm(sym, h, k, l, &hu, &ku, &lu);

		reference_version = find_refl(reference, hu, ku, lu);
		if ( reference_version == NULL ) continue;

		i1 = get_intensity(reference_version);
		i2 = get_intensity(refl);

		/* Calculate LSQ estimate of scaling factor */
		top += i1 * i2;
		bot += i2 * i2;

	}

	s = top / bot;

	return s;
}


static int merge_pattern(RefList *model, struct image *new, RefList *reference,
                         const SymOpList *sym,
                         double *hist_vals, signed int hist_h,
                         signed int hist_k, signed int hist_l, int *hist_n,
                         int config_nopolar)
{
	Reflection *refl;
	RefListIterator *iter;
	double scale;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	if ( reference != NULL ) {
		scale = scale_intensities(reference, new->reflections, sym);
	} else {
		scale = 1.0;
	}
	if ( isnan(scale) ) return 1;

	cell_get_reciprocal(new->indexed_cell, &asx, &asy, &asz,
	                                       &bsx, &bsy, &bsz,
	                                       &csx, &csy, &csz);

	for ( refl = first_refl(new->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double intensity;
		double xl, yl, zl;
		double pol, pa, pb, phi, tt, ool;
		signed int h, k, l;
		int cur_redundancy;
		double cur_intensity, cur_sumsq;
		Reflection *model_version;

		get_indices(refl, &h, &k, &l);

		/* Put into the asymmetric unit for the target group */
		get_asymm(sym, h, k, l, &h, &k, &l);

		model_version = find_refl(model, h, k, l);
		if ( model_version == NULL ) {
			model_version = add_refl(model, h, k, l);
		}

		intensity = scale * get_intensity(refl);

		if ( !config_nopolar ) {

			/* Polarisation correction assuming 100% polarisation
			 * along the x direction */
			xl = h*asx + k*bsx + l*csx;
			yl = h*asy + k*bsy + l*csy;
			zl = h*asz + k*bsz + l*csz;

			ool = 1.0 / new->lambda;
			tt = angle_between(0.0, 0.0, 1.0,  xl, yl, zl+ool);
			phi = atan2(yl, xl);
			pa = pow(sin(phi)*sin(tt), 2.0);
			pb = pow(cos(tt), 2.0);
			pol = 1.0 - 2.0*(1.0-pa) + (1.0+pb);
			intensity /= pol;

		}

		cur_intensity = get_intensity(model_version);
		set_intensity(model_version, cur_intensity + intensity);

		cur_redundancy = get_redundancy(model_version);
		set_redundancy(model_version, cur_redundancy+1);

		cur_sumsq = get_temp1(model_version);
		set_temp1(model_version, cur_sumsq + pow(intensity, 2.0));

		if ( hist_vals != NULL ) {

			if ( (h==hist_h) && (k==hist_k) && (l==hist_l) ) {
				hist_vals[*hist_n] = intensity;
				*hist_n += 1;
			}

		}

	}

	return 0;
}


static void merge_all(FILE *fh, RefList *model, RefList *reference,
                      int config_startafter, int config_stopafter,
                      const SymOpList *sym,
                      int n_total_patterns,
                      double *hist_vals, signed int hist_h,
                      signed int hist_k, signed int hist_l,
                      int *hist_i, int config_nopolar)
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

			int r;

			r = merge_pattern(model, &image, reference, sym,
			                  hist_vals, hist_h, hist_k, hist_l,
			                  hist_i, config_nopolar);

			if ( r == 0 ) n_used++;

		}

		free(image.filename);
		reflist_free(image.reflections);
		image_feature_list_free(image.features);
		cell_free(image.indexed_cell);

		progress_bar(n_patterns, n_total_patterns-config_startafter,
		             "Merging");

		if ( config_stopafter ) {
			if ( n_patterns == config_stopafter ) break;
		}

	} while ( rval == 0 );

	for ( refl = first_refl(model, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double intensity, sumsq, esd;
		int red;

		red = get_redundancy(refl);
		if ( red == 1 ) {
			set_redundancy(refl, 0);
			continue;
		}

		intensity = get_intensity(refl) / red;
		set_intensity(refl, intensity);

		sumsq = get_temp1(refl) / red;
		esd = sqrt(sumsq - pow(intensity, 2.0)) / sqrt(red);
		set_esd_intensity(refl, esd);

	}

	STATUS("%i of the patterns could be used.\n", n_used);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	FILE *fh;
	RefList *model;
	int config_maxonly = 0;
	int config_startafter = 0;
	int config_stopafter = 0;
	int config_sum = 0;
	int config_scale = 0;
	unsigned int n_total_patterns;
	char *sym_str = NULL;
	SymOpList *sym;
	char *pdb = NULL;
	char *histo = NULL;
	signed int hist_h, hist_k, hist_l;
	signed int hist_nbins=50;
	float hist_min=0.0, hist_max=0.0;
	double *hist_vals = NULL;
	int hist_i;
	int space_for_hist = 0;
	char *histo_params = NULL;
	int config_nopolar = 0;

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
		{"no-polarisation",    0, &config_nopolar,     1},
		{"no-polarization",    0, &config_nopolar,     1},
		{"symmetry",           1, NULL,               'y'},
		{"histogram",          1, NULL,               'g'},
		{"hist-parameters",    1, NULL,               'z'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:o:p:y:g:f:b:z:",
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
			sym_str = strdup(optarg);
			break;

		case 'g' :
			histo = strdup(optarg);
			break;

		case 'z' :
			histo_params = strdup(optarg);
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

	if ( sym_str == NULL ) sym_str = strdup("1");
	sym = get_pointgroup(sym_str);
	free(sym_str);

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

		space_for_hist = n_total_patterns * num_equivs(sym, NULL);
		hist_vals = malloc(space_for_hist * sizeof(double));
		free(histo);
		STATUS("Histogramming %i %i %i -> ", hist_h, hist_k, hist_l);

		/* Put into the asymmetric cell for the target group */
		get_asymm(sym, hist_h, hist_k, hist_l,
		          &hist_h, &hist_k, &hist_l);
		STATUS("%i %i %i\n", hist_h, hist_k, hist_l);

	}

	if ( histo_params != NULL ) {

		int rr;

		rr = sscanf(histo_params, "%f,%f,%i", &hist_min, &hist_max,
		                                      &hist_nbins);
		if ( rr != 3 ) {
			ERROR("Invalid parameters for '--hist-parameters'\n");
			return 1;
		}
		free(histo_params);
		if ( hist_max <= hist_min ) {
			ERROR("Invalid range for '--hist-parameters'. "
			      "Make sure that 'max' is greater than 'min'.\n");
			return 1;
		}

	}

	hist_i = 0;
	merge_all(fh, model, NULL, config_startafter, config_stopafter,
	          sym, n_total_patterns, hist_vals, hist_h, hist_k, hist_l,
	          &hist_i, config_nopolar);
	if ( ferror(fh) ) {
		ERROR("Stream read error.\n");
		return 1;
	}
	rewind(fh);

	if ( config_scale ) {

		RefList *reference;

		STATUS("Extra pass for scaling...\n");

		reference = copy_reflist(model);

		reflist_free(model);
		model = reflist_new();

		rewind(fh);

		merge_all(fh, model, reference,
			  config_startafter, config_stopafter, sym,
			  n_total_patterns,
			  hist_vals, hist_h, hist_k, hist_l, &hist_i,
			  config_nopolar);

		if ( ferror(fh) ) {
			ERROR("Stream read error.\n");
			return 1;
		}

		reflist_free(reference);

	}

	if ( space_for_hist && (hist_i >= space_for_hist) ) {
		ERROR("Histogram array was too small!\n");
	}

	if ( hist_vals != NULL ) {
		STATUS("%i %i %i was seen %i times.\n", hist_h, hist_k, hist_l,
		                                        hist_i);
		plot_histogram(hist_vals, hist_i, hist_min, hist_max,
		               hist_nbins);
	}

	write_reflist(output, model);

	fclose(fh);

	free(sym);
	reflist_free(model);
	free(output);

	return 0;
}
