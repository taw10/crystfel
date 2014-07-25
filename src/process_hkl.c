/*
 * process_hkl.c
 *
 * Assemble and process FEL Bragg intensities
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2014 Thomas White <taw@physics.org>
 *   2011      Andrew Martin <andrew.martin@desy.de>
 *   2012      Lorenzo Galli <lorenzo.galli@desy.de>
 *   2014      Chunhong Yoon <chun.hong.yoon@desy.de>
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
#include "statistics.h"
#include "reflist-utils.h"
#include "symmetry.h"
#include "stream.h"
#include "reflist.h"
#include "image.h"
#include "crystal.h"
#include "thread-pool.h"
#include "geometry.h"
#include "cell-utils.h"


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Assemble and process FEL Bragg intensities.\n"
"\n"
"  -h, --help                Display this help message.\n"
"      --version             Print CrystFEL version number and exit.\n"
"  -i, --input=<filename>    Specify input filename (\"-\" for stdin).\n"
"  -o, --output=<filename>   Specify output filename for merged intensities\n"
"                             Default: processed.hkl).\n"
"  -y, --symmetry=<sym>      Merge according to point group <sym>.\n"
"\n"
"      --start-after=<n>     Skip <n> crystals at the start of the stream.\n"
"      --stop-after=<n>      Stop after merging <n> crystals.\n"
"  -g, --histogram=<h,k,l>   Calculate the histogram of measurements for this\n"
"                             reflection.\n"
"  -z, --hist-parameters     Set the range for the histogram and the number of\n"
"          =<min,max,nbins>   bins. \n"
"\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model.\n"
"      --no-polarisation     Disable polarisation correction.\n"
"      --min-measurements=<n> Require at least <n> measurements before a\n"
"                             reflection appears in the output.  Default: 2\n"
"      --min-snr=<n>         Require individual intensity measurements to\n"
"                             have I > n * sigma(I).  Default: -infinity.\n"
"      --max-adu=<n>         Maximum peak value.  Default: infinity.\n"
"      --min-res=<n>         Merge only crystals which diffract above <n> A.\n"
"      --push-res=<n>        Integrate higher than apparent resolution cutoff.\n"
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


static double *check_hist_size(int n, double *hist_vals)
{
	int ns;
	double *tryMe;

	if ( n % 1000 ) return hist_vals;

	ns = n / 1000;
	ns = (ns+1)*1000;

	tryMe = realloc(hist_vals, ns*sizeof(double));
	if ( tryMe == NULL ) {
		ERROR("Failed to allocate space for histogram.\n");
	}
	return tryMe;
}


static int merge_crystal(RefList *model, struct image *image, Crystal *cr,
                         RefList *reference, const SymOpList *sym,
                         double **hist_vals, signed int hist_h,
                         signed int hist_k, signed int hist_l, int *hist_n,
                         int config_nopolar, double min_snr, double max_adu,
                         double push_res)
{
	Reflection *refl;
	RefListIterator *iter;
	double scale;

	/* First, correct for polarisation */
	if ( !config_nopolar ) {
		polarisation_correction(crystal_get_reflections(cr),
		                        crystal_get_cell(cr), image);
	}

	if ( reference != NULL ) {
		scale = scale_intensities(reference,
		                          crystal_get_reflections(cr), sym);
	} else {
		scale = 1.0;
	}
	if ( isnan(scale) ) return 1;

	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double refl_intensity, refl_sigma, refl_pk;
		signed int h, k, l;
		int model_redundancy;
		Reflection *model_version;
		double w;
		double temp, delta, R, mean, M2, sumweight;
		double res, max_res;

		refl_intensity = scale * get_intensity(refl);
		refl_sigma = scale * get_esd_intensity(refl);
		refl_pk = get_peak(refl);
		w = 1.0;//pow(refl_sigma, -2.0);

		if ( (min_snr > -INFINITY) && isnan(refl_sigma) ) continue;
		if ( refl_intensity < min_snr * refl_sigma ) continue;

		if ( refl_pk > max_adu ) continue;

		get_indices(refl, &h, &k, &l);

		max_res = push_res + crystal_get_resolution_limit(cr);
		res = 2.0*resolution(crystal_get_cell(cr), h, k, l);
		if ( res > max_res ) continue;

		/* Put into the asymmetric unit for the target group */
		get_asymm(sym, h, k, l, &h, &k, &l);

		model_version = find_refl(model, h, k, l);
		if ( model_version == NULL ) {
			model_version = add_refl(model, h, k, l);
		}

		mean = get_intensity(model_version);
		sumweight = get_temp1(model_version);
		M2 = get_temp2(model_version);

		temp = w + sumweight;
		delta = refl_intensity - mean;
		R = delta * w / temp;
		set_intensity(model_version, mean + R);
		set_temp2(model_version, M2 + sumweight * delta * R);
		set_temp1(model_version, temp);

		model_redundancy = get_redundancy(model_version);
		set_redundancy(model_version, ++model_redundancy);

		if ( *hist_vals != NULL ) {

			if ( (h==hist_h) && (k==hist_k) && (l==hist_l) ) {

				*hist_vals = check_hist_size(*hist_n,
							     *hist_vals);

				/* Check again because realloc might have
				 * failed */
				if ( *hist_vals != NULL ) {
					(*hist_vals)[*hist_n] = refl_intensity;
					*hist_n += 1;
				}

			}

		}

	}

	return 0;
}


static void display_progress(int n_images, int n_crystals, int n_crystals_used)
{
	if ( !isatty(STDERR_FILENO) ) return;
	if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

	pthread_mutex_lock(&stderr_lock);
	fprintf(stderr, "\r%i images processed, %i crystals, %i crystals used.",
	        n_images, n_crystals, n_crystals_used);
	pthread_mutex_unlock(&stderr_lock);

	fflush(stdout);
}



static int merge_all(Stream *st, RefList *model, RefList *reference,
                     const SymOpList *sym,
                     double **hist_vals, signed int hist_h,
                     signed int hist_k, signed int hist_l,
                     int *hist_i, int config_nopolar, int min_measurements,
                     double min_snr, double max_adu,
                     int start_after, int stop_after, double min_res,
                     double push_res)
{
	int rval;
	int n_images = 0;
	int n_crystals = 0;
	int n_crystals_used = 0;
	Reflection *refl;
	RefListIterator *iter;
	int n_crystals_seen = 0;

	do {

		struct image image;
		int i;

		image.det = NULL;

		/* Get data from next chunk */
		rval = read_chunk_2(st, &image, STREAM_READ_REFLECTIONS |
		                    STREAM_READ_UNITCELL);
		if ( rval ) break;

		n_images++;

		for ( i=0; i<image.n_crystals; i++ ) {

			int r;
			Crystal *cr = image.crystals[i];

			n_crystals_seen++;
			if ( n_crystals_seen <= start_after ) continue;

			if ( crystal_get_resolution_limit(cr) < min_res ) {
				continue;
			}

			n_crystals++;
			r = merge_crystal(model, &image, cr, reference, sym,
			                  hist_vals, hist_h, hist_k, hist_l,
			                  hist_i, config_nopolar, min_snr,
			                  max_adu, push_res);

			if ( r == 0 ) n_crystals_used++;

			reflist_free(crystal_get_reflections(cr));
			cell_free(crystal_get_cell(cr));
			crystal_free(cr);

			if ( n_crystals_used == stop_after ) break;

		}

		free(image.filename);
		image_feature_list_free(image.features);
		free(image.crystals);

		display_progress(n_images, n_crystals_seen, n_crystals_used);

		if ( (stop_after>0) && (n_crystals_used == stop_after) ) break;

	} while ( rval == 0 );

	for ( refl = first_refl(model, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double var;
		int red;

		red = get_redundancy(refl);
		if ( red < min_measurements ) {
			set_redundancy(refl, 0);
			continue;
		}

		var = get_temp2(refl) / get_temp1(refl);
		set_esd_intensity(refl, sqrt(var)/sqrt(red));
	}

	return 0;
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	Stream *st;
	RefList *model;
	int config_scale = 0;
	char *sym_str = NULL;
	SymOpList *sym;
	char *histo = NULL;
	signed int hist_h, hist_k, hist_l;
	signed int hist_nbins=50;
	float hist_min=0.0, hist_max=0.0;
	double *hist_vals = NULL;
	int hist_i;
	int space_for_hist = 0;
	char *histo_params = NULL;
	int config_nopolar = 0;
	char *rval;
	int min_measurements = 2;
	int r;
	int start_after = 0;
	int stop_after = 0;
	double min_snr = -INFINITY;
	double max_adu = +INFINITY;
	double min_res = 0.0;
	double push_res = +INFINITY;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"start-after",        1, NULL,               's'},
		{"stop-after",         1, NULL,               'f'},
		{"scale",              0, &config_scale,       1},
		{"no-polarisation",    0, &config_nopolar,     1},
		{"no-polarization",    0, &config_nopolar,     1},
		{"symmetry",           1, NULL,               'y'},
		{"histogram",          1, NULL,               'g'},
		{"hist-parameters",    1, NULL,               'z'},
		{"min-measurements",   1, NULL,                2},
		{"min-snr",            1, NULL,                3},
		{"max-adu",            1, NULL,                4},
		{"min-res",            1, NULL,                5},
		{"push-res",           1, NULL,                6},
		{"res-push",           1, NULL,                6}, /* compat */
		{"version",            0, NULL,                7},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:o:y:g:s:f:z:",
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
			errno = 0;
			start_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --start-after.\n");
				return 1;
			}
			break;

			case 'f' :
			errno = 0;
			stop_after = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --stop-after.\n");
				return 1;
			}
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

			case 2 :
			errno = 0;
			min_measurements = strtol(optarg, &rval, 10);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-measurements.\n");
				return 1;
			}
			break;

			case 3 :
			errno = 0;
			min_snr = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-snr.\n");
				return 1;
			}
			ERROR("WARNING: You have used --min-snr.\n");
			ERROR("WARNING: Please read the manual carefully to "
			      "learn about possible detrimental effects of this"
			      " option.\n");
			break;

			case 4 :
			errno = 0;
			max_adu = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --max-adu.\n");
				return 1;
			}
			break;

			case 5 :
			errno = 0;
			min_res = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --min-res.\n");
				return 1;
			}
			min_res = 1e10/min_res;
			break;

			case 6 :
			errno = 0;
			push_res = strtod(optarg, &rval);
			if ( *rval != '\0' ) {
				ERROR("Invalid value for --push-res.\n");
				return 1;
			}
			push_res = push_res*1e9;
			break;

			case 7 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case '?' :
			break;

			case 0 :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

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
	st = open_stream_for_read(filename);
	if ( st == NULL ) {
		ERROR("Failed to open stream.\n");
		return 1;
	}

	model = reflist_new();

	if ( histo != NULL ) {

		int r;

		r = sscanf(histo, "%i,%i,%i", &hist_h, &hist_k, &hist_l);
		if ( r != 3 ) {
			ERROR("Invalid indices for '--histogram'\n");
			return 1;
		}
		free(histo);

		/* Allocate enough space that hist_vals isn't NULL.
		 * check_hist_size will realloc it straight away */
		hist_vals = malloc(1*sizeof(double));
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
	r = merge_all(st, model, NULL, sym, &hist_vals, hist_h, hist_k, hist_l,
	              &hist_i, config_nopolar, min_measurements, min_snr,
	              max_adu, start_after, stop_after, min_res, push_res);
	fprintf(stderr, "\n");
	if ( r ) {
		ERROR("Error while reading stream.\n");
		return 1;
	}

	if ( config_scale ) {

		RefList *reference;

		if ( rewind_stream(st) ) {

			ERROR("Couldn't rewind stream - scaling cannot be "
			      "performed.\n");

		} else {

			int r;

			STATUS("Extra pass for scaling...\n");

			reference = model;
			model = reflist_new();

			free(hist_vals);
			hist_vals = malloc(1*sizeof(double));
			hist_i = 0;

			r = merge_all(st, model, reference, sym,
				     &hist_vals, hist_h, hist_k, hist_l, &hist_i,
				     config_nopolar, min_measurements, min_snr,
				     max_adu, start_after, stop_after, min_res,
				     push_res);
			fprintf(stderr, "\n");
			if ( r ) {
				ERROR("Error while reading stream.\n");
				return 1;
			}

			reflist_free(reference);

		}

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

	write_reflist_2(output, model, sym);

	close_stream(st);

	free_symoplist(sym);
	reflist_free(model);
	free(output);
	free(filename);

	return 0;
}
