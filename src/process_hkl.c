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
"  -c, --compare-with=<file> Compare with reflection intensities in this file\n"
"\n"
"  -e, --output-every=<n>    Analyse figures of merit after every n patterns\n"
"                             Default: 1000.  A value of zero means to do the\n"
"                             analysis only after reading all the patterns.\n"
"  -r, --rvsq                Output lists of R vs |q| (\"Luzzatti plots\")\n"
"                             when analysing figures of merit.\n"
"      --zone-axis           Output an [001] zone axis pattern each time the\n"
"                             figures of merit are analysed.\n"
"      --detwin              Correlate each new pattern with the current\n"
"                             model and choose the best fitting out of the\n"
"                             allowable twins.\n"
"      --scale               Scale each pattern for best fit with the current\n"
"                             model.\n"
);
}


static void write_RvsQ(const char *name, double *ref, double *trueref,
                       unsigned int *counts, double scale, UnitCell *cell)
{
	FILE *fh;
	double smax, sbracket;
	signed int h, k, l;

	fh = fopen(name, "w");

	smax = 0.0;
	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {
		double s = 2.0*resolution(cell, h, k, l);
		if ( (lookup_count(counts, h, k, l) > 0) && (s > smax) ) {
			smax = s;
		}
	}
	}
	}

	for ( sbracket=0.0; sbracket<smax; sbracket+=smax/RVQDV ) {

		double top = 0.0;
		double bot = 0.0;
		int nhits = 0;
		int nrefl = 0;
		double R;
		double hits_per_refl;

		for ( h=-INDMAX; h<INDMAX; h++ ) {
		for ( k=-INDMAX; k<INDMAX; k++ ) {
		for ( l=-INDMAX; l<INDMAX; l++ ) {

			double s;
			int c;

			if ( (h==0) && (k==0) && (l==0) ) continue;

			c = lookup_count(counts, h, k, l);
			s = 2.0*resolution(cell, h, k, l);
			if ((s>=sbracket) && (s<sbracket+smax/RVQDV) && (c>0)) {

				double obs, calc, obsi;

				obs = lookup_intensity(ref, h, k, l);
				calc = lookup_intensity(trueref, h, k, l);

				obsi = obs / (double)c;
				top += pow(obsi - scale*calc, 2.0);
				bot += pow(obsi, 2.0);
				nhits += c;
				nrefl++;

			}

		}
		}
		}

		R = sqrt(top/bot);
		hits_per_refl = nrefl ? (double)nhits/nrefl : 0;

		fprintf(fh, "%8.5f %8.5f %5.2f\n", sbracket+smax/(2.0*RVQDV),
		                                R, hits_per_refl);

	}
	fclose(fh);
}


static void process_reflections(double *ref, unsigned int *counts,
                                double *trueref, unsigned int *truecounts,
                                unsigned int n_patterns,
                                UnitCell *cell, int do_rvsq, int do_zoneaxis)
{
	int j;
	double mean_counts;
	int ctot = 0;
	int nmeas = 0;
	double R, scale;
	FILE *fh;

	for ( j=0; j<LIST_SIZE; j++ ) {
		ctot += counts[j];
		if ( counts[j] > 0 ) nmeas++;
	}
	mean_counts = (double)ctot/nmeas;

	R = stat_r2(ref, counts, trueref, truecounts, &scale);
	STATUS("%8u: R=%5.2f%% (sf=%7.4e)  mean meas/refl=%5.2f,"
	       " %i reflections measured\n",
	       n_patterns, R*100.0, scale, mean_counts, nmeas);

	if ( do_rvsq ) {
		/* Record graph of R against q for this N */
		char name[64];
		snprintf(name, 63, "R_vs_q-%u.dat", n_patterns);
		write_RvsQ(name, ref, trueref, counts, scale, cell);
	}

	if ( do_zoneaxis ) {
		char name[64];
		snprintf(name, 63, "ZA-%u.dat", n_patterns);
		write_reflections(name, counts, ref, 1, cell, 1);
	}

	fh = fopen("results/convergence.dat", "a");
	fprintf(fh, "%u %5.2f %5.2f\n", n_patterns, R*100.0, mean_counts);
	fclose(fh);
}


static void merge_pattern(double *model, const double *new,
                          unsigned int *model_counts,
                          const unsigned int *counts, int mo, int sum)
{
	signed int h, k, l;

	for ( l=-INDMAX; l<INDMAX; l++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( h=-INDMAX; h<INDMAX; h++ ) {

		double intensity;

		if ( lookup_count(counts, h, k, l) == 0 ) continue;

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
	}
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	FILE *fh;
	unsigned int n_patterns;
	double *model, *trueref = NULL;
	unsigned int *model_counts;
	char *rval;
	UnitCell *cell;
	int config_maxonly = 0;
	int config_every = 1000;
	int config_rvsq = 0;
	int config_stopafter = 0;
	int config_zoneaxis = 0;
	int config_sum = 0;
	int config_detwin = 0;
	int config_scale = 0;
	char *intfile = NULL;
	double *new_pattern = NULL;
	unsigned int *new_counts = NULL;
	unsigned int n_total_patterns;
	unsigned int *truecounts = NULL;
	float f0;
	int f0_valid;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"max-only",           0, &config_maxonly,     1},
		{"output-every",       1, NULL,               'e'},
		{"rvsq",               0, NULL,               'r'},
		{"stop-after",         1, NULL,               's'},
		{"zone-axis",          0, &config_zoneaxis,    1},
		{"compare-with",       0, NULL,               'c'},
		{"sum",                0, &config_sum,         1},
		{"detwin",             0, &config_detwin,      1},
		{"scale",              0, &config_scale,      1},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:ro:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'i' : {
			filename = strdup(optarg);
			break;
		}

		case 'o' : {
			output = strdup(optarg);
			break;
		}

		case 'r' : {
			config_rvsq = 1;
			break;
		}

		case 'e' : {
			config_every = atoi(optarg);
			break;
		}

		case 's' : {
			config_stopafter = atoi(optarg);
			break;
		}

		case 'c' : {
			intfile = strdup(optarg);
			break;
		}

		case 0 : {
			break;
		}

		default : {
			return 1;
		}
		}

	}

	if ( filename == NULL ) {
		ERROR("Please specify filename using the -i option\n");
		return 1;
	}

	if ( intfile != NULL ) {
		truecounts = new_list_count();
		STATUS("Comparing against '%s'\n", intfile);
		trueref = read_reflections(intfile, truecounts);
		free(intfile);
	} else {
		trueref = NULL;
	}

	model = new_list_intensity();
	model_counts = new_list_count();
	cell = load_cell_from_pdb("molecule.pdb");
	new_pattern = new_list_intensity();
	new_counts = new_list_count();

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
	n_total_patterns = 0;
	do {
		char line[1024];

		rval = fgets(line, 1023, fh);
		if ( (strncmp(line, "Reflections from indexing", 25) == 0)
		    || (strncmp(line, "New pattern", 11) == 0) ) {
		    n_total_patterns++;
		}
	} while ( rval != NULL );
	rewind(fh);
	STATUS("There are %i patterns to process\n", n_total_patterns);

	n_patterns = 0;
	f0_valid = 0;
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

			/* Detwin before scaling */
			if ( config_detwin ) {
				detwin_intensities(model, new_pattern,
				                   model_counts, new_counts);
			}

			/* Assume a default I0 if we don't have one by now */
			if ( config_scale && !f0_valid ) {
				ERROR("No f0 value.\n");
				f0 = 1.0;
			}

			/* Scale if requested */
			if ( config_scale ) {
				scale_intensities(model, new_pattern,
				                  model_counts, new_counts, f0,
				                  f0_valid);
			}

			/* Start of second or later pattern */
			merge_pattern(model, new_pattern, model_counts,
			              new_counts, config_maxonly, config_sum);

			if ( (trueref != NULL) && config_every
			    && (n_patterns % config_every == 0) ) {
				process_reflections(model, model_counts,
				                    trueref, truecounts,
				                    n_patterns, cell,
				                    config_rvsq,
				                    config_zoneaxis);
			}

			if ( n_patterns == config_stopafter ) break;

			zero_list_count(new_counts);

			n_patterns++;

			progress_bar(n_patterns, n_total_patterns, "Merging");

			f0_valid = 0;

		}

		if ( strncmp(line, "f0 = ", 5) == 0 ) {
			r = sscanf(line, "f0 = %f", &f0);
			if ( r != 1 ) {
				ERROR("Couldn't understand f0 line.\n");
				f0 = 1.0;
				f0_valid = 0;
				continue;
			}
			f0_valid = 1;
		}

		r = sscanf(line, "%i %i %i %f", &h, &k, &l, &intensity);
		if ( r != 4 ) continue;

		if ( (h==0) && (k==0) && (l==0) ) continue;

		if ( lookup_count(new_counts, h, k, l) != 0 ) {
			ERROR("More than one measurement for %i %i %i in"
			      " pattern number %i\n", h, k, l, n_patterns);
		}
		set_intensity(new_pattern, h, k, l, intensity);
		set_count(new_counts, h, k, l, 1);

	} while ( rval != NULL );

	fclose(fh);

	if ( trueref != NULL ) {
		process_reflections(model, model_counts, trueref, truecounts,
		                    n_patterns, cell, config_rvsq,
		                    config_zoneaxis);
	}

	if ( output != NULL ) {
		write_reflections(output, model_counts, model, 0, cell, 1);
	}

	if ( config_zoneaxis ) {
		char name[64];
		snprintf(name, 63, "ZA-%u.dat", n_patterns);
		write_reflections(name, model_counts, model, 1, cell, 10);
	}

	STATUS("There were %u patterns.\n", n_patterns);

	return 0;
}
