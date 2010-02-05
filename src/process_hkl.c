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


/* Number of divisions for R vs |q| graphs */
#define RVQDV (20)


static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Assemble and process FEL Bragg intensities.\n"
"\n"
"  -h, --help              Display this help message.\n"
"  -i, --input=<filename>  Specify input filename (\"-\" for stdin).\n"
"  -o, --output=<filename> Specify output filename for merged intensities\n"
"                           (don't specify for no output).\n"
"\n"
"      --max-only          Take the integrated intensity to be equal to the\n"
"                           maximum intensity measured for that reflection.\n"
"                           The default is to use the mean value from all\n"
"                           measurements.\n"
"  -e, --output-every=<n>  Analyse figures of merit after every n patterns\n"
"                           Default: 1000.  A value of zero means to do the\n"
"                           analysis only after reading all the patterns.\n"
"      --no-analyse        Don't perform any kind of analysis, just merge the\n"
"                           intensities.\n"
"  -r, --rvsq              Output lists of R vs |q| (\"Luzzatti plots\") when\n"
"                           analysing figures of merit.\n"
"      --stop-after=<n>    Stop after processing n patterns.  Zero means\n"
"                           keep going until the end of the input, and is the\n"
"                           default.\n"
"      --zone-axis         Output an [001] zone axis pattern each time the\n"
"                           figures of merit are analysed.\n");
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
		double s = resolution(cell, h, k, l);
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
			s = resolution(cell, h, k, l);
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


static void process_reflections(double *ref, double *trueref,
                                unsigned int *counts, unsigned int n_patterns,
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

	R = stat_r2(ref, trueref, counts, LIST_SIZE, &scale);
	STATUS("%8u: R=%5.2f%% (sf=%7.4e)  mean meas/refl=%5.2f,"
	       " %i reflections measured\n",
	       n_patterns, R*100.0, scale, mean_counts, nmeas);

	if ( do_rvsq ) {
		/* Record graph of R against q for this N */
		char name[64];
		snprintf(name, 63, "results/R_vs_q-%u.dat", n_patterns);
		write_RvsQ(name, ref, trueref, counts, scale, cell);
	}

	if ( do_zoneaxis ) {
		char name[64];
		snprintf(name, 63, "results/ZA-%u.dat", n_patterns);
		write_reflections(name, counts, ref, 1, cell);
	}

	fh = fopen("results/convergence.dat", "a");
	fprintf(fh, "%u %5.2f %5.2f\n", n_patterns, R*100.0, mean_counts);
	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	char *output = NULL;
	FILE *fh;
	unsigned int n_patterns;
	double *ref, *trueref = NULL;
	unsigned int *counts;
	char *rval;
	struct molecule *mol = NULL;
	int config_maxonly = 0;
	int config_every = 1000;
	int config_rvsq = 0;
	int config_stopafter = 0;
	int config_zoneaxis = 0;
	int config_noanalyse = 0;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"output",             1, NULL,               'o'},
		{"max-only",           0, &config_maxonly,     1},
		{"output-every",       1, NULL,               'e'},
		{"no-analyse",         0, &config_noanalyse,   1},
		{"rvsq",               0, NULL,               'r'},
		{"stop-after",         1, NULL,               's'},
		{"zone-axis",          0, &config_zoneaxis,    1},
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

	ref = new_list_intensity();
	counts = new_list_count();

	if ( !config_noanalyse ) {
		mol = load_molecule();
		get_reflections_cached(mol, eV_to_J(2.0e3));

		trueref = ideal_intensities(mol->reflections);
	}

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

	n_patterns = 0;
	do {

		char line[1024];
		signed int h, k, l, intensity;
		int r;

		rval = fgets(line, 1023, fh);
		if ( strncmp(line, "New pattern", 11) == 0 ) {

			if ( n_patterns == 0 ) {
				n_patterns++;
				continue;
			}

			if (config_every && (n_patterns % config_every == 0)) {
				process_reflections(ref, trueref, counts,
				                    n_patterns, mol->cell,
				                    config_rvsq,
				                    config_zoneaxis);
			}

			if ( n_patterns == config_stopafter ) break;

			n_patterns++;
		}

		r = sscanf(line, "%i %i %i %i", &h, &k, &l, &intensity);
		if ( r != 4 ) continue;

		if ( (h==0) && (k==0) && (l==0) ) continue;

		if ( !config_maxonly ) {
			integrate_intensity(ref, h, k, l, intensity);
			integrate_count(counts, h, k, l, 1);
		} else {
			if ( intensity > lookup_intensity(ref, h, k, l) ) {
				set_intensity(ref, h, k, l, intensity);
			}
			set_count(counts, h, k, l, 1);
		}

	} while ( rval != NULL );

	fclose(fh);

	if ( !config_noanalyse ) {
		process_reflections(ref, trueref, counts, n_patterns, mol->cell,
		                    config_rvsq, config_zoneaxis);
	}

	if ( output != NULL ) {
		write_reflections(output, counts, ref, 0, NULL);
	}

	STATUS("There were %u patterns.\n", n_patterns);

	return 0;
}
