/*
 * process_hkl.c
 *
 * Assemble and process FEL Bragg intensities
 *
 * (c) 2006-2009 Thomas White <thomas.white@desy.de>
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
"\n"
"      --max-only          Take the integrated intensity to be equal to the\n"
"                           maximum intensity measured for that reflection.\n"
"                           The default is to use the mean value from all\n"
"                           measurements.\n"
"  -e, --output-every=<n>  Analyse figures of merit after every n patterns.\n");
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


static void write_reflections(const char *filename, unsigned int *counts,
                              double *ref)
{
	FILE *fh;
	signed int h, k, l;

	fh = fopen(filename, "w");
	for ( h=-INDMAX; h<INDMAX; h++ ) {
	for ( k=-INDMAX; k<INDMAX; k++ ) {
	for ( l=-INDMAX; l<INDMAX; l++ ) {
		int N;
		N = lookup_count(counts, h, k, l);
		if ( N == 0 ) continue;
		double F = lookup_intensity(ref, h, k, l) / N;
		fprintf(fh, "%3i %3i %3i %f\n", h, k, l, F);
	}
	}
	}
	fclose(fh);
}


static double *ideal_intensities(double complex *sfac)
{
	double *ref;
	signed int h, k, l;

	ref = new_list_intensity();

	/* Generate ideal reflections from complex structure factors */
	for ( h=-INDMAX; h<=INDMAX; h++ ) {
	for ( k=-INDMAX; k<=INDMAX; k++ ) {
	for ( l=-INDMAX; l<=INDMAX; l++ ) {
		double complex F = lookup_sfac(sfac, h, k, l);
		double intensity = pow(cabs(F), 2.0);
		set_intensity(ref, h, k, l, intensity);
	}
	}
	}

	return ref;
}


static void process_reflections(double *ref, double *trueref,
                                unsigned int *counts, unsigned int n_patterns,
                                UnitCell *cell)
{
	int j;
	double mean_counts;
	int ctot = 0;
	int nmeas = 0;
	char name[64];
	double R, scale;
	double calc_222, obs_222;

	for ( j=0; j<LIST_SIZE; j++ ) {
		ctot += counts[j];
		if ( counts[j] > 0 ) nmeas++;
	}
	mean_counts = (double)ctot/nmeas;

	calc_222 = lookup_intensity(ref, 2, 2, 2) / lookup_count(counts, 2, 2, 2);
	obs_222 = lookup_intensity(trueref, 2, 2, 2);

	R = stat_r2(ref, trueref, counts, LIST_SIZE, &scale);
	STATUS("%8u: R=%5.2f%% (sf=%7.4e)  mean meas/refl=%5.2f,"
	       " %i reflections measured, %f\n",
	       n_patterns, R*100.0, scale, mean_counts, nmeas, calc_222/obs_222);

	/* Record graph of R against q for this N */
	snprintf(name, 63, "results/R_vs_q-%u.dat", n_patterns);
	write_RvsQ(name, ref, trueref, counts, scale, cell);
}


int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	FILE *fh;
	unsigned int n_patterns;
	double *ref, *trueref;
	unsigned int *counts;
	char *rval;
	struct molecule *mol;
	int config_maxonly = 0;
	int config_every = 1000;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{"max-only",           0, &config_maxonly,     1},
		{"output-every",       1, NULL,               'e'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:", longopts, NULL)) != -1) {

		switch (c) {
		case 'h' : {
			show_help(argv[0]);
			return 0;
		}

		case 'i' : {
			filename = strdup(optarg);
			break;
		}

		case 'e' : {
			config_every = atoi(optarg);
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

	if ( config_every <= 0 ) {
		ERROR("Invalid value for --output-every.\n");
		return 1;
	}

	mol = load_molecule();
	get_reflections_cached(mol, eV_to_J(2.0e3));

	ref = new_list_intensity();
	counts = new_list_count();
	trueref = ideal_intensities(mol->reflections);

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
		int h, k, l, intensity;
		int r;

		rval = fgets(line, 1023, fh);
		if ( strncmp(line, "New pattern", 11) == 0 ) {
			n_patterns++;
			if ( n_patterns % config_every == 0 ) {
				process_reflections(ref, trueref, counts,
				                    n_patterns, mol->cell);
			}
		}

		r = sscanf(line, "%i %i %i %i", &h, &k, &l, &intensity);
		if ( r != 4 ) continue;

		//if ( (abs(h)>3) || (abs(k)>3) || (abs(l)>3) ) continue;

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

	write_reflections("results/reflections.hkl", counts, ref);

	return 0;
}
