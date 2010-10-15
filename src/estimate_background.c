/*
 * estimate-background.c
 *
 * Like 'process_hkl', but for signal/noise ratios
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


/* Number of bins for plot of resolution shells */
#define NBINS (10)



static void show_help(const char *s)
{
	printf("Syntax: %s [options]\n\n", s);
	printf(
"Estimate peak SNR from FEL Bragg intensities.\n"
"\n"
"  -h, --help                Display this help message.\n"
"  -i, --input=<filename>    Specify input filename (\"-\" for stdin).\n"
);
}



int main(int argc, char *argv[])
{
	int c;
	char *filename = NULL;
	FILE *fh;
	int rval;
	double total_vol, vol_per_shell;
	double rmin, rmax;
	double rmins[NBINS];
	double rmaxs[NBINS];
	double snrs[NBINS];
	unsigned int counts[NBINS];
	UnitCell *real_cell;
	double overall_snr;
	unsigned int overall_counts;
	int i;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"input",              1, NULL,               'i'},
		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "hi:e:ro:p:y:g:f:a:r:",
	                        longopts, NULL)) != -1) {

		switch (c) {
		case 'h' :
			show_help(argv[0]);
			return 0;

		case 'i' :
			filename = strdup(optarg);
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

	/* FIXME: Fixed resolution shells */
	rmin = 0.120e9;
	rmax = 1.172e9;

	for ( i=0; i<NBINS; i++ ) {
		counts[i] = 0;
		snrs[i] = 0.0;
	}

	total_vol = pow(rmax, 3.0) - pow(rmin, 3.0);
	vol_per_shell = total_vol / NBINS;
	rmins[0] = rmin;
	for ( i=1; i<NBINS; i++ ) {

		double r;

		r = vol_per_shell + pow(rmins[i-1], 3.0);
		r = pow(r, 1.0/3.0);

		/* Shells of constant volume */
		rmaxs[i-1] = r;
		rmins[i] = r;

		/* Shells of constant thickness */
		//rmins[i] = rmins[i-1] + (rmax-rmin)/NBINS;
		//rmaxs[i-1] = rmins[i-1] + (rmax-rmin)/NBINS;

		STATUS("Shell %i: %f to %f\n", i-1,
		       rmins[i-1]/1e9, rmaxs[i-1]/1e9);

	}
	rmaxs[NBINS-1] = rmax;
	STATUS("Shell %i: %f to %f\n", NBINS-1,
	       rmins[NBINS-1]/1e9, rmaxs[NBINS-1]/1e9);

	real_cell = load_cell_from_pdb("../../1JB0.pdb");

	do {

		char *rval2;
		UnitCell *cell;
		char *filename;
		char line[1024];
		int done = 0;

		rval = find_chunk(fh, &cell, &filename);
		if ( rval != 0 ) break;

		do {

			signed int h, k, l;
			float intensity, x, y, max, bg;
			int r, j;
			double d, snr;
			int bin;

			rval2 = fgets(line, 1024, fh);

			if ( strncmp(line, "Peak statistics:", 16) == 0 ) {
				done = 1;
			}

			r = sscanf(line, "%i %i %i %f (at %f,%f) max=%f bg=%f",
			           &h, &k, &l, &intensity, &x, &y, &max, &bg);

			if ( r != 8 ) {
				continue;
			}

			d = resolution(real_cell, h, k, l) * 2.0;
			//STATUS("'%s'", line);
			//STATUS("-> %i %i %i %f\n", h, k, l, d/1e9);

			bin = -1;
			for ( j=0; j<NBINS; j++ ) {
				if ( (d>rmins[j]) && (d<=rmaxs[j]) ) {
					bin = j;
					break;
				}
			}
			if ( bin == -1 ) {
				ERROR("Warnung! %i %i %i %f\n", h, k, l, d/1e9);
				continue;
			}

			snr = max/bg;
			snrs[bin] += snr;
			counts[bin]++;

		} while ( !done );

	} while ( !rval );

	/* Print out results */
	overall_snr = 0.0;
	overall_counts = 0;
	for ( i=0; i<NBINS; i++ ) {

		double snr, cen;

		cen = rmins[i] + (rmaxs[i] - rmins[i])/2.0;
		snr = snrs[i] / (double)counts[i];
		printf("%f %f (%f / %i)\n", cen*1.0e-9, snr, snrs[i], counts[i]);
		if ( i>0 ) {
			overall_snr += snrs[i];
			overall_counts += counts[i];
		}

	}

	printf("Overall: %f (%f / %i)\n", overall_snr/(double)overall_counts,
	                                  overall_snr, overall_counts);


	fclose(fh);

	return 0;
}
