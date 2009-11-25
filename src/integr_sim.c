/*
 * integr_sim.c
 *
 * Test integration of intensities
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

#include "cell.h"
#include "image.h"
#include "utils.h"
#include "statistics.h"
#include "sfac.h"


static void main_show_help(const char *s)
{
	printf("Syntax: %s\n\n", s);
	printf("Test relrod integration\n\n");
	printf("  -h              Display this help message\n");
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

	for ( sbracket=0.0; sbracket<smax; sbracket+=smax/10.0 ) {

		double top = 0.0;
		double bot = 0.0;
		int n = 0;

		for ( h=-INDMAX; h<INDMAX; h++ ) {
		for ( k=-INDMAX; k<INDMAX; k++ ) {
		for ( l=-INDMAX; l<INDMAX; l++ ) {

			double s;
			int c;
			c = lookup_count(counts, h, k, l);
			s = resolution(cell, h, k, l);
			if ((s>=sbracket) && (s<sbracket+smax/10.0) && (c>0)) {

				double obs, calc, obsi;

				obs = lookup_intensity(ref, h, k, l);
				calc = lookup_intensity(trueref, h, k, l);

				obsi = obs / (double)c;
				top += fabs(obsi/scale - calc);
				bot += obsi/scale;
				n++;
			}

		}
		}
		}

		fprintf(fh, "%8.5f %8.5f %i\n", sbracket+smax/20.0, top/bot, n);

	}
	fclose(fh);
}


int main(int argc, char *argv[])
{
	int c;
	int i;
	struct image image;
	double *ref, *trueref;
	unsigned int *counts;
	size_t ref_size;
	signed int h, k, l;
	double R, scale;
	FILE *fh;

	while ((c = getopt(argc, argv, "h")) != -1) {

		switch ( c ) {

			case 'h' : {
				main_show_help(argv[0]);
				return 0;
			}

		}

	}

	/* Define image parameters */
	image.width = 1024;
	image.height = 1024;
	image.fmode = FORMULATION_CLEN;
	image.x_centre = 512.5;
	image.y_centre = 512.5;
	image.camera_len = 0.05;  /* 5 cm (front CCD can move from 5cm-20cm) */
	image.resolution = 13333.3; /* 75 micron pixel size */
	image.xray_energy = eV_to_J(2.0e3); /* 2 keV energy */
	image.lambda = ph_en_to_lambda(image.xray_energy);  /* Wavelength */
	image.molecule = load_molecule();

	/* Prepare array for integrated intensities */
	ref = new_list_intensity();

	/* Array for sample counts */
	counts = new_list_count();

	/* Calculate true intensities */
	get_reflections_cached(image.molecule, image.xray_energy);
	/* Complex structure factors now in image.molecule->reflections */

	for ( i=1; i<=10e3; i++ ) {

		#if 0
		image.orientation = random_quaternion();

		/* Calculate reflections using large smax
		 * (rather than the actual value) */
		//get_reflections(&image, cell, 1.0e9);

		nrefl = image_feature_count(image.rflist);
		for ( j=0; j<nrefl; j++ ) {

			struct imagefeature *f;
			double t, s, intensity, F;

			f = image_get_feature(image.rflist, j);

			t = 100.0e-9;  /* Thickness 100 nm */
			s = f->s;      /* Get excitation error */
			F = structure_factor(f->h, f->k, f->l);

			/* Calculate intensity from this reflection */
			intensity = pow( F * SINC(M_PI*t*s), 2.0);

			if ( intensity < 0.1 ) continue;

			if ( (f->h == 2) && (f->k == 2) && (f->l == 2) ) {
				fprintf(fh1, "%f %f\n", s, intensity);
			}

			if ( (f->h == 15) && (f->k == 15) && (f->l == 15) ) {
				fprintf(fh2, "%f %f\n", s, intensity);
			}

			integrate_reflection(ref, f->h, f->k, f->l, intensity);
			add_count(counts, f->h, f->k, f->l, 1);
		}
		#endif

		if ( i % 1000 == 0 ) {

			int j;
			double mean_counts;
			int ctot = 0;
			int nmeas = 0;
			double ff;
			char name[64];

			for ( j=0; j<ref_size; j++ ) {
				ctot += counts[j];
				if ( counts[j] > 0 ) nmeas++;
			}
			mean_counts = (double)ctot/nmeas;

			ff = lookup_intensity(ref, 2, 2, 2)
			                        / lookup_count(counts, 2, 2, 2);

			R = stat_r2(ref, trueref, counts, ref_size, &scale);
			printf("%8i: R=%5.2f%% (sf=%7.4f)"
			       " mean meas/refl=%5.2f,"
			       " %i reflections measured. %f\n",
			       i, R*100.0, scale, mean_counts, nmeas,
			       ff);

			/* Record graph of R against q for this N */
			snprintf(name, 63, "R_vs_q-%i.dat", i);
			write_RvsQ(name, ref, trueref, counts,
			           scale, image.molecule->cell);

		}

	}

	fh = fopen("reflections.dat", "w");
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

	return 0;
}
