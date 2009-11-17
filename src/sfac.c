/*
 * sfac.c
 *
 * Scattering factors
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "utils.h"
#include "sfac.h"


/* Look up f1 and f2 for this atom at this energy (in J/photon) */
static double complex get_f1f2(const char *n, double en)
{
	FILE *fh;
	char filename[64];
	char line[1024];
	char *rval;
	double last_E, last_f1, last_f2;

	snprintf(filename, 63, "scattering-factors/%s.nff", n);
	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		fprintf(stderr, "Couldn't open file '%s'\n", filename);
		return 0.0;
	}

	en = J_to_eV(en);

	/* Discard first line */
	fgets(line, 1023, fh);

	last_E = 0.0;
	last_f1 = 0.0;
	last_f2 = 0.0;
	do {

		int r;
		double E, f1, f2;
		float E_f, f1_f, f2_f;

		rval = fgets(line, 1023, fh);

		r = sscanf(line, "%f %f %f", &E_f, &f1_f, &f2_f);
		if ( r != 3 ) {
			fprintf(stderr, "WTF?\n");
			abort();
		}
		/* Promote to double precision */
		E = E_f;  f1 = f1_f;  f2 = f2_f;

		/* Find the first energy greater than the required value */
		if ( E < en ) {
			/* Store old values ready for interpolation*/
			last_E = E;
			last_f1 = f1;
			last_f2 = f2;
		} else {

			/* Perform (linear) interpolation */
			double f;
			double actual_f1, actual_f2;

			f = (en - last_E) / (E - last_E);

			actual_f1 = last_f1 + f * (f1 - last_f1);
			actual_f2 = last_f2 + f * (f2 - last_f2);

			fclose(fh);
			return actual_f1 + I*actual_f2;

		}

	} while ( rval != NULL );

	fclose(fh);

	fprintf(stderr, "Couldn't find scattering factors for '%s' at %f eV!\n",
	        n, en);

	return 0.0;
}


/* s = sin(theta)/lambda */
static double get_waas_kirf(const char *n, double s)
{
	FILE *fh;
	char *rval;
	double f;
	float a1, a2, a3, a4, a5, c, b1, b2, b3, b4, b5;
	double s2;
	
	fh = fopen("scattering-factors/f0_WaasKirf.dat", "r");
	if ( fh == NULL ) {
		fprintf(stderr, "Couldn't open f0_WaasKirf.dat\n");
		return 0.0;
	}

	do {

		int r;
		char line[1024];
		char sp[1024];
		int Z;

		rval = fgets(line, 1023, fh);
		
		if ( (line[0] != '#') || (line[1] != 'S') ) continue;

		r = sscanf(line, "#S  %i  %s", &Z, sp);
		if ( (r != 2) || (strcmp(sp, n) != 0) ) continue;
		
		/* Skip two lines */
		fgets(line, 1023, fh);
		fgets(line, 1023, fh);
		
		/* Read scattering coefficients */
		rval = fgets(line, 1023, fh);
		r = sscanf(line, "  %f  %f  %f  %f  %f  %f  %f %f %f %f %f",
		          &a1, &a2, &a3, &a4, &a5, &c, &b1, &b2, &b3, &b4, &b5);
		if ( r != 11 ) {
			fprintf(stderr, "Couldn't read scattering factors\n");
			return 0.0;
		}
		
		break;
		
	} while ( rval != NULL );
	
	fclose(fh);
	
	s2 = pow(s/1e10, 2.0);  /* s2 is s squared in Angstroms squared */
	f = c + a1*exp(-b1*s2) + a2*exp(-b2*s2) + a3*exp(-b3*s2)
	      + a4*exp(-b4*s2) + a5*exp(-b5*s2);

	return f;
}


/* Get complex scattering factors for element 'n' at energy 'en' (J/photon),
 * at resolution 's' = sin(theta)/lambda (in m^-1) */
double complex get_sfac(const char *n, double s, double en)
{
	double complex f1f2;
	double fq, fq0;
	
	f1f2 = get_f1f2(n, en);
	fq = get_waas_kirf(n, s);
	fq0 = get_waas_kirf(n, 0.0);
	
	return fq - fq0 + f1f2;
}
