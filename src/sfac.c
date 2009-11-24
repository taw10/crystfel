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


#define N_MEMO 1024


/* Look up f1 and f2 for this atom at this energy (in J/photon) */
static double complex get_f1f2(const char *n, double en)
{
	FILE *fh;
	char filename[64];
	char line[1024];
	char *rval;
	double last_E, last_f1, last_f2;
	static char *memo_n[N_MEMO];
	static double memo_en[N_MEMO];
	static double complex memo_res[N_MEMO];
	static int n_memo = 0;
	int i;

	for ( i=0; i<n_memo; i++ ) {
		if ( (memo_en[i] == en) && (strcmp(memo_n[i], n) == 0) ) {
			return memo_res[i];
		}
	}

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
			double complex res;

			f = (en - last_E) / (E - last_E);

			actual_f1 = last_f1 + f * (f1 - last_f1);
			actual_f2 = last_f2 + f * (f2 - last_f2);

			fclose(fh);

			res = actual_f1 + I*actual_f2;

			memo_n[n_memo] = strdup(n);
			memo_en[n_memo] = en;
			memo_res[n_memo++] = res;
			n_memo = n_memo % N_MEMO;

			return res;

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
	static char *memo_n[N_MEMO];
	static double memo_s[N_MEMO];
	static double memo_res[N_MEMO];
	static int n_memo = 0;
	int i;

	for ( i=0; i<n_memo; i++ ) {
		if ( (memo_s[i] == s) && (strcmp(memo_n[i], n) == 0) ) {
			return memo_res[i];
		}
	}

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

	memo_n[n_memo] = strdup(n);
	memo_s[n_memo] = s;
	memo_res[n_memo++] = f;
	n_memo = n_memo % N_MEMO;

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


static void centre_molecule(struct molecule *mol)
{
	int i;
	double tx = 0.0;
	double ty = 0.0;
	double tz = 0.0;
	double total = 0.0;

	/* Atoms are grouped by species for faster calculation */
	for ( i=0; i<mol->n_species; i++ ) {

		struct mol_species *spec;
		int j;

		spec = mol->species[i];

		for ( j=0; j<spec->n_atoms; j++ ) {

			double sfac = get_waas_kirf(spec->species, 0.0);

			tx += spec->x[j] * sfac;
			ty += spec->y[j] * sfac;
			tz += spec->z[j] * sfac;
			total += sfac;

		}

	}

	mol->xc = tx / total;
	mol->yc = ty / total;
	mol->zc = tz / total;

	for ( i=0; i<mol->n_species; i++ ) {

		struct mol_species *spec;
		int j;

		spec = mol->species[i];

		for ( j=0; j<spec->n_atoms; j++ ) {

			spec->x[j] -= mol->xc;
			spec->y[j] -= mol->yc;
			spec->z[j] -= mol->zc;

		}

	}

	printf("Molecule was shifted by %f,%f,%f nm\n",
	       mol->xc*1e9, mol->yc*1e9, mol->zc*1e9);

}


/* Load PDB file into a memory format suitable for efficient(ish) structure
 * factor calculation */
struct molecule *load_molecule()
{
	struct molecule *mol;
	FILE *fh;
	char line[1024];
	char *rval;
	int i;

	mol = malloc(sizeof(struct molecule));
	if ( mol == NULL ) return NULL;
	mol->n_species = 0;
	mol->reflections = NULL;

	/* FIXME: Read cell from file */
	mol->cell = cell_new_from_parameters(28.10e-9,
	                                     28.10e-9,
	                                     16.52e-9,
	                                     deg2rad(90.0),
	                                     deg2rad(90.0),
	                                     deg2rad(120.0));

	fh = fopen("molecule.pdb", "r");
	if ( fh == NULL ) {
                fprintf(stderr, "Couldn't open file\n");
                return NULL;
        }

	do {

		char el[4];
		int j, r;
		int done = 0;
		float xf, yf, zf, occf, Bf;
		double x, y, z, occ, B;
		char *coords;

		rval = fgets(line, 1023, fh);

		/* Only interested in atoms */
		if ( strncmp(line, "HETATM", 6) != 0 ) continue;

		/* The following crimes against programming style
		 * were brought to you by Wizbit Enterprises, Inc. */
		if ( line[76] == ' ' ) {
			el[0] = line[77];
			el[1] = '\0';
		} else {
			el[0] = line[76];
			el[1] = line[77];
			el[2] = '\0';
		}

		coords = line + 29;
		r = sscanf(coords, "%f %f %f %f %f", &xf, &yf, &zf, &occf, &Bf);
		if ( r != 5 ) {
			fprintf(stderr, "WTF?\n");
			abort();
		}
		/* Promote to double precision */
		x = xf;  y = yf;  z = zf;  occ = occf;  B = Bf;

		for ( j=0; j<mol->n_species; j++ ) {

			struct mol_species *spec;
			int n;

			spec = mol->species[j];

			if ( strcmp(spec->species, el) != 0 ) continue;

			n = mol->species[j]->n_atoms;

			spec->x[n] = x*1.0e-10;  /* Convert to m */
			spec->y[n] = y*1.0e-10;
			spec->z[n] = z*1.0e-10;
			spec->occ[n] = occ;
			spec->B[n] = B*1.0e-20;  /* Convert to m^2 */
			mol->species[j]->n_atoms++;

			done = 1;

		}

		if ( !done ) {

			/* Need to create record for this species */
			struct mol_species *spec;

			spec = malloc(sizeof(struct mol_species));

			memcpy(spec->species, el, 4);
			spec->x[0] = x*1.0e-10;  /* Convert to nm */
			spec->y[0] = y*1.0e-10;
			spec->z[0] = z*1.0e-10;
			spec->occ[0] = occ;
			spec->B[0] = B*1.0e-20;  /* Convert to nm^2 */
			spec->n_atoms = 1;

			mol->species[mol->n_species] = spec;
			mol->n_species++;

		}


	} while ( rval != NULL );

	fclose(fh);

	centre_molecule(mol);

	printf("There are %i species\n", mol->n_species);
	for ( i=0; i<mol->n_species; i++ ) {
		printf("%3s : %6i\n", mol->species[i]->species,
		       mol->species[i]->n_atoms);
	}

	return mol;
}


double complex *get_reflections(struct molecule *mol, double en)
{
	double complex *reflections;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	signed int h, k, l;

	cell_get_reciprocal(mol->cell, &asx, &asy, &asz,
	                               &bsx, &bsy, &bsz,
	                               &csx, &csy, &csz);

	reflections = reflist_new();

	for ( h=-INDMAX; h<=INDMAX; h++ ) {
	for ( k=-INDMAX; k<=INDMAX; k++ ) {
	for ( l=-INDMAX; l<=INDMAX; l++ ) {

		double complex F;
		int i;
		double s;

		/* Atoms are grouped by species for faster calculation */
		for ( i=0; i<mol->n_species; i++ ) {

			double complex sfac;
			double complex contrib = 0.0;
			struct mol_species *spec;
			int j;

			spec = mol->species[i];

			for ( j=0; j<spec->n_atoms; j++ ) {

				double ph, u, v, w;

				u = h*asx + k*bsx + l*csx;
				v = h*asy + k*bsy + l*csy;
				w = h*asz + k*bsz + l*csz;

				ph = u*spec->x[j] + v*spec->y[j] + w*spec->z[j];

				/* Conversion from revolutions to radians */
				contrib += cos(-2.0*M_PI*ph)
				                          + I*sin(-2.0*M_PI*ph);

			}

			sfac = get_sfac(spec->species, s, en);
			F += sfac * contrib * exp(-2.0 * spec->B[j] * s);

		}

		integrate_reflection(reflections, h, k, l, F);

	}
	}
	progress_bar((h+INDMAX+1), 2*INDMAX);
	}
	printf("\n");

	return reflections;
}
