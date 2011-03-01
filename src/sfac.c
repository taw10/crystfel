/*
 * sfac.c
 *
 * Scattering factors
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <ctype.h>

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
	static int memo_eV[N_MEMO];
	static double complex memo_res[N_MEMO];
	static int n_memo = 0;
	int eV;
	int i;

	eV = (int)rint(J_to_eV(en));

	for ( i=0; i<n_memo; i++ ) {
		if ( (memo_eV[i] == eV) && (strcmp(memo_n[i], n) == 0) ) {
			return memo_res[i];
		}
	}

	snprintf(filename, 63, DATADIR"/crystfel/%s.nff", n);
	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Couldn't open file '%s'\n", filename);
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
			STATUS("I couldn't understand a line in the f1f2 "
			       "tables\n");
			continue;
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
			memo_eV[n_memo] = eV;
			memo_res[n_memo++] = res;
			n_memo = n_memo % N_MEMO;

			return res;

		}

	} while ( rval != NULL );

	fclose(fh);

	ERROR("Couldn't find scattering factors for '%s' at %f eV!\n", n, en);

	return 0.0;
}


struct waas_kirf_factors {
	char *n;
	float a1;
	float a2;
	float a3;
	float a4;
	float a5;
	float b1;
	float b2;
	float b3;
	float b4;
	float b5;
	float c;
};

/* s = sin(theta)/lambda in metres^-1*/
static double get_waas_kirf(const char *n, double s)
{
	FILE *fh;
	char *rval;
	double f;
	float a1, a2, a3, a4, a5, c, b1, b2, b3, b4, b5;
	double s2;
	int i;
	static struct waas_kirf_factors waaskirfcache[N_MEMO];
	int found;
	static int n_waaskirfcache = 0;

	found = 0;
	for ( i=0; i<n_waaskirfcache; i++ ) {
		if ( strcmp(waaskirfcache[i].n, n) == 0 ) {

			found = 1;

			a1 = waaskirfcache[n_waaskirfcache].a1;
			a2 = waaskirfcache[n_waaskirfcache].a2;
			a3 = waaskirfcache[n_waaskirfcache].a3;
			a4 = waaskirfcache[n_waaskirfcache].a4;
			a5 = waaskirfcache[n_waaskirfcache].a5;
			b1 = waaskirfcache[n_waaskirfcache].b1;
			b2 = waaskirfcache[n_waaskirfcache].b2;
			b3 = waaskirfcache[n_waaskirfcache].b3;
			b4 = waaskirfcache[n_waaskirfcache].b4;
			b5 = waaskirfcache[n_waaskirfcache].b5;
			c = waaskirfcache[n_waaskirfcache].c;

			break;
		}
	}

	if ( !found ) {

		fh = fopen(DATADIR"/crystfel/f0_WaasKirf.dat", "r");
		if ( fh == NULL ) {
			ERROR("Couldn't open f0_WaasKirf.dat\n");
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
			r = sscanf(line,
			             "  %f  %f  %f  %f  %f  %f  %f %f %f %f %f",
			            &a1, &a2, &a3, &a4, &a5, &c,
			            &b1, &b2, &b3, &b4, &b5);
			if ( r != 11 ) {
				ERROR("Couldn't read scattering "
				      "factors (WaasKirf)\n");
				return 0.0;
			}

			break;

		} while ( rval != NULL );

		fclose(fh);

		waaskirfcache[n_waaskirfcache].a1 = a1;
		waaskirfcache[n_waaskirfcache].a2 = a2;
		waaskirfcache[n_waaskirfcache].a3 = a3;
		waaskirfcache[n_waaskirfcache].a4 = a4;
		waaskirfcache[n_waaskirfcache].a5 = a5;
		waaskirfcache[n_waaskirfcache].c = c;
		waaskirfcache[n_waaskirfcache].b1 = b1;
		waaskirfcache[n_waaskirfcache].b2 = b2;
		waaskirfcache[n_waaskirfcache].b3 = b3;
		waaskirfcache[n_waaskirfcache].b4 = b4;
		waaskirfcache[n_waaskirfcache].b5 = b5;
		waaskirfcache[n_waaskirfcache++].n = strdup(n);
		n_waaskirfcache = n_waaskirfcache % N_MEMO;  /* Unlikely */

	}

	s2 = pow(s/1e10, 2.0);  /* s2 is s squared in Angstroms^-2 */
	f = c + a1*exp(-b1*s2) + a2*exp(-b2*s2) + a3*exp(-b3*s2)
	      + a4*exp(-b4*s2) + a5*exp(-b5*s2);

	return f;
}


/* Get complex scattering factors for element 'n' at energy 'en' (J/photon),
 * at resolution 's' = sin(theta)/lambda (in m^-1) */
static double complex get_sfac(const char *n, double s, double en)
{
	double complex f1f2;
	double fq;

	/* Anomalous part (point-like K-shell electrons only) */
	f1f2 = get_f1f2(n, en);

	/* Atomic form factor part (not complex-valued) */
	fq = get_waas_kirf(n, s) - get_waas_kirf(n, 0.0);

	return fq + f1f2;
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

	STATUS("The atoms were shifted by %5.3f, %5.3f, %5.3f nm\n",
	       mol->xc*1e9, mol->yc*1e9, mol->zc*1e9);
}


/* Load PDB file into a memory format suitable for efficient(ish) structure
 * factor calculation */
struct molecule *load_molecule(const char *filename)
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
	mol->cell = NULL;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
                ERROR("Couldn't open PDB file\n");
                return NULL;
        }

	do {

		int j, r;
		int done = 0;
		char xs[8];
		char ys[8];
		char zs[8];
		char occs[6];
		char Bs[6];
		char el[3];
		float xf, yf, zf, occf, Bf;
		double x, y, z, occ, B;

		rval = fgets(line, 1023, fh);

		/* Only interested in atoms */
		if ( (strncmp(line, "HETATM", 6) != 0)
		  && (strncmp(line, "ATOM  ", 6) != 0) ) continue;

		chomp(line);
		if ( strlen(line) != 80 ) {
			STATUS("Line does not have the correct length (%i):\n",
			       (int)strlen(line));
			STATUS("'%s'\n", line);
			continue;
		}

		/* Separate out fixed-width fields */
		memcpy(xs, line+30, 5);    xs[7] = '\0';
		memcpy(ys, line+38, 5);    ys[7] = '\0';
		memcpy(zs, line+46, 5);    zs[7] = '\0';
		memcpy(occs, line+54, 5);  occs[5] = '\0';
		memcpy(Bs, line+60, 5);    Bs[5] = '\0';
		memcpy(el, line+76, 2);    el[2] = '\0';

		/* Convert fields */
		r = sscanf(xs, "%f", &xf);
		r += sscanf(ys, "%f", &yf);
		r += sscanf(zs, "%f", &zf);
		r += sscanf(occs, "%f", &occf);
		r += sscanf(Bs, "%f", &Bf);
		if ( el[0] == ' ' ) {
			el[0] = el[1];
			el[1] = '\0';
		} else {
			el[1] = tolower(el[1]);
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
			if ( mol->species[j]->n_atoms == MAX_ATOMS ) {
				ERROR("Too many atoms of type '%s'!\n", el);
				exit(1);
			}

			done = 1;

		}

		if ( !done ) {

			/* Need to create record for this species */
			struct mol_species *spec;

			spec = malloc(sizeof(struct mol_species));

			memcpy(spec->species, el, 1+strlen(el));
			spec->x[0] = x*1.0e-10;  /* Convert to m */
			spec->y[0] = y*1.0e-10;
			spec->z[0] = z*1.0e-10;
			spec->occ[0] = occ;
			spec->B[0] = B*1.0e-20;  /* Convert to m^2 */
			spec->n_atoms = 1;

			mol->species[mol->n_species] = spec;
			mol->n_species++;

		}

	} while ( rval != NULL );

	fclose(fh);

	centre_molecule(mol);

	STATUS("There are %i species:\n", mol->n_species); fflush(stderr);
	for ( i=0; i<mol->n_species; i++ ) {
		STATUS("%3s : %6i\n", mol->species[i]->species,
		       mol->species[i]->n_atoms);
	}

	mol->cell = load_cell_from_pdb(filename);
	if ( mol->cell == NULL ) {
		ERROR("No unit cell found in PDB file\n");
		return NULL;
	}

	return mol;
}


void free_molecule(struct molecule *mol)
{
	int i;

	for ( i=0; i<mol->n_species; i++ ) {
		free(mol->species[i]);
	}

	free(mol->cell);
	free(mol);
}


double *get_reflections(struct molecule *mol, double en, double res,
                        double *phases, ReflItemList *items)
{
	double *reflections;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	signed int h, k, l;
	const int do_thermal = 1;

	cell_get_reciprocal(mol->cell, &asx, &asy, &asz,
	                               &bsx, &bsy, &bsz,
	                               &csx, &csy, &csz);

	reflections = new_list_intensity();

	for ( h=-INDMAX; h<=INDMAX; h++ ) {
	for ( k=-INDMAX; k<=INDMAX; k++ ) {
	for ( l=-INDMAX; l<=INDMAX; l++ ) {

		double complex F = 0.0;
		int i;
		double s, oneoverd;

		/* We need sin(theta)/lambda = 1/2d */
		s = resolution(mol->cell, h, k, l);
		oneoverd = 2.0 * s;
		if ( oneoverd > res ) continue;

		/* Atoms are grouped by species for faster calculation */
		for ( i=0; i<mol->n_species; i++ ) {

			double complex sfac;
			double complex contrib = 0.0;
			struct mol_species *spec;
			int j;

			spec = mol->species[i];

			for ( j=0; j<spec->n_atoms; j++ ) {

				double ph, u, v, w;
				double complex cpart;

				u = h*asx + k*bsx + l*csx;
				v = h*asy + k*bsy + l*csy;
				w = h*asz + k*bsz + l*csz;

				ph = u*spec->x[j] + v*spec->y[j] + w*spec->z[j];

				/* Conversion from revolutions to radians */
				cpart = cexp(-2.0*M_PI*I*ph);
				if ( do_thermal ) {
					cpart *= exp(-2.0 * spec->B[j] * s * s);
				}
				contrib += cpart;

			}

			sfac = get_sfac(spec->species, s, en);
			F += sfac * contrib;

		}

		set_intensity(reflections, h, k, l, pow(cabs(F), 2.0));
		if ( phases != NULL ) set_phase(phases, h, k, l, carg(F));
		if ( items != NULL ) add_item(items, h, k, l);

	}
	progress_bar((k+INDMAX)+IDIM*(h+INDMAX), IDIM*IDIM-1,
	             "Calculating reflection intensities");
	}
	}

	return reflections;
}
