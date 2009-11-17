/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "ewald.h"
#include "diffraction.h"
#include "sfac.h"


static double lattice_factor(struct threevec q, double ax, double ay, double az,
                                                double bx, double by, double bz,
                                                double cx, double cy, double cz)
{
	struct threevec Udotq;
	double f1, f2, f3;
	int na = 8;
	int nb = 8;
	int nc = 8;

	Udotq.u = (ax*q.u + ay*q.v + az*q.w)/2.0;
	Udotq.v = (bx*q.u + by*q.v + bz*q.w)/2.0;
	Udotq.w = (cx*q.u + cy*q.v + cz*q.w)/2.0;

	if ( na > 1 ) {
		f1 = sin(2.0*M_PI*(double)na*Udotq.u) / sin(2.0*M_PI*Udotq.u);
	} else {
		f1 = 1.0;
	}

	if ( nb > 1 ) {
		f2 = sin(2.0*M_PI*(double)nb*Udotq.v) / sin(2.0*M_PI*Udotq.v);
	} else {
		f2 = 1.0;
	}

	if ( nc > 1 ) {
		f3 = sin(2.0*M_PI*(double)nc*Udotq.w) / sin(2.0*M_PI*Udotq.w);
	} else {
		f3 = 1.0;
	}

	return f1 * f2 * f3;
}


/* Return structure factor for molecule 'mol' at energy 'en' (J/photon) at
 * scattering vector 'q' */
static double complex molecule_factor(struct molecule *mol, struct threevec q,
                                      double en)
{
	int i;
	double F = 0.0;
	double s;

	/* s = sin(theta)/lambda = 1/2d = (1/d)/2.0 */
	s = modulus(q.u, q.v, q.w) / 2.0;

	for ( i=0; i<mol->n_species; i++ ) {

		double complex sfac;
		double complex contrib = 0.0;
		struct mol_species *spec;
		int j;

		spec = mol->species[i];

		for ( j=0; j<spec->n_atoms; j++ ) {

			double ph;

			ph= q.u*spec->x[j] + q.v*spec->y[j] + q.w*spec->z[j];
			
			/* Conversion from revolutions to radians is required */
			contrib += cos(2.0*M_PI*ph) + I*sin(2.0*M_PI*ph);
			
		}

		sfac = get_sfac(spec->species, s, en);
		F += sfac * contrib * exp(-2.0 * spec->B[j] * s);
		
	}

	return F;
}


/* Load PDB file into a memory format suitable for efficient(ish) structure
 * factor calculation */
static struct molecule *load_molecule()
{
	struct molecule *mol;
	FILE *fh;
	char line[1024];
	char *rval;
	int i;

	mol = malloc(sizeof(struct molecule));
	if ( mol == NULL ) return NULL;
	mol->n_species = 0;

	fh = fopen("molecule.pdb", "r");
	if ( fh == NULL ) {
                fprintf(stderr, "Couldn't open file\n");
                return NULL;
        }

	do {

		char el[4];
		int j, r;
		int done = 0;
		float x, y, z, occ, B;
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
		r = sscanf(coords, "%f %f %f %f %f", &x, &y, &z, &occ, &B);
		if ( r != 5 ) {
			fprintf(stderr, "WTF?\n");
			abort();
		}

		for ( j=0; j<mol->n_species; j++ ) {

			struct mol_species *spec;
			int n;

			spec = mol->species[j];

			if ( strcmp(spec->species, el) != 0 ) continue;

			n = mol->species[j]->n_atoms;

			spec->x[n] = x;
			spec->y[n] = y;
			spec->z[n] = z;
			spec->occ[n] = occ;
			spec->B[n] = B;
			mol->species[j]->n_atoms++;

			done = 1;

		}

		if ( !done ) {

			/* Need to create record for this species */
			struct mol_species *spec;

			spec = malloc(sizeof(struct mol_species));

			memcpy(spec->species, el, 4);
			spec->x[0] = x;
			spec->y[0] = y;
			spec->z[0] = z;
			spec->occ[0] = occ;
			spec->B[0] = B;
			spec->n_atoms = 1;

			mol->species[mol->n_species] = spec;
			mol->n_species++;

		}


	} while ( rval != NULL );

	fclose(fh);

	printf("There are %i species\n", mol->n_species);
	for ( i=0; i<mol->n_species; i++ ) {
		printf("'%s': %i\n", mol->species[i]->species,
		       mol->species[i]->n_atoms);
	}

	return mol;
}


void get_diffraction(struct image *image, UnitCell *cell)
{
	int x, y;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	/* Generate the array of reciprocal space vectors in image->qvecs */
	get_ewald(image);
	image->molecule = load_molecule();
	if ( image->molecule == NULL ) return;

	cell_get_cartesian(cell, &ax, &ay, &az,
		                 &bx, &by, &bz,
		                 &cx, &cy, &cz);

	image->sfacs = malloc(image->width * image->height
	                      * sizeof(double complex));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double f_lattice;
		double complex f_molecule;
		struct threevec q;

		q = image->qvecs[x + image->width*y];

		f_lattice = lattice_factor(q, ax,ay,az,bx,by,bz,cx,cy,cz);
		f_molecule = molecule_factor(image->molecule, q,
		                             image->xray_energy);

		image->sfacs[x + image->width*y] = f_lattice * f_molecule;

	}
	printf("x=%i\n", x);
	}
}
