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


/* Density of water in kg/m^3 */
#define WATER_DENSITY (1.0e6)

/* Molar mass of water, in kg/mol */
#define WATER_MOLAR_MASS (18.01528e3)

/* Avogadro's number */
#define AVOGADRO (6.022e23)

/* Number of slices to divide water column into, where the
 * slices are coplanar with the beam and the cylinder axis */
#define N_WATER_SLICES 100


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
	double complex F = 0.0;
	double s;

	/* s = sin(theta)/lambda = 1/2d = (1/d)/2.0 */
	s = modulus(q.u, q.v, q.w) / 2.0;

	/* Atoms are grouped by species for faster calculation */
	for ( i=0; i<mol->n_species; i++ ) {

		double complex sfac;
		double complex contrib = 0.0;
		struct mol_species *spec;
		int j;

		spec = mol->species[i];

		for ( j=0; j<spec->n_atoms; j++ ) {

			double ph;

			ph = q.u*spec->x[j] + q.v*spec->y[j] + q.w*spec->z[j];

			/* Conversion from revolutions to radians is required */
			contrib += cos(2.0*M_PI*ph) + I*sin(2.0*M_PI*ph);

		}

		sfac = get_sfac(spec->species, s, en);
		F += sfac * contrib * exp(-2.0 * spec->B[j] * s);

	}

	return F;
}


static double complex water_factor(struct threevec q, double en)
{
	double x;
	double complex res = 0.0;
	int n = 0;
	double s;
	double molecules_per_m3;
	double molecules_per_m;
	const double rc = 0.5e-6; /* Radius of cylinder */
	const double rb = 1.5e-6; /* Radius of beam */

	/* Density of water molecules */
	molecules_per_m3 = WATER_DENSITY * (AVOGADRO / WATER_MOLAR_MASS);
	
	/* Number of water molecules per slice */
	molecules_per_m = (2*rb*2*rc) * molecules_per_m3 / N_WATER_SLICES;

	/* s = sin(theta)/lambda = 1/2d = (1/d)/2.0 */
	s = modulus(q.u, q.v, q.w) / 2.0;

	for ( x=-rc; x<rc; x+=rc/(N_WATER_SLICES/2) ) {

		double ph;
		double thickness;

		/* The thickness of a cylinder in a direction parallel
		 * to the beam.  The cylinder axis is perpendicular to
		 * the beam */
		thickness = 2.0 * sqrt((rc*rc)-(x*x));

		ph = q.u*x;
		res += thickness * (cos(2.0*M_PI*ph) + I*sin(2.0*M_PI*ph));
		n++;

	}
	res *= (get_sfac("O", s, en) + 2.0*get_sfac("H", s, en))
	       * molecules_per_m;

	/* Correct for sampling of integration */
	return (res/n) * (2*rc);
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

	progress_bar(0, image->width-1);
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double f_lattice;
		double complex f_molecule;
		double complex f_water;
		struct threevec q;

		q = image->qvecs[x + image->width*y];

		f_lattice = lattice_factor(q, ax,ay,az,bx,by,bz,cx,cy,cz);
		f_molecule = molecule_factor(image->molecule, q,
		                             image->xray_energy);

		/* Nasty approximation follows */
		if ( y == image->height/2 ) {
			f_water = water_factor(q, image->xray_energy);
		} else {
			f_water = 0.0;
		}

		image->sfacs[x + image->width*y] = (f_lattice * f_molecule)
		                                                  + f_water;

	}
	progress_bar(x, image->width-1);
	}
}
