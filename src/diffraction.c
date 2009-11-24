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


static double lattice_factor(struct threevec q, double ax, double ay, double az,
                                                double bx, double by, double bz,
                                                double cx, double cy, double cz)
{
	struct threevec Udotq;
	double f1, f2, f3;
	int na = 4;
	int nb = 4;
	int nc = 30;

	Udotq.u = ax*q.u + ay*q.v + az*q.w;
	Udotq.v = bx*q.u + by*q.v + bz*q.w;
	Udotq.w = cx*q.u + cy*q.v + cz*q.w;

	/* At exact Bragg condition, f1 = na */
	if ( na > 1 ) {
		f1 = sin(M_PI*(double)na*Udotq.u) / sin(M_PI*Udotq.u);
	} else {
		f1 = 1.0;
	}

	/* At exact Bragg condition, f2 = nb */
	if ( nb > 1 ) {
		f2 = sin(M_PI*(double)nb*Udotq.v) / sin(M_PI*Udotq.v);
	} else {
		f2 = 1.0;
	}

	/* At exact Bragg condition, f3 = nc */
	if ( nc > 1 ) {
		f3 = sin(M_PI*(double)nc*Udotq.w) / sin(M_PI*Udotq.w);
	} else {
		f3 = 1.0;
	}

	/* At exact Bragg condition, this will multiply the molecular
	 * part of the structure factor by the number of unit cells,
	 * as desired (more scattering from bigger crystal!) */
	return f1 * f2 * f3;
}


/* Look up the structure factor for the nearest Bragg condition */
static double complex molecule_factor(struct molecule *mol, struct threevec q,
                                      double ax, double ay, double az,
                                      double bx, double by, double bz,
                                      double cx, double cy, double cz)
{
	double hd, kd, ld;
	signed int h, k, l;
	double complex r;

	hd = q.u * ax + q.v * ay + q.w * az;
	kd = q.u * bx + q.v * by + q.w * bz;
	ld = q.u * cx + q.v * cy + q.w * cz;
	h = (signed int)nearbyint(hd);
	k = (signed int)nearbyint(kd);
	l = (signed int)nearbyint(ld);

	r = get_integral(mol->reflections, h, k, l);

	return r;
}


static double complex water_factor(struct threevec q, double en)
{
	return 0.0;
}


void get_diffraction(struct image *image)
{
	int x, y;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	/* Generate the array of reciprocal space vectors in image->qvecs */
	get_ewald(image);

	if ( image->molecule == NULL ) {
		image->molecule = load_molecule();
		if ( image->molecule == NULL ) return;
	}

	cell_get_cartesian(image->molecule->cell, &ax, &ay, &az,
		                                  &bx, &by, &bz,
		                                  &cx, &cy, &cz);

	image->sfacs = malloc(image->width * image->height
	                      * sizeof(double complex));

	if ( image->molecule->reflections == NULL ) {
		get_reflections_cached(image->molecule, image->xray_energy);
	}

	progress_bar(0, image->width-1);
	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		double f_lattice;
		double complex f_molecule;
		double complex f_water;
		struct threevec q;
		double complex val;

		q = image->qvecs[x + image->width*y];

		f_lattice = lattice_factor(q, ax,ay,az,bx,by,bz,cx,cy,cz);
		f_molecule = molecule_factor(image->molecule, q,
		                             ax,ay,az,bx,by,bz,cx,cy,cz);

		f_water = water_factor(q, image->xray_energy);
		val = (f_molecule * f_lattice) + f_water;
		image->sfacs[x + image->width*y] = val;

	}
	progress_bar(x, image->width-1);
	}
}
