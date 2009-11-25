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
	h = (signed int)rint(hd);
	k = (signed int)rint(kd);
	l = (signed int)rint(ld);

	r = lookup_sfac(mol->reflections, h, k, l);

	return r;
}


double water_intensity(struct threevec q, double en)
{
	double complex fH, fO;
	double s, modq;
	double intensity;

	/* Interatomic distances in water molecule */
	const double rOH = 0.09584e-9;
	const double rHH = 0.1515e-9;

	/* Dimensions of water column */
	const double water_r = 0.5e-6;
	const double beam_r = 1.5e-6;

	/* Volume of water column */
	const double water_v = M_PI*pow(water_r, 2.0) * 2.0 * beam_r;

	/* Number of water molecules */
	const double n_water = water_v * WATER_DENSITY
	                                        * (AVOGADRO / WATER_MOLAR_MASS);

	/* s = sin(theta)/lambda = 1/2d = |q|/2 */
	modq = modulus(q.u, q.v, q.w);
	s = modq / 2.0;

	fH = get_sfac("H", s, en);
	fO = get_sfac("O", s, en);

	/* Four O-H cross terms */
	intensity = 4.0*fH*fO * sin(modq*rOH)/(modq*rOH);

	/* Three H-H cross terms */
	intensity += 3.0*fH*fH * sin(modq*rHH)/(modq*rHH);

	/* Three diagonal terms */
	intensity += 2.0*fH*fH + fO*fO;

	return intensity * n_water;
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
		struct threevec q;
		double complex val;

		q = image->qvecs[x + image->width*y];

		f_lattice = lattice_factor(q, ax,ay,az,bx,by,bz,cx,cy,cz);
		f_molecule = molecule_factor(image->molecule, q,
		                             ax,ay,az,bx,by,bz,cx,cy,cz);

		val = f_molecule * f_lattice;
		image->sfacs[x + image->width*y] = val;

	}
	progress_bar(x, image->width-1);
	}
}
