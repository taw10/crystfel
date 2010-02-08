/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
#include "diffraction.h"
#include "sfac.h"


#define SAMPLING (2)
#define BWSAMPLING (10)
#define BANDWIDTH (0.015)


static double lattice_factor(struct rvec q, double ax, double ay, double az,
                                            double bx, double by, double bz,
                                            double cx, double cy, double cz,
                                            int na, int nb, int nc)
{
	struct rvec Udotq;
	double f1, f2, f3;

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
static double complex molecule_factor(struct molecule *mol, struct rvec q,
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


double water_intensity(struct rvec q, double en,
                       double beam_r, double water_r)
{
	double complex fH, fO;
	double s, modq;
	double width;
	double complex ifac;

	/* Interatomic distances in water molecule */
	const double rOH = 0.09584e-9;
	const double rHH = 0.1515e-9;

	/* Volume of water column, approximated as:
	 * (2water_r) * (2beam_r) * smallest(2beam_r, 2water_r)
	 * neglecting the curvature of the faces of the volume */
	if ( beam_r > water_r ) {
		width = 2.0 * water_r;
	} else {
		width = 2.0 * beam_r;
	}
	const double water_v = 2.0*beam_r * 2.0*water_r * width;

	/* Number of water molecules */
	const double n_water = water_v * WATER_DENSITY
	                                        * (AVOGADRO / WATER_MOLAR_MASS);

	/* s = sin(theta)/lambda = 1/2d = |q|/2 */
	modq = modulus(q.u, q.v, q.w);
	s = modq / 2.0;

	fH = get_sfac("H", s, en);
	fO = get_sfac("O", s, en);

	/* Four O-H cross terms */
	ifac = 4.0*fH*fO * sin(2.0*M_PI*modq*rOH)/(2.0*M_PI*modq*rOH);

	/* Three H-H cross terms */
	ifac += 3.0*fH*fH * sin(2.0*M_PI*modq*rHH)/(2.0*M_PI*modq*rHH);

	/* Three diagonal terms */
	ifac += 2.0*fH*fH + fO*fO;

	return cabs(ifac) * n_water;
}


struct rvec get_q(struct image *image, unsigned int xs, unsigned int ys,
                  unsigned int sampling, float *ttp, float k)
{
	struct rvec q;
	float twothetax, twothetay, twotheta, r;
	float rx = 0.0;
	float ry = 0.0;
	int p;

	const unsigned int x = xs / sampling;
	const unsigned int y = ys / sampling; /* Integer part only */

	for ( p=0; p<image->det.n_panels; p++ ) {
		if ( (x >= image->det.panels[p].min_x)
		  && (x <= image->det.panels[p].max_x)
		  && (y >= image->det.panels[p].min_y)
		  && (y <= image->det.panels[p].max_y) ) {
			rx = ((float)xs - (sampling*image->det.panels[p].cx))
			               / (sampling * image->resolution);
			ry = ((float)ys - (sampling*image->det.panels[p].cy))
			               / (sampling * image->resolution);
			break;
		}
	}

	/* Calculate q-vector for this sub-pixel */
	r = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
	twothetax = atan2(rx, image->camera_len);
	twothetay = atan2(ry, image->camera_len);
	twotheta = atan2(r, image->camera_len);

	if ( ttp != NULL ) *ttp = twotheta;

	q.u = k * sin(twothetax);
	q.v = k * sin(twothetay);
	q.w = k - k * cos(twotheta);

	return quat_rot(q, image->orientation);
}


void get_diffraction(struct image *image, int na, int nb, int nc, int no_sfac)
{
	unsigned int xs, ys;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double a, b, c, d;
	float kc;

	if ( image->molecule == NULL ) return;

	cell_get_cartesian(image->molecule->cell, &ax, &ay, &az,
		                                  &bx, &by, &bz,
		                                  &cx, &cy, &cz);

	cell_get_parameters(image->molecule->cell,
	                    &a, &b, &c, &d, &d, &d);
	STATUS("Particle size = %i x %i x %i (=%5.2f x %5.2f x %5.2f nm)\n",
	       na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

	/* Allocate (and zero) the "diffraction array" */
	image->sfacs = calloc(image->width * image->height,
	                      sizeof(double complex));

	if ( !no_sfac ) {
		if ( image->molecule->reflections == NULL ) {
			get_reflections_cached(image->molecule,
			                       ph_lambda_to_en(image->lambda));
		}
	}

	/* Needed later for Lorentz calculation */
	image->twotheta = malloc(image->width * image->height * sizeof(double));

	kc = 1.0/image->lambda;  /* Centre value */

	for ( xs=0; xs<image->width*SAMPLING; xs++ ) {
	for ( ys=0; ys<image->height*SAMPLING; ys++ ) {

		double f_lattice;
		double complex f_molecule;
		struct rvec q;
		float twotheta;
		double sw = 1.0/(SAMPLING*SAMPLING);  /* Sample weight */

		const unsigned int x = xs / SAMPLING;
		const unsigned int y = ys / SAMPLING; /* Integer part only */

		int kstep;

		for ( kstep=0; kstep<BWSAMPLING; kstep++ ) {

			float k;

			/* Calculate k this time round */
			k = kc + (kstep-(BWSAMPLING/2)) *
			                              kc*(BANDWIDTH/BWSAMPLING);

			q = get_q(image, xs, ys, SAMPLING, &twotheta, k);
			image->twotheta[x + image->width*y] = twotheta;

			f_lattice = lattice_factor(q, ax, ay, az,
			                              bx, by, bz,
			                              cx, cy, cz,
			                              na, nb, nc);
			if ( no_sfac ) {
				f_molecule = 10000.0;
			} else {
				f_molecule = molecule_factor(image->molecule, q,
			                            ax,ay,az,bx,by,bz,cx,cy,cz);
			}

			image->sfacs[x + image->width*y] += sw * f_molecule * f_lattice;

		}

	}
	progress_bar(xs, SAMPLING*image->width-1, "Calculating lattice factors");
	}
}
