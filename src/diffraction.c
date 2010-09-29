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
#include <assert.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "diffraction.h"
#include "sfac.h"
#include "parameters-lcls.tmp"


#define SAMPLING (4)
#define BWSAMPLING (1)
#define BANDWIDTH (0.0 / 100.0)


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


static double interpolate_linear(const double *ref,
                                 float hd, signed int k, signed int l)
{
	signed int h;
	double val1, val2;
	float f;

	h = (signed int)hd;
	if ( hd < 0.0 ) h -= 1;
	f = hd - (float)h;
	assert(f >= 0.0);

	val1 = lookup_intensity(ref, h, k, l);
	val2 = lookup_intensity(ref, h+1, k, l);

	val1 = val1;
	val2 = val2;

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_bilinear(const double *ref,
                                   float hd, float kd, signed int l)
{
	signed int k;
	double val1, val2;
	float f;

	k = (signed int)kd;
	if ( kd < 0.0 ) k -= 1;
	f = kd - (float)k;
	assert(f >= 0.0);

	val1 = interpolate_linear(ref, hd, k, l);
	val2 = interpolate_linear(ref, hd, k+1, l);

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_intensity(const double *ref,
                                    float hd, float kd, float ld)
{
	signed int l;
	double val1, val2;
	float f;

	l = (signed int)ld;
	if ( ld < 0.0 ) l -= 1;
	f = ld - (float)l;
	assert(f >= 0.0);

	val1 = interpolate_bilinear(ref, hd, kd, l);
	val2 = interpolate_bilinear(ref, hd, kd, l+1);

	return (1.0-f)*val1 + f*val2;
}


static double complex interpolate_phased_linear(const double *ref,
                                                const double *phases,
                                                float hd,
                                                signed int k, signed int l)
{
	signed int h;
	double val1, val2;
	float f;
	double ph1, ph2;
	double re1, re2, im1, im2;
	double re, im;

	h = (signed int)hd;
	if ( hd < 0.0 ) h -= 1;
	f = hd - (float)h;
	assert(f >= 0.0);

	val1 = lookup_intensity(ref, h, k, l);
	val2 = lookup_intensity(ref, h+1, k, l);
	ph1 = lookup_phase(phases, h, k, l);
	ph2 = lookup_phase(phases, h+1, k, l);

	val1 = val1;
	val2 = val2;

	/* Calculate real and imaginary parts */
	re1 = val1 * cos(ph1);
	im1 = val1 * sin(ph1);
	re2 = val2 * cos(ph2);
	im2 = val2 * sin(ph2);

	re = (1.0-f)*re1 + f*re2;
	im = (1.0-f)*im1 + f*im2;

	return re + im*I;
}


static double complex interpolate_phased_bilinear(const double *ref,
                                                  const double *phases,
                                                  float hd, float kd,
                                                  signed int l)
{
	signed int k;
	double complex val1, val2;
	float f;

	k = (signed int)kd;
	if ( kd < 0.0 ) k -= 1;
	f = kd - (float)k;
	assert(f >= 0.0);

	val1 = interpolate_phased_linear(ref, phases, hd, k, l);
	val2 = interpolate_phased_linear(ref, phases, hd, k+1, l);

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_phased_intensity(const double *ref,
                                           const double *phases,
                                           float hd, float kd, float ld)
{
	signed int l;
	double complex val1, val2;
	float f;

	l = (signed int)ld;
	if ( ld < 0.0 ) l -= 1;
	f = ld - (float)l;
	assert(f >= 0.0);

	val1 = interpolate_phased_bilinear(ref, phases, hd, kd, l);
	val2 = interpolate_phased_bilinear(ref, phases, hd, kd, l+1);

	return cabs((1.0-f)*val1 + f*val2);
}


/* Look up the structure factor for the nearest Bragg condition */
static double molecule_factor(const double *intensities,const double *phases,
                              struct rvec q,
                              double ax, double ay, double az,
                              double bx, double by, double bz,
                              double cx, double cy, double cz,
                              GradientMethod m)
{
	float hd, kd, ld;
	signed int h, k, l;
	double r;

	hd = q.u * ax + q.v * ay + q.w * az;
	kd = q.u * bx + q.v * by + q.w * bz;
	ld = q.u * cx + q.v * cy + q.w * cz;

	switch ( m ) {
	case GRADIENT_MOSAIC :
		h = (signed int)rint(hd);
		k = (signed int)rint(kd);
		l = (signed int)rint(ld);
		r = lookup_intensity(intensities, h, k, l);
		break;
	case GRADIENT_INTERPOLATE :
		r = interpolate_intensity(intensities, hd, kd, ld);
		break;
	case GRADIENT_PHASED :
		r = interpolate_phased_intensity(intensities, phases,
		                                 hd, kd, ld);
		break;
	default:
		ERROR("This gradient method not implemented yet.\n");
		exit(1);
	}

	return r;
}


double water_diffraction(struct rvec q, double en,
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


void get_diffraction(struct image *image, int na, int nb, int nc,
                     const double *intensities, const double *phases,
                     UnitCell *cell, int do_water, GradientMethod m)
{
	unsigned int xs, ys;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	float k, klow, bwstep;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* Allocate (and zero) the "diffraction array" */
	image->data = calloc(image->width * image->height, sizeof(float));

	/* Needed later for Lorentz calculation */
	image->twotheta = malloc(image->width * image->height * sizeof(double));

	k = 1.0/image->lambda;  /* Centre value */
	klow = k - k*(BANDWIDTH/2.0);  /* Lower value */
	bwstep = k * BANDWIDTH / BWSAMPLING;

	for ( xs=0; xs<image->width*SAMPLING; xs++ ) {
	for ( ys=0; ys<image->height*SAMPLING; ys++ ) {

		struct rvec q;
		float twotheta;
		double sw = 1.0/(SAMPLING*SAMPLING);  /* Sample weight */

		const unsigned int x = xs / SAMPLING;
		const unsigned int y = ys / SAMPLING; /* Integer part only */

		int kstep;

		for ( kstep=0; kstep<BWSAMPLING; kstep++ ) {

			float k;
			double kw = 1.0/BWSAMPLING;
			double intensity;
			double f_lattice, I_lattice;
			double I_molecule;

			/* Calculate k this time round */
			k = klow + kstep * bwstep;

			q = get_q(image, xs, ys, SAMPLING, &twotheta, k);
			image->twotheta[x + image->width*y] = twotheta;

			f_lattice = lattice_factor(q, ax, ay, az,
			                              bx, by, bz,
			                              cx, cy, cz,
			                              na, nb, nc);

			if ( intensities == NULL ) {
				I_molecule = 1.0e10;
			} else {
				I_molecule = molecule_factor(intensities,
				                             phases, q,
			                                     ax,ay,az,
			                                     bx,by,bz,cx,cy,cz,
				                             m);
			}

			I_lattice = pow(f_lattice, 2.0);

			intensity = sw * kw * I_lattice * I_molecule;
			image->data[x + image->width*y] += intensity;

		}

		if ( do_water ) {

			/* Bandwidth not simulated for water */
			struct rvec q;

			q = get_q(image, x, y, 1, NULL, 1.0/image->lambda);

			/* Add intensity contribution from water */
			image->data[x + image->width*y] += water_diffraction(q,
			                        ph_lambda_to_en(image->lambda),
			                        BEAM_RADIUS, WATER_RADIUS) * sw;

		}


	}
	progress_bar(xs, SAMPLING*image->width-1, "Calculating diffraction");
	}
}
