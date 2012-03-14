/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <assert.h>
#include <fenv.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "diffraction.h"
#include "beam-parameters.h"
#include "symmetry.h"
#include "pattern_sim.h"


#define SINC_LUT_ELEMENTS (4096)


static double *get_sinc_lut(int n)
{
	int i;
	double *lut;

	lut = malloc(SINC_LUT_ELEMENTS*sizeof(double));
	lut[0] = n;
	if ( n == 1 ) {
		for ( i=1; i<SINC_LUT_ELEMENTS; i++ ) {
			lut[i] = 1.0;
		}
	} else {
		for ( i=1; i<SINC_LUT_ELEMENTS; i++ ) {
			double x, val;
			x = (double)i/SINC_LUT_ELEMENTS;
			val = fabs(sin(M_PI*n*x)/sin(M_PI*x));
			lut[i] = val;
		}
	}

	return lut;
}


static double interpolate_lut(double *lut, double val)
{
	double i, pos, f;
	unsigned int low, high;

	pos = SINC_LUT_ELEMENTS * modf(fabs(val), &i);
	low = (int)pos;  /* Discard fractional part */
	high = low + 1;
	f = modf(pos, &i);  /* Fraction */
	if ( high == SINC_LUT_ELEMENTS ) high = 0;

	return (1.0-f)*lut[low] + f*lut[high];
}


static double lattice_factor(struct rvec q, double ax, double ay, double az,
                                            double bx, double by, double bz,
                                            double cx, double cy, double cz,
                                            double *lut_a, double *lut_b,
                                            double *lut_c)
{
	struct rvec Udotq;
	double f1, f2, f3;

	Udotq.u = ax*q.u + ay*q.v + az*q.w;
	Udotq.v = bx*q.u + by*q.v + bz*q.w;
	Udotq.w = cx*q.u + cy*q.v + cz*q.w;

	f1 = interpolate_lut(lut_a, Udotq.u);
	f2 = interpolate_lut(lut_b, Udotq.v);
	f3 = interpolate_lut(lut_c, Udotq.w);

	return f1 * f2 * f3;
}


static double sym_lookup_intensity(const double *intensities,
                                   const unsigned char *flags,
                                   const SymOpList *sym,
                                   signed int h, signed int k, signed int l)
{
	int i;
	double ret = 0.0;

	for ( i=0; i<num_equivs(sym, NULL); i++ ) {

		signed int he;
		signed int ke;
		signed int le;
		double f, val;

		get_equiv(sym, NULL, i, h, k, l, &he, &ke, &le);

		f = (double)lookup_arr_flag(flags, he, ke, le);
		val = lookup_arr_intensity(intensities, he, ke, le);

		ret += f*val;

	}

	return ret;
}


static double sym_lookup_phase(const double *phases,
                               const unsigned char *flags, const SymOpList *sym,
                               signed int h, signed int k, signed int l)
{
	int i;
	double ret = 0.0;

	for ( i=0; i<num_equivs(sym, NULL); i++ ) {

		signed int he;
		signed int ke;
		signed int le;
		double f, val;

		get_equiv(sym, NULL, i, h, k, l, &he, &ke, &le);

		f = (double)lookup_arr_flag(flags, he, ke, le);
		val = lookup_arr_phase(phases, he, ke, le);

		ret += f*val;

	}

	return ret;
}


static double interpolate_linear(const double *ref, const unsigned char *flags,
                                 const SymOpList *sym, float hd,
                                 signed int k, signed int l)
{
	signed int h;
	double val1, val2;
	float f;

	h = (signed int)hd;
	if ( hd < 0.0 ) h -= 1;
	f = hd - (float)h;
	assert(f >= 0.0);

	val1 = sym_lookup_intensity(ref, flags, sym, h, k, l);
	val2 = sym_lookup_intensity(ref, flags, sym, h+1, k, l);

	val1 = val1;
	val2 = val2;

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_bilinear(const double *ref,
                                   const unsigned char *flags,
                                   const SymOpList *sym,
                                   float hd, float kd, signed int l)
{
	signed int k;
	double val1, val2;
	float f;

	k = (signed int)kd;
	if ( kd < 0.0 ) k -= 1;
	f = kd - (float)k;
	assert(f >= 0.0);

	val1 = interpolate_linear(ref, flags, sym, hd, k, l);
	val2 = interpolate_linear(ref, flags, sym, hd, k+1, l);

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_intensity(const double *ref,
                                    const unsigned char *flags,
                                    const SymOpList *sym,
                                    float hd, float kd, float ld)
{
	signed int l;
	double val1, val2;
	float f;

	l = (signed int)ld;
	if ( ld < 0.0 ) l -= 1;
	f = ld - (float)l;
	assert(f >= 0.0);

	val1 = interpolate_bilinear(ref, flags, sym, hd, kd, l);
	val2 = interpolate_bilinear(ref, flags, sym, hd, kd, l+1);

	return (1.0-f)*val1 + f*val2;
}


static double complex interpolate_phased_linear(const double *ref,
                                                const double *phases,
                                                const unsigned char *flags,
                                                const SymOpList *sym,
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

	val1 = sym_lookup_intensity(ref, flags, sym, h, k, l);
	val2 = sym_lookup_intensity(ref, flags, sym, h+1, k, l);
	ph1 = sym_lookup_phase(phases, flags, sym, h, k, l);
	ph2 = sym_lookup_phase(phases, flags, sym, h+1, k, l);

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
                                                  const unsigned char *flags,
                                                  const SymOpList *sym,
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

	val1 = interpolate_phased_linear(ref, phases, flags, sym, hd, k, l);
	val2 = interpolate_phased_linear(ref, phases, flags, sym, hd, k+1, l);

	return (1.0-f)*val1 + f*val2;
}


static double interpolate_phased_intensity(const double *ref,
                                           const double *phases,
                                           const unsigned char *flags,
                                           const SymOpList *sym,
                                           float hd, float kd, float ld)
{
	signed int l;
	double complex val1, val2;
	float f;

	l = (signed int)ld;
	if ( ld < 0.0 ) l -= 1;
	f = ld - (float)l;
	assert(f >= 0.0);

	val1 = interpolate_phased_bilinear(ref, phases, flags, sym,
	                                   hd, kd, l);
	val2 = interpolate_phased_bilinear(ref, phases, flags, sym,
	                                   hd, kd, l+1);

	return cabs((1.0-f)*val1 + f*val2);
}


/* Look up the structure factor for the nearest Bragg condition */
static double molecule_factor(const double *intensities, const double *phases,
                              const unsigned char *flags, struct rvec q,
                              double ax, double ay, double az,
                              double bx, double by, double bz,
                              double cx, double cy, double cz,
                              GradientMethod m, const SymOpList *sym)
{
	float hd, kd, ld;
	signed int h, k, l;
	double r;

	hd = q.u * ax + q.v * ay + q.w * az;
	kd = q.u * bx + q.v * by + q.w * bz;
	ld = q.u * cx + q.v * cy + q.w * cz;

	/* No flags -> flat intensity distribution */
	if ( flags == NULL ) return 1.0e5;

	switch ( m ) {

		case GRADIENT_MOSAIC :
		fesetround(1);  /* Round to nearest */
		h = (signed int)rint(hd);
		k = (signed int)rint(kd);
		l = (signed int)rint(ld);
		if ( abs(h) > INDMAX ) r = 0.0;
		else if ( abs(k) > INDMAX ) r = 0.0;
		else if ( abs(l) > INDMAX ) r = 0.0;
		else r = sym_lookup_intensity(intensities, flags, sym, h, k, l);
		break;

		case GRADIENT_INTERPOLATE :
		r = interpolate_intensity(intensities, flags, sym, hd, kd, ld);
		break;

		case GRADIENT_PHASED :
		r = interpolate_phased_intensity(intensities, phases, flags,
		                                 sym, hd, kd, ld);
		break;

		default:
		ERROR("This gradient method not implemented yet.\n");
		exit(1);
	}

	return r;
}


void get_diffraction(struct image *image, int na, int nb, int nc,
                     const double *intensities, const double *phases,
                     const unsigned char *flags, UnitCell *cell,
                     GradientMethod m, const SymOpList *sym)
{
	unsigned int fs, ss;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double *lut_a;
	double *lut_b;
	double *lut_c;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* Allocate (and zero) the "diffraction array" */
	image->data = calloc(image->width * image->height, sizeof(float));

	/* Needed later for Lorentz calculation */
	image->twotheta = malloc(image->width * image->height * sizeof(double));

	lut_a = get_sinc_lut(na);
	lut_b = get_sinc_lut(nb);
	lut_c = get_sinc_lut(nc);

	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		int idx;
		double k;
		double f_lattice, I_lattice;
		double I_molecule;
		struct rvec q;
		double twotheta;

		/* Calculate k this time round */
		k = 1.0/image->lambda;

		q = get_q(image, fs, ss, &twotheta, k);

		f_lattice = lattice_factor(q, ax, ay, az,
		                              bx, by, bz,
		                              cx, cy, cz,
		                              lut_a, lut_b, lut_c);

		I_molecule = molecule_factor(intensities,
		                             phases, flags, q,
		                             ax,ay,az,bx,by,bz,cx,cy,cz,
		                             m, sym);

		I_lattice = pow(f_lattice, 2.0);

		idx = fs + image->width*ss;
		image->data[idx] = I_lattice * I_molecule;
		image->twotheta[idx] = twotheta;


	}
	progress_bar(fs, image->width-1, "Calculating diffraction");
	}

	free(lut_a);
	free(lut_b);
	free(lut_c);
}
