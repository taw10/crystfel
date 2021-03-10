/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2020 Thomas White <taw@physics.org>
 *   2013-2014 Chun Hong Yoon <chun.hong.yoon@desy.de>
 *   2013      Alexandra Tolstikova
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
#include "symmetry.h"
#include "pattern_sim.h"


#define SINC_LUT_ELEMENTS (4096)


static double *get_sinc_lut(int n, int no_fringes, int flat)
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
			if ( (flat || no_fringes) && (x > 1.0/n) && (1.0-x > 1.0/n) ) {
				val = 0.0;
			} else if ( flat ) {
				val = n;
			} else {
				val = fabs(sin(M_PI*n*x)/sin(M_PI*x));
			}
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

	for ( i=0; i<num_equivs(sym, NULL); i++ ) {

		signed int he;
		signed int ke;
		signed int le;
		int f;

		get_equiv(sym, NULL, i, h, k, l, &he, &ke, &le);

		f = lookup_arr_flag(flags, he, ke, le);

		if ( f ) return lookup_arr_phase(phases, he, ke, le);

	}

	return 0.0;
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
	if ( flags == NULL ) return 100.0;

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


static void diffraction_panel(struct image *image, const double *intensities,
                              const double *phases, const unsigned char *flags,
                              UnitCell *cell, GradientMethod m,
                              const SymOpList *sym, double k,
                              double ax, double ay, double az,
                              double bx, double by, double bz,
                              double cx, double cy, double cz,
                              double *lut_a, double *lut_b, double *lut_c,
                              int pn, double weight)
{
	int fs, ss;
	const int nxs = 4;
	const int nys = 4;
	struct detgeom_panel *p = &image->detgeom->panels[pn];

	weight /= nxs*nys;

	for ( ss=0; ss<p->h; ss++ ) {
	for ( fs=0; fs<p->w; fs++ ) {

		int idx;
		double f_lattice, I_lattice;
		double I_molecule;
		int xs, ys;
		float xo, yo;

		for ( xs=0; xs<nxs; xs++ ) {
		for ( ys=0; ys<nys; ys++ ) {

			double qv[3];
			struct rvec q;

			xo = (1.0/nxs) * xs;
			yo = (1.0/nys) * ys;

			detgeom_transform_coords(p, fs+xo, ss+yo,
			                         0.0, 0.0, 1.0/k, qv);

			q.u = qv[0]; q.v = qv[1]; q.w = qv[2];

			f_lattice = lattice_factor(q, ax, ay, az,
					           bx, by, bz,
					           cx, cy, cz,
					           lut_a, lut_b, lut_c);

			I_molecule = molecule_factor(intensities,
					             phases, flags, q,
					             ax, ay, az,
					             bx, by, bz,
					             cx, cy, cz,
					             m, sym);

			I_lattice = pow(f_lattice, 2.0);

			idx = fs + p->w*ss;
			image->dp[pn][idx] += I_lattice * I_molecule * weight;

		}
		}
	}
	progress_bar(ss, p->h-1, "Calculating diffraction");
	}
}


static void diffraction_at_k(struct image *image, const double *intensities,
                             const double *phases, const unsigned char *flags,
                             UnitCell *cell, GradientMethod m,
                             const SymOpList *sym, double k,
                             double ax, double ay, double az,
                             double bx, double by, double bz,
                             double cx, double cy, double cz,
                             double *lut_a, double *lut_b, double *lut_c,
                             double weight)
{
	int i;

	for ( i=0; i<image->detgeom->n_panels; i++ ) {
		diffraction_panel(image, intensities, phases, flags, cell, m,
		                  sym, k, ax, ay, az, bx, by, bz, cx, cy, cz,
		                  lut_a, lut_b, lut_c, i, weight);
	}
}


void get_diffraction(struct image *image, int na, int nb, int nc,
                     const double *intensities, const double *phases,
                     const unsigned char *flags, UnitCell *cell,
                     GradientMethod m, const SymOpList *sym,
                     int no_fringes, int flat, int n_samples)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double *lut_a;
	double *lut_b;
	double *lut_c;
	int i;
	double kmin, kmax, step;
	double norm = 0.0;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	lut_a = get_sinc_lut(na, no_fringes, flat);
	lut_b = get_sinc_lut(nb, no_fringes, flat);
	lut_c = get_sinc_lut(nc, no_fringes, flat);

	spectrum_get_range(image->spectrum, &kmin, &kmax);
	step = (kmax-kmin)/(n_samples+1);

	/* Determine normalisation factor such that weights add up to 1 after
	 * sampling (bins must have constant width) */
	for ( i=1; i<=n_samples; i++ ) {
		double k = kmin + i*step;
		norm += spectrum_get_density_at_k(image->spectrum, k);
	}
	for ( i=1; i<=n_samples; i++ ) {

		double k = kmin + i*step;
		double prob;

		/* Probability = p.d.f. times step width */
		prob = spectrum_get_density_at_k(image->spectrum, k)/norm;

		STATUS("Wavelength: %e m, weight = %.5f\n", 1.0/k, prob);

		diffraction_at_k(image, intensities, phases,
		                flags, cell, m, sym, k,
		                ax, ay, az, bx, by, bz, cx, cy, cz,
		                lut_a, lut_b, lut_c, prob);

	}

	free(lut_a);
	free(lut_b);
	free(lut_c);
}
