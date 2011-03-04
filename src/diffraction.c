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
#include <fenv.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "diffraction.h"
#include "sfac.h"
#include "beam-parameters.h"
#include "symmetry.h"


#define SAMPLING (4)
#define BWSAMPLING (10)


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


static double sym_lookup_intensity(const double *intensities,
                                   const unsigned char *flags, const char *sym,
                                   signed int h, signed int k, signed int l)
{
	int i;
	double ret = 0.0;

	for ( i=0; i<num_general_equivs(sym); i++ ) {

		signed int he;
		signed int ke;
		signed int le;
		double f, val;

		get_general_equiv(h, k, l, &he, &ke, &le, sym, i);

		f = (double)lookup_flag(flags, he, ke, le);
		val = lookup_intensity(intensities, he, ke, le);

		ret += f*val;

	}

	return ret;
}


static double sym_lookup_phase(const double *phases,
                               const unsigned char *flags, const char *sym,
                               signed int h, signed int k, signed int l)
{
	int i;
	double ret = 0.0;

	for ( i=0; i<num_general_equivs(sym); i++ ) {

		signed int he;
		signed int ke;
		signed int le;
		double f, val;

		get_general_equiv(h, k, l, &he, &ke, &le, sym, i);

		f = (double)lookup_flag(flags, he, ke, le);
		val = lookup_phase(phases, he, ke, le);

		ret += f*val;

	}

	return ret;
}


static double interpolate_linear(const double *ref, const unsigned char *flags,
                                 const char *sym, float hd,
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
                                   const unsigned char *flags, const char *sym,
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
                                    const unsigned char *flags, const char *sym,
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
                                                const char *sym,
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
                                                  const char *sym,
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
                                           const char *sym,
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
                              GradientMethod m, const char *sym)
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
                     GradientMethod m, const char *sym)
{
	unsigned int fs, ss;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	float klow, khigh, bwstep;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* Allocate (and zero) the "diffraction array" */
	image->data = calloc(image->width * image->height, sizeof(float));

	/* Needed later for Lorentz calculation */
	image->twotheta = malloc(image->width * image->height * sizeof(double));

	klow = 1.0/(image->lambda*(1.0 + image->beam->bandwidth/2.0));
	khigh = 1.0/(image->lambda*(1.0 - image->beam->bandwidth/2.0));
	bwstep = (khigh-klow) / BWSAMPLING;

	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		int fs_step, ss_step, kstep;
		int idx = fs + image->width*ss;

		for ( fs_step=0; fs_step<SAMPLING; fs_step++ ) {
		for ( ss_step=0; ss_step<SAMPLING; ss_step++ ) {
		for ( kstep=0; kstep<BWSAMPLING; kstep++ ) {

			double k;
			double intensity;
			double f_lattice, I_lattice;
			double I_molecule;
			struct rvec q;
			double twotheta;
			const double dfs = fs + (fs_step / SAMPLING);
			const double dss = ss + (ss_step / SAMPLING);

			/* Calculate k this time round */
			k = klow + kstep * bwstep;

			q = get_q(image, dfs, dss, &twotheta, k);

			f_lattice = lattice_factor(q, ax, ay, az,
			                              bx, by, bz,
			                              cx, cy, cz,
			                              na, nb, nc);

			I_molecule = molecule_factor(intensities,
			                             phases, flags, q,
			                             ax,ay,az,bx,by,bz,cx,cy,cz,
			                             m, sym);

			I_lattice = pow(f_lattice, 2.0);
			intensity = I_lattice * I_molecule;

			image->data[idx] += intensity;

			if ( fs_step + ss_step + kstep == 0 ) {
				image->twotheta[idx] = twotheta;
			}

		}
		}
		}

		image->data[idx] /= SAMPLING*SAMPLING*BWSAMPLING;


	}
	progress_bar(fs, image->width-1, "Calculating diffraction");
	}
}
