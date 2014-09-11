/*
 * diffraction.c
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2014 Thomas White <taw@physics.org>
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
#include "beam-parameters.h"
#include "symmetry.h"
#include "pattern_sim.h"


#define SINC_LUT_ELEMENTS (4096)


static double *get_sinc_lut(int n, int no_fringes)
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
			if ( no_fringes && (x > 1.0/n) && (1.0-x > 1.0/n) ) {
				val = 0.0;
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
	unsigned int fs, ss;
	const int nxs = 4;
	const int nys = 4;

	weight /= nxs*nys;

	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		int idx;
		double f_lattice, I_lattice;
		double I_molecule;
		struct rvec q;
		double twotheta;
		int xs, ys;
		float xo, yo;

		for ( xs=0; xs<nxs; xs++ ) {
		for ( ys=0; ys<nys; ys++ ) {

			xo = (1.0/nxs) * xs;
			yo = (1.0/nys) * ys;

			q = get_q(image, fs+xo, ss+yo, &twotheta, k);

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

			idx = fs + image->width*ss;
			image->data[idx] += I_lattice * I_molecule * weight;
			image->twotheta[idx] = twotheta;

		}
		}
	}
	progress_bar(fs, image->width-1, "Calculating diffraction");
	}
}


static int compare_samples(const void *a, const void *b)
{
	struct sample *sample1 = (struct sample *)a;
	struct sample *sample2 = (struct sample *)b;
	if ( sample1->weight < sample2->weight ) {
		return 1;
	}
	return -1;
}


static struct sample *get_gaussian_spectrum(double eV_cen, double eV_step,
                                            double sigma, int spec_size,
                                            double eV_start)
{
	struct sample *spectrum;
	int i;
	double eV;

	spectrum = malloc(spec_size * sizeof(struct sample));
	if ( spectrum == NULL ) return NULL;

	if (eV_start == 0) { /* eV starts at 3 sigma below the mean*/
		eV = eV_cen - (spec_size/2)*eV_step;
	} else {
		eV = eV_start;
	}

	for ( i=0; i<spec_size; i++ ) {

		spectrum[i].k = 1.0/ph_eV_to_lambda(eV);
		spectrum[i].weight = exp(-(pow(eV - eV_cen, 2.0)
		                         / (2.0*sigma*sigma)));
		eV += eV_step;
	}

	return spectrum;
}


static int add_sase_noise(struct sample *spectrum, int nsteps, gsl_rng *rng)
{
	struct sample *noise;
	int i, j;
	double *gaussianNoise;
	int shiftLim = 5;
	double noise_mean = 0.0;
	double noise_sigma = 1.0;

	if ( shiftLim > nsteps ) shiftLim = nsteps;

	noise = malloc(nsteps * sizeof(struct sample));
	if ( noise == NULL ) return 1;

	gaussianNoise = malloc(3 * nsteps * sizeof(double));
	if ( gaussianNoise == NULL ) {
		free(noise);
		return 1;
	}

	/* Generate Gaussian noise of length of spectrum
	 * (replicate on both ends for circshift below) */
	for ( i=0; i<nsteps; i++) {

		noise[i].weight = 0.0;

		/* Gaussian noise with mean = 0, std = 1 */
		gaussianNoise[i] = gaussian_noise(rng, noise_mean, noise_sigma);
		gaussianNoise[i+nsteps] = gaussianNoise[i];
		gaussianNoise[i+2*nsteps] = gaussianNoise[i];
	}

	/* Sum Gaussian noise by circshifting by +/- shiftLim */
	for ( i=nsteps; i<2*nsteps; i++ ) {
		for ( j=-shiftLim; j<=shiftLim; j++ ) {
			noise[i-nsteps].weight += gaussianNoise[i+j];
		}
	}

	/* Normalize the number of circshift sum */
	for ( i=0; i<nsteps; i++) {
		noise[i].weight = noise[i].weight/(2*shiftLim+1);
	}

	/* Noise amplitude should have a 2 x Gaussian distribution */
	for ( i=0; i<nsteps; i++ ) {
		noise[i].weight = 2.0 * spectrum[i].weight * noise[i].weight;
	}

	/* Add noise to spectrum */
	for ( i=0; i<nsteps; i++ ) {

		spectrum[i].weight += noise[i].weight;

		/* The final spectrum can not be negative */
		if ( spectrum[i].weight < 0.0 ) spectrum[i].weight = 0.0;

	}

	return 0;
}


struct sample *generate_tophat(struct image *image)
{
	struct sample *spectrum;
	int i;
	double k, k_step;

	double halfwidth = image->bw * image->lambda / 2.0;  /* m */
	double mink = 1.0/(image->lambda + halfwidth);
	double maxk = 1.0/(image->lambda - halfwidth);

	spectrum = malloc(image->nsamples * sizeof(struct sample));
	if ( spectrum == NULL ) return NULL;

	k = mink;
	k_step = (maxk-mink)/(image->nsamples-1);
	for ( i=0; i<image->nsamples; i++ ) {
		spectrum[i].k = k;
		spectrum[i].weight = 1.0/(double)image->nsamples;
		k += k_step;
	}

	image->spectrum_size = image->nsamples;

	return spectrum;
}


struct sample *generate_SASE(struct image *image, gsl_rng *rng)
{
	struct sample *spectrum;
	int i;
	const int spec_size = 1024;
	double eV_cen;  /* Central photon energy for this spectrum */
	const double jitter_sigma_eV = 8.0;

	/* Central wavelength jitters with Gaussian distribution */
	eV_cen = gaussian_noise(rng, ph_lambda_to_eV(image->lambda),
	                        jitter_sigma_eV);

	/* Convert FWHM to standard deviation.  Note that bandwidth is taken to
	 * be "delta E over E" (E = photon energy), not the bandwidth in terms
	 * of wavelength, but the difference should be very small */
	double sigma = (image->bw*eV_cen) / (2.0*sqrt(2.0*log(2.0)));

	/* The spectrum will be calculated to a resolution which spreads six
	 * sigmas of the original (no SASE noise) Gaussian pulse over spec_size
	 * points */
	double eV_step = 6.0*sigma/(spec_size-1);

	spectrum = get_gaussian_spectrum(eV_cen, eV_step, sigma, spec_size,0);

	/* Add SASE-type noise to Gaussian spectrum */
	add_sase_noise(spectrum, spec_size, rng);

	/* Normalise intensity (before taking restricted number of samples) */
	double total_weight = 0.0;
	for ( i=0; i<spec_size; i++ ) {
		total_weight += spectrum[i].weight;
	}
	for ( i=0; i<spec_size; i++ ) {
		spectrum[i].weight /= total_weight;
	}

	/* Sort samples in spectrum by weight.  Diffraction calculation will
	 * take the requested number, starting from the brightest */
	qsort(spectrum, spec_size, sizeof(struct sample), compare_samples);

	image->spectrum_size = spec_size;

	return spectrum;
}


struct sample *generate_twocolour(struct image *image)
{
	struct sample *spectrum;
	struct sample *spectrum1;
	struct sample *spectrum2;
	int i;
	double eV_cen1;  /* Central photon energy for first colour */
	double eV_cen2;  /* Central photon energy for second colour */
	double eV_cen;  /* Central photon energy for this spectrum */
	const int spec_size = 1024;

	eV_cen = ph_lambda_to_eV(image->lambda);

	double halfwidth = eV_cen*image->bw/2.0; /* eV */

	eV_cen1 = eV_cen - halfwidth;
	eV_cen2 = eV_cen + halfwidth;

	/* Hard-code sigma to be 1/5 of bandwidth */
	double sigma = eV_cen*image->bw/5.0; /* eV */

	/* The spectrum will be calculated to a resolution which spreads six
	* sigmas of the original (no SASE noise) Gaussian pulse over spec_size
	* points */
	double eV_start = eV_cen1 - 3*sigma;
	double eV_end = eV_cen2 + 3*sigma;
	double eV_step = (eV_end - eV_start)/(spec_size-1);

	spectrum1 = get_gaussian_spectrum(eV_cen1, eV_step, sigma, spec_size,
	                                  eV_start);
	spectrum2 = get_gaussian_spectrum(eV_cen2, eV_step, sigma, spec_size,
	                                  eV_start);

	spectrum = malloc(spec_size * sizeof(struct sample));
	if ( spectrum == NULL ) return NULL;

	for ( i=0; i<spec_size; i++ ) {
		spectrum[i].weight = spectrum1[i].weight + spectrum2[i].weight;
		spectrum[i].k = spectrum1[i].k;
		if ( spectrum2[i].k != spectrum1[i].k ) {
			printf("%e %e\n", spectrum1[i].k, spectrum2[i].k);
		}
	}

	/* Normalise intensity (before taking restricted number of samples) */
	double total_weight = 0.0;
	for ( i=0; i<spec_size; i++ ) {
		total_weight += spectrum[i].weight;
	}

	for ( i=0; i<spec_size; i++ ) {
		spectrum[i].weight /= total_weight;
	}

	/* Sort samples in spectrum by weight.  Diffraction calculation will
	* take the requested number, starting from the brightest */
	qsort(spectrum, spec_size, sizeof(struct sample), compare_samples);

	image->spectrum_size = spec_size;

	return spectrum;
}


void get_diffraction(struct image *image, int na, int nb, int nc,
                     const double *intensities, const double *phases,
                     const unsigned char *flags, UnitCell *cell,
                     GradientMethod m, const SymOpList *sym, int no_fringes)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double *lut_a;
	double *lut_b;
	double *lut_c;
	int i;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* Allocate (and zero) the "diffraction array" */
	image->data = calloc(image->width * image->height, sizeof(float));

	/* Needed later for Lorentz calculation */
	image->twotheta = malloc(image->width * image->height * sizeof(double));

	lut_a = get_sinc_lut(na, no_fringes);
	lut_b = get_sinc_lut(nb, no_fringes);
	lut_c = get_sinc_lut(nc, no_fringes);

	for ( i=0; i<image->nsamples; i++ ) {

		printf("%.1f eV, weight = %.5f\n",
		       ph_lambda_to_eV(1.0/image->spectrum[i].k),
		       image->spectrum[i].weight);

		diffraction_at_k(image, intensities, phases,
		                flags, cell, m, sym, image->spectrum[i].k,
		                ax, ay, az, bx, by, bz, cx, cy, cz,
		                lut_a, lut_b, lut_c, image->spectrum[i].weight);


	}

	free(lut_a);
	free(lut_b);
	free(lut_c);
}
