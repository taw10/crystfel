/*
 * spectrum.c
 *
 * A class representing a radiation spectrum
 *
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "spectrum.h"
#include "utils.h"


/**
 * \file spectrum.h
 */


enum spectrumrep
{
	SPEC_HISTOGRAM,
	SPEC_GAUSSIANS
};

struct _spectrum
{
	enum spectrumrep rep;

	/* Gaussian representation */
	struct gaussian *gaussians;
	int n_gaussians;

	/* Histogram representation */
	double *k;
	double *pdf;
	int n_samples;
};


/**
 * Create a new \ref Spectrum.
 *
 * \returns The new spectrum, or NULL on failure.
 *
 */
Spectrum *spectrum_new()
{
	Spectrum *s;

	s = malloc(sizeof(Spectrum));
	if ( s == NULL ) return NULL;

	s->rep = SPEC_GAUSSIANS;

	s->gaussians = NULL;
	s->n_gaussians = 0;

	s->k = NULL;
	s->pdf = NULL;
	s->n_samples = 0;

	return s;
}


/**
 * \param s A \ref Spectrum
 *
 * Frees a \ref Spectrum.
 */
void spectrum_free(Spectrum *s)
{
	if ( s == NULL ) return;
	free(s->gaussians);
	free(s->k);
	free(s->pdf);
	free(s);
}


/**
 * \param s A \ref Spectrum
 *
 * \returns The number of Gaussians in the spectrum, or zero if \p s is not
 * currently represented as Gaussians.
 */
int spectrum_get_num_gaussians(Spectrum *s)
{
	if ( s->rep == SPEC_GAUSSIANS ) return s->n_gaussians;
	return 0;
}


/**
 * \param s A \ref Spectrum
 * \param n The index number of the required Gaussian
 *
 * Returns The \p n-th Gaussian in the spectrum.  The Gaussians are
 * returned in descending order of integrated intensity, indexed from zero.
 *
 * If \p n is greater than or equal to the number of Gaussians in the spectrum,
 * or if the spectrum is not represented as Gaussians, the returned Gaussian
 * will have zero height.
 *
 * \returns The \p n-th Gaussian.
 */
struct gaussian spectrum_get_gaussian(Spectrum *s, int n)
{
	struct gaussian g;
	if ( (s->rep != SPEC_GAUSSIANS) || (n >= s->n_gaussians) ) {
		g.kcen = 0.0;
		g.sigma = 0.0;
		g.height = 0.0;
	} else {
		g = s->gaussians[n];
	}
	return g;
}


/**
 * \param s A \ref Spectrum
 * \param k A wavenumber (in 1/metres)
 *
 * Retrieves the spectral density at wavenumber \p k.
 * This is a sample from a probability density function, so to calculate the
 * "amount of intensity" from this, you'll need to multiply the value by a
 * small width of k.
 *
 * \returns The density at \p k.
 */
double spectrum_get_density_at_k(Spectrum *s, double k)
{
	return 0.0;
}


static double smallest(double *vals, int n_vals)
{
	int i;
	double v = +INFINITY;
	for ( i=0; i<n_vals; i++ ) {
		if ( vals[i] < v ) v = vals[i];
	}
	return v;
}


static double largest(double *vals, int n_vals)
{
	int i;
	double v = -INFINITY;
	for ( i=0; i<n_vals; i++ ) {
		if ( vals[i] > v ) v = vals[i];
	}
	return v;
}


static double gauss_low(struct gaussian *gauss, int n_gauss)
{
	int i;
	double v = +INFINITY;
	for ( i=0; i<n_gauss; i++ ) {
		double gv = gauss[i].kcen - 5.0*gauss[i].sigma;
		if ( gv > v ) v = gv;
	}
	return v;
}


static double gauss_high(struct gaussian *gauss, int n_gauss)
{
	int i;
	double v = -INFINITY;
	for ( i=0; i<n_gauss; i++ ) {
		double gv = gauss[i].kcen + 5.0*gauss[i].sigma;
		if ( gv > v ) v = gv;
	}
	return v;
}


/**
 * \param s A \ref Spectrum
 * \param kmin Location to store minimum k value
 * \param kmax Location to store maximum k value
 *
 * Sets \p kmin and \p kmax to the range of k values in the spectrum.  If the
 * spectrum is represented as Gaussians, the range returned will be enough to
 * contain at least 5 sigmas of all the Gaussians.
 */
void spectrum_get_range(Spectrum *s, double *kmin, double *kmax)
{
	if ( s->rep == SPEC_HISTOGRAM ) {
		*kmin = smallest(s->k, s->n_samples);
		*kmax = largest(s->k, s->n_samples);
	} else {
		assert(s->rep == SPEC_GAUSSIANS);
		*kmin = gauss_low(s->gaussians, s->n_gaussians);
		*kmax = gauss_high(s->gaussians, s->n_gaussians);
	}
}


/**
 * \param s A \ref Spectrum
 * \param gs Pointer to array of \ref gaussian structures
 * \param n_gauss Number of Gaussians in \p gs
 *
 * Sets the spectrum in terms of a sum of Gaussians.
 */
void spectrum_set_gaussians(Spectrum *s, struct gaussian *gs, int n_gauss)
{
	/* FIXME: sort them */
}


/**
 * \param s A \ref Spectrum
 * \param kcens Pointer to array of k values in the centres of the bins
 * \param heights Pointer to array of spectral density values
 * \param nbins Number of bins
 *
 * Sets the spectrum in terms of a histogram.
 */
void spectrum_set_histogram(Spectrum *s, double *kcens, double *heights,
                            int nbins)
{
}


/**
 * \param filename Filename for the input file
 *
 * Loads the spectrum from \s filename.
 *
 * \returns A newly allocated \ref Spectrum, or NULL on error.
 */
Spectrum *spectrum_load(const char *filename)
{
	return NULL;
}
