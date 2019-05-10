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

	return s;
}


/**
 * \param s A \ref Spectrum
 *
 * Frees a \ref Spectrum.
 */
void spectrum_free(Spectrum *s)
{
	free(s);
}


/**
 * \param s A \ref Spectrum
 * \param tol Fraction of spectrum required
 *
 * \returns The number of Gaussians needed to represent at least \p tol
 * fraction of the total intensity in the spectrum.
 */
int spectrum_get_num_gaussians(Spectrum *s, double tol)
{
	return 0;
}


/**
 * \param s A \ref Spectrum
 * \param n The index number of the required Gaussian
 *
 * Returns The \p n-th Gaussian to represent the spectrum.  The Gaussians are
 * returned in descending order of integrated intensity, indexed from zero.
 *
 * If \p n is greater than or equal to the number of Gaussians needed to account
 * for 100% of the spectrum, the returned Gaussian will have zero height.
 *
 * \returns The \p n-th Gaussian.
 */
struct gaussian spectrum_get_gaussian(Spectrum *s, int n)
{
	struct gaussian g;
	g.kcen = 0.0;
	g.sigma = 0.0;
	g.height = 0.0;
	return g;
}


/**
 * \param s A \ref Spectrum
 * \param k A wavenumber (in 1/metres)
 *
 * Returns the spectral density at wavenumber \p k.
 * To calculate the "amount of intensity" from this, you'll need to multiply
 * the value by a small width of k.
 *
 * \returns The density at \p k.
 */
double spectrum_get_density_at_k(Spectrum *s, double k)
{
	return 0.0;
}


/**
 * \param s A \ref Spectrum
 * \param tol Fraction of spectrum required
 * \param kmin Location to store minimum k value
 * \param kmax Location to store maximum k value
 *
 * Returns the range of k values needed to capture \p tol of the spectrum.
 */
void spectrum_get_range(Spectrum *s, double tol, double *kmin, double *kmax)
{
}


/**
 * \param s A \ref Spectrum
 * \param gs Pointer to array of \ref gaussian structures
 * \param n_gauss Number of Gaussians in \p gs
 *
 * Sets the spectrum in terms of a number of Gaussians.
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
