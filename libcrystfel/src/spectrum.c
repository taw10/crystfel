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

#include <assert.h>
#include <gsl/gsl_sort.h>

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
	if ( s->rep == SPEC_HISTOGRAM ) {
		int i = 0;
		double frac;
		if ( k <= s->k[0] ) return 0.0;
		if ( k >= s->k[s->n_samples-1] ) return 0.0;
		/* k is definitely after the first sample, and definitely
		 * before the last one */
		while ( (s->k[i] < k) && (i<s->n_samples) ) i++;
		assert(i < s->n_samples);
		frac = (k - s->k[i-1]) / (s->k[i] - s->k[i-1]);
		return s->pdf[i-1] + frac * (s->pdf[i] - s->pdf[i-1]);
	}

	if ( s->rep == SPEC_GAUSSIANS ) {
		double total = 0.0;
		int i;
		for ( i=0; i<s->n_gaussians; i++ ) {
			double a = s->gaussians[i].height;
			double b = s->gaussians[i].kcen;
			double c = s->gaussians[i].sigma;
			total += a*exp(-(k-b)*(k-b)/(2.0*c*c));
		}
		return total;
	}

	return 0.0;
}


static double smallest_in_list(double *vals, int n_vals)
{
	int i;
	double v = +INFINITY;
	for ( i=0; i<n_vals; i++ ) {
		if ( vals[i] < v ) v = vals[i];
	}
	return v;
}


static double largest_in_list(double *vals, int n_vals)
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
		if ( gv < v ) v = gv;
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
 *
 * The values will be returned in units of 1/m.
 */
void spectrum_get_range(Spectrum *s, double *kmin, double *kmax)
{
	if ( s->rep == SPEC_HISTOGRAM ) {
		*kmin = smallest_in_list(s->k, s->n_samples);
		*kmax = largest_in_list(s->k, s->n_samples);
	} else {
		assert(s->rep == SPEC_GAUSSIANS);
		*kmin = gauss_low(s->gaussians, s->n_gaussians);
		*kmax = gauss_high(s->gaussians, s->n_gaussians);
	}
}


static signed int cmp_gauss(const void *va, const void *vb)
{
	const struct gaussian *a = va;
	const struct gaussian *b = vb;
	/* Integral of Gaussian = height * std deviation */
	if ( a->height * a->sigma > b->height * b->sigma ) return +1;
	return -1;
}


static void normalise_gaussians(struct gaussian *gauss, int n_gauss)
{
	int i;
	double total_area = 0.0;
	for ( i=0; i<n_gauss; i++ ) {
		total_area += gauss[i].height * gauss[i].sigma;
	}
	total_area *= sqrt(2.0 * M_PI);
	for ( i=0; i<n_gauss; i++ ) {
		gauss[i].height /= total_area;
	}
}


/**
 * \param s A \ref Spectrum
 * \param gs Pointer to array of \ref gaussian structures
 * \param n_gauss Number of Gaussians in \p gs
 *
 * Sets the spectrum in terms of a sum of Gaussians.
 *
 * The spectral density function will be normalised, so only the relative areas
 * under the curves are relevant.
 *
 * The input array will be copied, so you can safely free it after calling this
 * function.
 */
void spectrum_set_gaussians(Spectrum *s, struct gaussian *gs, int n_gauss)
{
	/* Free old contents (if any - may be NULL) */
	free(s->gaussians);
	free(s->k);
	free(s->pdf);

	s->gaussians = malloc(n_gauss * sizeof(struct gaussian));
	if ( s->gaussians == NULL ) return;

	memcpy(s->gaussians, gs, n_gauss*sizeof(struct gaussian));
	s->n_gaussians = n_gauss;
	s->rep = SPEC_GAUSSIANS;

	qsort(s->gaussians, s->n_gaussians, sizeof(struct gaussian), cmp_gauss);
	normalise_gaussians(s->gaussians, s->n_gaussians);
}


/* Samples must already have been sorted */
static void normalise_pdf(double *k, double *pdf, int n)
{
	int i;
	double total_area = 0.0;
	double old_k = k[0];
	double old_pdf = pdf[0];
	for ( i=1; i<n; i++ ) {
		total_area += (pdf[i]+old_pdf)*(k[i]-old_k)/2.0;
		old_k = k[i];
		old_pdf = pdf[i];
	}

	printf("total area under PDF = %f\n", total_area);

	for ( i=0; i<n; i++ ) {
		pdf[i] /= total_area;
	}
}


/**
 * \param s A \ref Spectrum
 * \param kvals Pointer to array of k values (in 1/m);
 * \param heights Pointer to array of spectral density samples
 * \param n Number of samples
 *
 * Sets the spectrum in terms of samples of the probability density function.
 * The spectral density function will be normalised.  The spectrum can have a
 * non-zero value only for k values in the range [kmin,kmax] (exclusive
 * interval), where kmin and kmax are the smallest and largest values in
 * \p kvals.
 *
 * The input arrays will be copied, so you can safely free them after calling
 * this function.
 */
void spectrum_set_pdf(Spectrum *s, double *kvals, double *heights, int n)
{
	/* Free old contents (if any - may be NULL) */
	free(s->gaussians);
	free(s->k);
	free(s->pdf);

	s->k = malloc(n * sizeof(double));
	if ( s->k == NULL ) return;

	s->pdf = malloc(n * sizeof(double));
	if ( s->pdf == NULL ) return;

	memcpy(s->k, kvals, n*sizeof(double));
	memcpy(s->pdf, heights, n*sizeof(double));
	s->n_samples = n;
	s->rep = SPEC_HISTOGRAM;

	gsl_sort2(s->k, 1, s->pdf, 1, n);
	normalise_pdf(s->k, s->pdf, s->n_samples);
}


static int read_esrf_spectrum(FILE *fh, Spectrum *s)
{
	double *k = NULL;
	double *samp = NULL;
	int n_bins = 0;
	int max_bins = 0;

	while ( !feof(fh) ) {

		float energy, weight;
		if ( fscanf(fh, "%e %e\n", &energy, &weight) != 2 ) return 1;

		if ( n_bins == max_bins ) {
			max_bins += 64;
			k = realloc(k, max_bins*sizeof(double));
			samp = realloc(samp, max_bins*sizeof(double));
			if ( (k==NULL) || (samp==NULL) ) return 1;
		}

		k[n_bins] = ph_eV_to_k(energy*1000.0);
		samp[n_bins] = weight;
		n_bins++;

	}

	spectrum_set_pdf(s, k, samp, n_bins);
	free(k);
	free(samp);

	return 0;
}


/**
 * \param filename Filename for the input file
 *
 * Loads the spectrum from \s filename.  Currently, only the ESRF spectrum file
 * format is supported.
 *
 * \returns A newly allocated \ref Spectrum, or NULL on error.
 */
Spectrum *spectrum_load(const char *filename)
{
	FILE *fh;
	Spectrum *s;
	char line[1024];

	fh = fopen(filename, "r");
	if ( fh == NULL ) return NULL;

	s = spectrum_new();
	if ( s == NULL ) {
		fclose(fh);
		return NULL;
	}

	if ( fgets(line, 1024, fh) != line ) {
		ERROR("Failed to read '%s'\n", filename);
		spectrum_free(s);
		return NULL;
	}

	chomp(line);
	if ( strcmp(line, "# energy/keV current/A") == 0 ) {
		if ( read_esrf_spectrum(fh, s) ) {
			ERROR("Failed to read ESRF spectrum from %s\n",
			      filename);
			spectrum_free(s);
			return NULL;
		}
	} else {
		ERROR("Spectrum format not recognised: %s\n", filename);
		fclose(fh);
		spectrum_free(s);
		return NULL;
	}

	fclose(fh);
	return s;
}
