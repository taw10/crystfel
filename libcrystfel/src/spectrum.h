/*
 * spectrum.h
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

#ifndef SPECTRUM_H
#define SPECTRUM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/**
 * \file spectrum.h
 * Data structure representing a radiation spectrum
 */

/**
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 **/
typedef struct _spectrum Spectrum;


/**
 * Structure representing a Gaussian distribution of spectral density.
 */
struct gaussian
{
	double kcen;    /**< k value at centre of Gaussian (in 1/m) */
	double sigma;   /**< Standard deviation of Gaussian (in 1/m) */
	double height;  /**< Height of Gaussian (arbitrary units) */
};


#ifdef __cplusplus
extern "C" {
#endif

/* Alloc/free */
extern Spectrum *spectrum_new(void);
extern void spectrum_free(Spectrum *s);
extern Spectrum *spectrum_load(const char *filename);

/* Representation as Gaussians */
extern void spectrum_set_gaussians(Spectrum *s, struct gaussian *gs,
                                   int n_gauss);
extern int spectrum_get_num_gaussians(Spectrum *s);
extern struct gaussian spectrum_get_gaussian(Spectrum *s, int n);

/* Representation as histogram */
extern void spectrum_set_histogram(Spectrum *s, double *kcens, double *heights,
                                   int nbins);
extern void spectrum_get_range(Spectrum *s, double *kmin, double *kmax);
extern double spectrum_get_density_at_k(Spectrum *s, double k);


#ifdef __cplusplus
}
#endif

#endif	/* SPECTRUM_H */
