/*
 * beam-parameters.h
 *
 * Beam parameters
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010,2012-2014 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani <valerio.mariani@desy.de>
 *   2012      Chunhong Yoon
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

#ifndef BEAM_PARAMETERS_H
#define BEAM_PARAMETERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct beam_params;
struct event;
struct hdfile;

#include "events.h"
#include "hdf5-file.h"

struct beam_params
{
	double fluence;        /* photons per pulse */
	double beam_radius;    /* metres */
	double photon_energy;  /* eV per photon */
	char *photon_energy_from; /* HDF5 dataset name */
	double photon_energy_scale;  /* Scale factor for photon energy, if the
	                              * energy is to be from the HDF5 file */
	double bandwidth;      /* FWHM(wavelength) over wavelength.
	                        *  Note: current simulation code just uses
	                        *        a rectangular distribution with this as
	                        *        its (full) width. */
	double divergence;     /* divergence (radians) */

	double profile_radius; /* Reciprocal space size of a reflection */
};

#ifdef __cplusplus
extern "C" {
#endif

extern struct beam_params *get_beam_parameters(const char *filename);
extern void free_beam_parameters(struct beam_params *beam);

extern void fill_in_beam_parameters(struct beam_params *beam, struct hdfile *f,
                             struct event* ev);

#ifdef __cplusplus
}
#endif

#endif	/* BEAM_PARAMETERS_H */
