/*
 * beam-parameters.h
 *
 * Beam parameters
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010,2012 Thomas White <taw@physics.org>
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


struct beam_params
{
	double fluence;        /* photons per pulse */
	double beam_radius;    /* metres */
	double photon_energy;  /* eV per photon */
	double bandwidth;      /* FWHM(wavelength) over wavelength.
	                        *  Note: current simulation code just uses
	                        *        a rectangular distribution with this as
	                        *        its (full) width. */
	double divergence;     /* divergence (radians) */
};


extern struct beam_params *get_beam_parameters(const char *filename);


#endif	/* BEAM_PARAMETERS_H */
