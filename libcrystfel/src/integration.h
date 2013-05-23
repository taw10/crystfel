/*
 * integration.h
 *
 * Integration of intensities
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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

#ifndef INTEGRATION_H
#define INTEGRATION_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define INTEGRATION_DEFAULTS_RINGS (INTEGRATION_RINGS | INTEGRATION_SATURATED)
#define INTEGRATION_DEFAULTS_REFINE (INTEGRATION_REFINE | INTEGRATION_SATURATED)

/**
 * IntegrationMethod:
 * @INTEGRATION_NONE: No integration at all
 * @INTEGRATION_RINGS: Predict reflections and integrate them all
 * @INTEGRATION_REFINE: As @INTEGRATION_RINGS, but
 *
 * @INTEGRATION_SATURATED: Integrate saturated reflections
 *
 * An enumeration of all the available integration methods.
 **/
typedef enum {

	INTEGRATION_NONE   = 0,

	/* The core integration methods themselves */
	INTEGRATION_RINGS  = 1,
	INTEGRATION_REFINE = 2,

	/* Bits at the top of the IntegrationMethod are flags which modify the
	 * behaviour of the integration. */
	INTEGRATION_SATURATED = 256,
	INTEGRATION_CLOSER = 512,

} IntegrationMethod;

/* This defines the bits in "IntegrationMethod" which are used to represent the
 * core of the integration method */
#define INTEGRATION_METHOD_MASK (0xff)

extern IntegrationMethod integration_method(const char *t, int *err);

extern void integrate_all(struct image *image, IntegrationMethod meth,
	                  int use_closer, double min_snr,
	                  double ir_inn, double ir_mid, double ir_out,
	                  int integrate_saturated);

#endif	/* INTEGRATION_H */
