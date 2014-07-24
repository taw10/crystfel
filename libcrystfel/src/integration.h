/*
 * integration.h
 *
 * Integration of intensities
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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


#include "geometry.h"

typedef enum {

	INTDIAG_NONE,
	INTDIAG_RANDOM,
	INTDIAG_ALL,
	INTDIAG_INDICES,
	INTDIAG_NEGATIVE,
	INTDIAG_IMPLAUSIBLE,
	INTDIAG_STRONG

} IntDiag;

#define INTEGRATION_DEFAULTS_RINGS (INTEGRATION_RINGS)
#define INTEGRATION_DEFAULTS_PROF2D (INTEGRATION_PROF2D | INTEGRATION_CENTER)

/**
 * IntegrationMethod:
 * @INTEGRATION_NONE: No integration at all
 * @INTEGRATION_RINGS: Summation of pixel values inside ring, minus background
 * @INTEGRATION_PROF2D: Two dimensional profile fitting
 * @INTEGRATION_SATURATED: Integrate saturated reflections
 * @INTEGRATION_CENTER: Center the peak in the box prior to integration
 * @INTEGRATION_RESCUT: Stop integrating at the diffraction limit of the crystal
 *
 * An enumeration of all the available integration methods.
 **/
typedef enum {

	INTEGRATION_NONE   = 0,

	/* The core integration methods themselves */
	INTEGRATION_RINGS  = 1,
	INTEGRATION_PROF2D = 2,

	/* Bits at the top of the IntegrationMethod are flags which modify the
	 * behaviour of the integration. */
	INTEGRATION_SATURATED = 256,
	INTEGRATION_CENTER = 512,
	INTEGRATION_RESCUT = 1024,

} IntegrationMethod;

/* This defines the bits in "IntegrationMethod" which are used to represent the
 * core of the integration method */
#define INTEGRATION_METHOD_MASK (0xff)

#ifdef __cplusplus
extern "C" {
#endif

extern IntegrationMethod integration_method(const char *t, int *err);

extern void integrate_all(struct image *image, IntegrationMethod meth,
                          double ir_inn, double ir_mid, double ir_out,
                          IntDiag int_diag,
                          signed int idh, signed int idk, signed int idl);

extern void integrate_all_2(struct image *image, IntegrationMethod meth,
                            double push_res,
                            double ir_inn, double ir_mid, double ir_out,
                            IntDiag int_diag,
                            signed int idh, signed int idk, signed int idl);

extern void integrate_all_3(struct image *image, IntegrationMethod meth,
                            PartialityModel pmodel, double push_res,
                            double ir_inn, double ir_mid, double ir_out,
                            IntDiag int_diag,
                            signed int idh, signed int idk, signed int idl);

extern void integrate_all_4(struct image *image, IntegrationMethod meth,
                            PartialityModel pmodel, double push_res,
                            double ir_inn, double ir_mid, double ir_out,
                            IntDiag int_diag,
                            signed int idh, signed int idk, signed int idl,
                            int results_pipe);


#ifdef __cplusplus
}
#endif

#endif	/* INTEGRATION_H */
