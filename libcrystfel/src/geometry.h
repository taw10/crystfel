/*
 * geometry.h
 *
 * Geometry of diffraction
 *
 * Copyright © 2013-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
 *   2012      Richard Kirian
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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "reflist.h"
#include "cell.h"
#include "crystal.h"
#include "detgeom.h"
#include "image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file geometry.h
 * Geometry of diffraction
 *
 * This contains the prediction and partiality calculation functions.
 */

/**
 * A PartialityModel describes a geometrical model which can be used to
 * calculate spot partialities and Lorentz correction factors.
 **/
typedef enum {

	PMODEL_UNITY,   /**< Set all partialities and Lorentz factors to 1. */
	PMODEL_XSPHERE, /**< Flat sphere model with super-Gaussian spectrum */
	PMODEL_OFFSET,  /**< Ewald offset model for monochromatic beam */
	PMODEL_RANDOM,  /**< Randomly assigned partialities */
	PMODEL_GGPM     /**< Generalised gaussian prediction model */

} PartialityModel;


/**
 * This structure represents the polarisation of the incident radiation
 */
struct polarisation
{
	double fraction;   /**< Polarisation fraction (0 to 1) */

	double angle;      /**< Angle of electron beam, radians, clockwise from
	                    * horizontal when looking along beam */

	int disable;       /**< If set, don't even correct for the polarisation
	                    * of the diffracted beam (1 - 0.5sin^2(2theta)) */
};


extern RefList *predict_to_res(Crystal *cryst, struct image *image, double max_res);

extern void calculate_partialities(RefList *list,
                                   Crystal *cryst,
                                   struct image *image,
                                   PartialityModel pmodel);

extern void update_predictions(RefList *list, Crystal *cryst, struct image *image);
extern struct polarisation parse_polarisation(const char *text);
extern void polarisation_correction(RefList *list, UnitCell *cell,
                                    struct polarisation p);

extern double sphere_fraction(double rlow, double rhigh, double pr);
extern double gaussian_fraction(double rlow, double rhigh, double pr);

#ifdef __cplusplus
}
#endif

#endif	/* GEOMETRY_H */
