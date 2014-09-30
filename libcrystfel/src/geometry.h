/*
 * geometry.h
 *
 * Geometry of diffraction
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "reflist.h"
#include "cell.h"
#include "crystal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * PartialityModel:
 * @PMODEL_UNITY : Set all all partialities and Lorentz factors to 1.
 * @PMODEL_SCSPHERE : Sphere model with source coverage factor included
 * @PMODEL_SCGAUSSIAN : Gaussian model with source coverage factor included
 *
 * A %PartialityModel describes a geometrical model which can be used to
 * calculate spot partialities and Lorentz correction factors.
 **/
typedef enum {

	PMODEL_UNITY,
	PMODEL_SCSPHERE,
	PMODEL_SCGAUSSIAN,

} PartialityModel;

extern RefList *find_intersections(struct image *image, Crystal *cryst,
                                   PartialityModel pmodel);

/* Deprecated: select reflections using Kirian-style pixel proximity */
extern RefList *select_intersections(struct image *image, Crystal *cryst);

extern void update_partialities(Crystal *cryst, PartialityModel pmodel);
extern void update_partialities_2(Crystal *cryst, PartialityModel pmodel,
                                  int *n_gained, int *n_lost,
                                  double *mean_p_change);
extern void polarisation_correction(RefList *list, UnitCell *cell,
                                    struct image *image);

extern double sphere_fraction(double rlow, double rhigh, double pr);
extern double gaussian_fraction(double rlow, double rhigh, double pr);

#ifdef __cplusplus
}
#endif

#endif	/* GEOMETRY_H */
