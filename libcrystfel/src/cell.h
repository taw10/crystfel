/*
 * cell.h
 *
 * A class representing a unit cell
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2012,2014,2017 Thomas White <taw@physics.org>
 *   2010,2012           Richard Kirian
 *   2012                Lorenzo Galli
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

#ifndef CELL_H
#define CELL_H

#include "utils.h"
#include "integer_matrix.h"

/**
 * \file cell.h
 * Unit cell structure
 */

/**
 * Structure representing a 3D vector in reciprocal space.
 *
 * Note: Heavily abused to serve as a real space vector as well.
 **/
struct rvec
{
	double   u;  /**< x component (in reciprocal space) */
	double   v;  /**< y component (in reciprocal space) */
	double   w;  /**< z component (in reciprocal space) */
};


/**
 * An enumeration of the possible lattice types: triclinic, monoclinic,
 * orthorhombic, tetragonal, rhombohedral, hexagonal and cubic.
 **/
typedef enum
{
	L_TRICLINIC,
	L_MONOCLINIC,
	L_ORTHORHOMBIC,
	L_TETRAGONAL,
	L_RHOMBOHEDRAL,
	L_HEXAGONAL,
	L_CUBIC
} LatticeType;


/**
 * Opaque data structure representing a unit cell.
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 **/
typedef struct _unitcell UnitCell;


#ifdef __cplusplus
extern "C" {
#endif

extern UnitCell *cell_new(void);
extern UnitCell *cell_new_from_cell(const UnitCell *orig);
extern void cell_free(UnitCell *cell);

/* Lengths in m, angles in radians */
extern UnitCell *cell_new_from_parameters(double a, double b, double c,
				double alpha, double beta, double gamma);

extern UnitCell *cell_new_from_reciprocal_axes(struct rvec as, struct rvec bs,
                                               struct rvec cs);

extern UnitCell *cell_new_from_direct_axes(struct rvec as, struct rvec bs,
                                           struct rvec cs);

extern int cell_has_parameters(const UnitCell *cell);

extern void cell_set_cartesian(UnitCell *cell,
                               double ax, double ay, double az,
                               double bx, double by, double bz,
                               double cx, double cy, double cz);

extern void cell_set_parameters(UnitCell *cell, double a, double b, double c,
				double alpha, double beta, double gamma);

extern int cell_get_parameters(UnitCell *cell, double *a, double *b, double *c,
                               double *alpha, double *beta, double *gamma);

extern int cell_get_cartesian(UnitCell *cell,
                              double *ax, double *ay, double *az,
                              double *bx, double *by, double *bz,
                              double *cx, double *cy, double *cz);

extern int cell_get_reciprocal(UnitCell *cell,
                               double *asx, double *asy, double *asz,
                               double *bsx, double *bsy, double *bsz,
                               double *csx, double *csy, double *csz);

extern void cell_set_reciprocal(UnitCell *cell,
                                double asx, double asy, double asz,
                                double bsx, double bsy, double bsz,
                                double csx, double csy, double csz);

extern LatticeType cell_get_lattice_type(const UnitCell *cell);
extern void cell_set_lattice_type(UnitCell *cell, LatticeType lattice_type);

struct g6
{
	double A;
	double B;
	double C;
	double D;
	double E;
	double F;
};

extern struct g6 cell_get_G6(UnitCell *cell);

extern char cell_get_centering(const UnitCell *cell);
extern void cell_set_centering(UnitCell *cell, char centering);

extern char cell_get_unique_axis(const UnitCell *cell);
extern void cell_set_unique_axis(UnitCell *cell, char unique_axis);

extern UnitCell *cell_transform_gsl_direct(UnitCell *in, gsl_matrix *m);

extern UnitCell *cell_transform_rational(UnitCell *cell, RationalMatrix *m);
extern UnitCell *cell_transform_rational_inverse(UnitCell *cell, RationalMatrix *m);

extern UnitCell *cell_transform_intmat(UnitCell *cell, IntegerMatrix *m);
extern UnitCell *cell_transform_intmat_inverse(UnitCell *cell, IntegerMatrix *m);

#ifdef __cplusplus
}
#endif

#endif	/* CELL_H */
