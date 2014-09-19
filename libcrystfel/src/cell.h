/*
 * cell.h
 *
 * A class representing a unit cell
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2012,2014 Thomas White <taw@physics.org>
 *   2010,2012      Richard Kirian
 *   2012           Lorenzo Galli
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "utils.h"
#include "integer_matrix.h"

/**
 * rvec:
 * @u: x component (in reciprocal space)
 * @v: y component (in reciprocal space)
 * @w: z component (in reciprocal space)
 *
 * Structure representing a 3D vector in reciprocal space.
 * Note: Heavily abused to serve as a real space vector as well.
 **/
struct rvec
{
	double   u;
	double   v;
	double   w;
};


/**
 * LatticeType:
 * @L_TRICLINIC: Triclinic lattice
 * @L_MONOCLINIC: Monoclinic lattice
 * @L_ORTHORHOMBIC: Orthorhombic lattice
 * @L_TETRAGONAL: Tetragonal lattice
 * @L_RHOMBOHEDRAL: Rhombohedral lattice
 * @L_HEXAGONAL: Hexagonal lattice
 * @L_CUBIC: Cubic lattice
 *
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
 * UnitCell:
 *
 * This data structure is opaque.  You must use the available accessor functions
 * to read and write its contents.
 **/
typedef struct _unitcell UnitCell;


/**
 * UnitCellTransformation:
 *
 * This opaque data structure represents a tranformation of a unit cell, such
 * as a rotation or a centering operation.
 **/
typedef struct _unitcelltransformation UnitCellTransformation;

#ifdef __cplusplus
extern "C" {
#endif

extern UnitCell *cell_new(void);
extern UnitCell *cell_new_from_cell(UnitCell *orig);
extern void cell_free(UnitCell *cell);

/* Lengths in m, angles in radians */
extern UnitCell *cell_new_from_parameters(double a, double b, double c,
				double alpha, double beta, double gamma);

extern UnitCell *cell_new_from_reciprocal_axes(struct rvec as, struct rvec bs,
                                               struct rvec cs);

extern UnitCell *cell_new_from_direct_axes(struct rvec as, struct rvec bs,
                                           struct rvec cs);

extern int cell_has_parameters(UnitCell *cell);

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

extern LatticeType cell_get_lattice_type(UnitCell *cell);
extern void cell_set_lattice_type(UnitCell *cell, LatticeType lattice_type);

extern char cell_get_centering(UnitCell *cell);
extern void cell_set_centering(UnitCell *cell, char centering);

extern char cell_get_unique_axis(UnitCell *cell);
extern void cell_set_unique_axis(UnitCell *cell, char unique_axis);

extern const char *cell_rep(UnitCell *cell);

extern UnitCell *cell_transform(UnitCell *cell, UnitCellTransformation *t);
extern UnitCell *cell_transform_inverse(UnitCell *cell,
                                        UnitCellTransformation *t);

extern UnitCellTransformation *tfn_identity(void);
extern UnitCellTransformation *tfn_from_intmat(IntegerMatrix *m);
extern void tfn_combine(UnitCellTransformation *t,
                        double *na, double *nb, double *nc);
extern void tfn_print(UnitCellTransformation *t);
extern UnitCellTransformation *tfn_inverse(UnitCellTransformation *t);
extern double *tfn_vector(double a, double b, double c);
extern void tfn_free(UnitCellTransformation *t);

#ifdef __cplusplus
}
#endif

#endif	/* CELL_H */
