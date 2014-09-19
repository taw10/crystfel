/*
 * cell-utils.h
 *
 * Unit Cell utility functions
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2013,2014 Thomas White <taw@physics.org>
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

#ifndef CELL_UTILS_H
#define CELL_UTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cell.h"

#ifdef __cplusplus
extern "C" {
#endif

extern double resolution(UnitCell *cell,
                         signed int h, signed int k, signed int l);

extern UnitCell *cell_rotate(UnitCell *in, struct quaternion quat);
extern UnitCell *rotate_cell(UnitCell *in, double omega, double phi,
                             double rot);

extern void cell_print(UnitCell *cell);

extern UnitCell *match_cell(UnitCell *cell, UnitCell *tempcell, int verbose,
                            const float *ltl, int reduce);

extern UnitCell *match_cell_ab(UnitCell *cell, UnitCell *tempcell);

extern UnitCell *load_cell_from_pdb(const char *filename);
extern UnitCell *load_cell_from_file(const char *filename);

extern int cell_is_sensible(UnitCell *cell);

extern int validate_cell(UnitCell *cell);

extern UnitCell *uncenter_cell(UnitCell *in, UnitCellTransformation **t);

extern int bravais_lattice(UnitCell *cell);

extern int right_handed(UnitCell *cell);

extern const char *str_lattice(LatticeType l);
extern LatticeType lattice_from_str(const char *s);

extern int forbidden_reflection(UnitCell *cell,
                                signed int h, signed int k, signed int l);

#ifdef __cplusplus
}
#endif

#endif	/* CELL_UTILS_H */
