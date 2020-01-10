/*
 * cell-utils.h
 *
 * Unit Cell utility functions
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2018 Thomas White <taw@physics.org>
 *   2012      Lorenzo Galli
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

#include <gsl/gsl_matrix.h>

#include "cell.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file cell-utils.h
 * Unit cell utility functions.
 **/


extern double resolution(UnitCell *cell,
                         signed int h, signed int k, signed int l);

extern UnitCell *cell_rotate(UnitCell *in, struct quaternion quat);
extern UnitCell *rotate_cell(UnitCell *in, double omega, double phi,
                             double rot);

extern void cell_print(UnitCell *cell);
extern void cell_print_full(UnitCell *cell);

extern UnitCell *load_cell_from_pdb(const char *filename);
extern UnitCell *load_cell_from_file(const char *filename);
extern void write_cell(UnitCell *cell, FILE *fh);

extern int cell_is_sensible(UnitCell *cell);

extern int validate_cell(UnitCell *cell);

extern UnitCell *uncenter_cell(UnitCell *in, IntegerMatrix **pC,
                               RationalMatrix **pCi);

extern int bravais_lattice(UnitCell *cell);

extern int right_handed(UnitCell *cell);

extern const char *str_lattice(LatticeType l);
extern LatticeType lattice_from_str(const char *s);

extern int forbidden_reflection(UnitCell *cell,
                                signed int h, signed int k, signed int l);

extern double cell_get_volume(UnitCell *cell);

extern int compare_cell_parameters(UnitCell *cell, UnitCell *reference,
                                   const double *tols);

extern int compare_cell_parameters_and_orientation(UnitCell *cell,
                                                   UnitCell *reference,
                                                   const double *tols);

extern int compare_permuted_cell_parameters_and_orientation(UnitCell *cell,
                                                            UnitCell *reference,
                                                            const double *tols,
                                                            IntegerMatrix **pmb);

extern int compare_derivative_cell_parameters(UnitCell *cell, UnitCell *reference,
                                              const double *tols, int csl,
                                              RationalMatrix **pmb);

extern UnitCell *compare_reindexed_cell_parameters(UnitCell *cell_in,
                                                   UnitCell *reference_in,
                                                   const double *tols,
                                                   RationalMatrix **pmb);

#ifdef __cplusplus
}
#endif

#endif	/* CELL_UTILS_H */
