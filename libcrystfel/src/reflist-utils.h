/*
 * reflist-utils.h
 *
 * Utilities to complement the core reflist.c
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2014 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef REFLIST_UTILS_H
#define REFLIST_UTILS_H

#include "image.h"
#include "reflist.h"
#include "cell.h"
#include "symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif

#define REFLECTION_END_MARKER "End of reflections"


extern void write_reflections_to_file(FILE *fh, RefList *list);

extern int write_reflist(const char *filename, RefList *list);
extern int write_reflist_2(const char *filename, RefList *list, SymOpList *sym);

extern RefList *read_reflections_from_file(FILE *fh);
extern RefList *read_reflections(const char *filename);

extern int check_list_symmetry(RefList *list, const SymOpList *sym);
extern int find_equiv_in_list(RefList *list, signed int h, signed int k,
                              signed int l, const SymOpList *sym, signed int *hu,
                              signed int *ku, signed int *lu);

extern RefList *asymmetric_indices(RefList *in, const SymOpList *sym);

extern void resolution_limits(RefList *list, UnitCell *cell,
                              double *rmin, double *rmax);

extern double max_intensity(RefList *list);

extern RefList *res_cutoff(RefList *list, UnitCell *cell,
                           double min, double max);

extern RefList *copy_reflist(RefList *list);

#ifdef __cplusplus
}
#endif

#endif	/* REFLIST_UTILS_H */
