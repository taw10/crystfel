/*
 * symmetry.h
 *
 * Symmetry
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
 *   2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
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

#ifndef SYMMETRY_H
#define SYMMETRY_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "integer_matrix.h"

/**
 * SymOpList
 *
 * The SymOpList is an opaque data structure containing a list of point symmetry
 * operations.  It could represent an point group or a list of indexing
 * ambiguities (twin laws), or similar.
 **/
typedef struct _symoplist SymOpList;

/**
 * SymOpMask
 *
 * The SymOpMask is an opaque data structure containing a list of flags
 * associated with point symmetry operations in a specific %SymOpList.  It is
 * used to filter the operations in the %SymOpList to avoid duplicating
 * equivalent reflections when the reflection is somehow special (e.g. 'hk0').
 **/
typedef struct _symopmask SymOpMask;

#ifdef __cplusplus
extern "C" {
#endif

extern void free_symoplist(SymOpList *ops);
extern SymOpList *get_pointgroup(const char *sym);

extern SymOpMask *new_symopmask(const SymOpList *list);
extern void free_symopmask(SymOpMask *m);

extern void special_position(const SymOpList *ops, SymOpMask *m,
                             signed int h, signed int k, signed int l);
extern void get_asymm(const SymOpList *ops,
                      signed int h, signed int k, signed int l,
                      signed int *hp, signed int *kp, signed int *lp);
extern int num_equivs(const SymOpList *ops, const SymOpMask *m);
extern void get_equiv(const SymOpList *ops, const SymOpMask *m, int idx,
                      signed int h, signed int k, signed int l,
                      signed int *he, signed int *ke, signed int *le);
extern IntegerMatrix *get_symop(const SymOpList *ops, const SymOpMask *m,
                                int idx);

extern SymOpList *get_ambiguities(const SymOpList *source,
                                  const SymOpList *target);
extern int is_subgroup(const SymOpList *source, const SymOpList *target);

extern int is_centrosymmetric(const SymOpList *s);
extern const char *symmetry_name(const SymOpList *ops);
extern void set_symmetry_name(SymOpList *ops, const char *name);
extern void describe_symmetry(const SymOpList *s);

extern int is_centric(signed int h, signed int k, signed int l,
                      const SymOpList *ops);

extern void add_symop(SymOpList *ops, IntegerMatrix *m);
extern SymOpList *parse_symmetry_operations(const char *s);
extern char *get_matrix_name(const IntegerMatrix *m, int row);

#ifdef __cplusplus
}
#endif

#endif	/* SYMMETRY_H */
