/*
 * reflist-utils.h
 *
 * Utilities to complement the core reflist.c
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef REFLIST_UTILS_H
#define REFLIST_UTILS_H


#include "reflist.h"
#include "cell.h"
#include "symmetry.h"


#define REFLECTION_END_MARKER "End of reflections"


extern void write_reflections_to_file(FILE *fh, RefList *list, UnitCell *cell);

extern int write_reflist(const char *filename, RefList *list, UnitCell *cell);

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

#endif	/* REFLIST_UTILS_H */
