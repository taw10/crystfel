/*
 * symmetry.h
 *
 * Symmetry
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef SYMMETRY_H
#define SYMMETRY_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * SymOpList
 *
 * Opaque type.
 **/
typedef struct _symoplist SymOpList;

extern void free_symoplist(SymOpList *ops);

extern SymOpList *get_pointgroup(const char *sym);

extern SymOpList *special_position(SymOpList *ops,
                                   signed int h, signed int k, signed int l);

extern void get_asymm(signed int h, signed int k, signed int l,
                      signed int *hp, signed int *kp, signed int *lp,
                      const char *sym);

extern int num_equivs(signed int h, signed int k, signed int l,
                      const char *sym);

extern int num_general_equivs(const char *sym);

extern void get_equiv(signed int h, signed int k, signed int l,
                      signed int *he, signed int *ke, signed int *le,
                      const char *sym, int idx);

extern void get_general_equiv(signed int h, signed int k, signed int l,
                              signed int *he, signed int *ke, signed int *le,
                              const char *sym, int idx);

extern SymOpList *get_twins(SymOpList *source, SymOpList *target);

#endif	/* SYMMETRY_H */
