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

/**
 * SymOpMask
 *
 * Opaque type.
 **/
typedef struct _symopmask SymOpMask;

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

extern SymOpList *get_ambiguities(const SymOpList *source,
                                  const SymOpList *target);
extern int is_subgroup(const SymOpList *source, const SymOpList *target);

extern int is_centrosymmetric(const SymOpList *s);
extern const char *symmetry_name(const SymOpList *ops);
extern void describe_symmetry(const SymOpList *s);

#endif	/* SYMMETRY_H */
