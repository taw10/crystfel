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

extern const char *get_holohedral(const char *sym);

extern ReflItemList *get_twins(ReflItemList *items,
                               const char *holo, const char *mero);


#endif	/* SYMMETRY_H */
