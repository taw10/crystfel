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

#endif	/* SYMMETRY_H */
