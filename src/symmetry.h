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

extern void apply_symmetry(double ref, unsigned int *counts, const char *sym);

#endif	/* SYMMETRY_H */
