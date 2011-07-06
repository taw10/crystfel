/*
 * geometry.h
 *
 * Geometry of diffraction
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "reflist.h"

extern RefList *find_intersections(struct image *image, UnitCell *cell);

extern void predict_corresponding_reflections(struct image *image,
                                              const char *sym, int *n_expected,
                                              int *n_found, int *n_notfound);

extern void update_partialities(struct image *image);


#endif	/* GEOMETRY_H */
