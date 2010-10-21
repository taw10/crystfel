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

extern struct cpeak *find_intersections(struct image *image, UnitCell *cell,
                                        double divergence, double bandwidth,
                                        int *n, int output);

extern double partiality(struct image *image, signed int h,
                         signed int k, signed int l);
extern double integrate_all(struct image *image, struct cpeak *cpeaks, int n);


#endif	/* GEOMETRY_H */
