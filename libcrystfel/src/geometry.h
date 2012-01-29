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
#include "cell.h"

#ifdef __cplusplus
extern "C" {
#endif

RefList *find_intersections(struct image *image, UnitCell *cell);

void update_partialities(struct image *image);

#ifdef __cplusplus
}
#endif

#endif	/* GEOMETRY_H */
