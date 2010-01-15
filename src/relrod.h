/*
 * relrod.h
 *
 * Calculate reflection positions via line-sphere intersection test
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef RELROD_H
#define RELROD_H

#include "image.h"
#include "cell.h"

extern void get_reflections(struct image *image, UnitCell *cell, double smax);

#endif	/* RELROD_H */
