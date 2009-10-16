/*
 * relrod.h
 *
 * Calculate reflection positions via line-sphere intersection test
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
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
