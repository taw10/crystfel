/*
 * intensities.h
 *
 * Extract intensities from patterns
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef INTENSITIES_H
#define INTENSITIES_H

#include "image.h"
#include "cell.h"

extern void output_intensities(struct image *image, UnitCell *cell);

#endif	/* INTENSITIES_H */
