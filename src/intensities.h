/*
 * intensities.h
 *
 * Extract intensities from patterns
 *
 * (c) 2007-2009 Thomas White <taw@physics.org>
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

extern void output_intensities(struct image *image);

#endif	/* INTENSITIES_H */
