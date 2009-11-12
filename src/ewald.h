/*
 * ewald.h
 *
 * Calculate q-vector arrays
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef EWALD_H
#define EWALD_H

#include "image.h"

extern void get_ewald(struct image *image);

#endif	/* EWALD_H */
