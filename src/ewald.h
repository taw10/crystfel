/*
 * ewald.h
 *
 * Calculate q-vector arrays
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
