/*
 * reflections.h
 *
 * Utilities for handling reflections
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef REFLECTIONS_H
#define REFLECTIONS_H


#include "cell.h"


extern void write_reflections(const char *filename, unsigned int *counts,
                              double *ref, int zone_axis, UnitCell *cell);

extern double *read_reflections(const char *filename, unsigned int *counts);

extern double *ideal_intensities(double complex *sfac);


#endif	/* REFLECTIONS_H */
