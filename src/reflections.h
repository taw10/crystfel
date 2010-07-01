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
                              double *ref, double *phases, int zone_axis,
                              UnitCell *cell, unsigned int min_counts);

extern double *read_reflections(const char *filename, unsigned int *counts,
                                double *phases);

extern double *ideal_intensities(double complex *sfac);

extern void divide_down(double *intensities, unsigned int *counts);


#endif	/* REFLECTIONS_H */
