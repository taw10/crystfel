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
#include "utils.h"


extern void write_reflections(const char *filename, ReflItemList *items,
                              double *intensities, double *phases,
                              unsigned int *counts, UnitCell *cell,
                              double adu_per_photon);

extern ReflItemList *read_reflections(const char *filename,
                                      double *intensities, double *phases,
                                      unsigned int *counts, double *esds);

extern double *ideal_intensities(double complex *sfac);


#endif	/* REFLECTIONS_H */
