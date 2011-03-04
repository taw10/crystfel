/*
 * diffraction.h
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DIFFRACTION_H
#define DIFFRACTION_H

#include "image.h"
#include "cell.h"

typedef enum {
	GRADIENT_MOSAIC,
	GRADIENT_INTERPOLATE,
	GRADIENT_PHASED
} GradientMethod;

extern void get_diffraction(struct image *image, int na, int nb, int nc,
                            const double *intensities, const double *phases,
                            const unsigned char *flags, UnitCell *cell,
                            GradientMethod m, const char *sym);

#endif	/* DIFFRACTION_H */
