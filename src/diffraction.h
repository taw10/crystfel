/*
 * diffraction.h
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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
                            const double *intensities,
                            const unsigned int *counts, UnitCell *cell,
                            int do_water, GradientMethod m);
extern struct rvec get_q(struct image *image, unsigned int xs, unsigned int ys,
                         unsigned int sampling, float *ttp, float k);

#endif	/* DIFFRACTION_H */
