/*
 * diffraction.h
 *
 * Calculate diffraction patterns by Fourier methods
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DIFFRACTION_H
#define DIFFRACTION_H

#include "image.h"
#include "cell.h"
#include "symmetry.h"

typedef enum {
	GRADIENT_MOSAIC,
	GRADIENT_INTERPOLATE,
	GRADIENT_PHASED
} GradientMethod;

extern void get_diffraction(struct image *image, int na, int nb, int nc,
                            const double *intensities, const double *phases,
                            const unsigned char *flags, UnitCell *cell,
                            GradientMethod m, const SymOpList *sym);

#endif	/* DIFFRACTION_H */
