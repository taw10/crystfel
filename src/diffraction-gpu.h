/*
 * diffraction-gpu.h
 *
 * Calculate diffraction patterns by Fourier methods (GPU version)
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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

#ifndef DIFFRACTION_GPU_H
#define DIFFRACTION_GPU_H

#include "image.h"
#include "cell.h"

struct gpu_context;

#if HAVE_OPENCL

extern int get_diffraction_gpu(struct gpu_context *gctx, struct image *image,
                               int na, int nb, int nc, UnitCell *ucell,
                               int no_fringes);
extern struct gpu_context *setup_gpu(int no_sfac,
                                     const double *intensities,
                                     const unsigned char *flags,
                                     const char *sym, int dev_num);
extern void cleanup_gpu(struct gpu_context *gctx);

#else

static int get_diffraction_gpu(struct gpu_context *gctx, struct image *image,
                               int na, int nb, int nc, UnitCell *ucell,
                               int no_fringes)
{
	/* Do nothing */
	ERROR("This copy of CrystFEL was not compiled with OpenCL support.\n");
	return 1;
}

static struct gpu_context *setup_gpu(int no_sfac,
                                     const double *intensities,
                                     const unsigned char *flags,
                                     const char *sym, int dev_num)
{
	return NULL;
}

static void cleanup_gpu(struct gpu_context *gctx)
{
}

#endif

#endif	/* DIFFRACTION_GPU_H */
