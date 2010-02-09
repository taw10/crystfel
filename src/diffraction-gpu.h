/*
 * diffraction-gpu.h
 *
 * Calculate diffraction patterns by Fourier methods (GPU version)
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DIFFRACTION_GPU_H
#define DIFFRACTION_GPU_H

#include "image.h"
#include "cell.h"

#if HAVE_OPENCL
extern void get_diffraction_gpu(struct image *image, int na, int nb, int nc,
                                int nosfac);
#else
static void get_diffraction_gpu(struct image *image, int na, int nb, int nc,
                                int nosfac)
{
	/* Do nothing */
	ERROR("This copy of CrystFEL was not compiled with OpenCL support.\n");
}
#endif

#endif	/* DIFFRACTION_GPU_H */
