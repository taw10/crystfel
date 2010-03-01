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

struct gpu_context;

#if HAVE_OPENCL

extern void get_diffraction_gpu(struct gpu_context *gctx, struct image *image);
extern struct gpu_context *setup_gpu(int no_sfac, struct image *image,
                                     struct molecule *molecule,
                                     int na, int nb, int nc);
extern void cleanup_gpu(struct gpu_context *gctx);

#else

static void get_diffraction_gpu(struct gpu_context *gctx, struct image *image)
{
	/* Do nothing */
	ERROR("This copy of CrystFEL was not compiled with OpenCL support.\n");
}

static struct gpu_context *setup_gpu(int no_sfac, struct image *image,
                                     struct molecule *molecule,
                                     int na, int nb, int nc)
{
	return NULL;
}

static void cleanup_gpu(struct gpu_context *gctx)
{
}

#endif

#endif	/* DIFFRACTION_GPU_H */
