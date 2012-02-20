/*
 * diffraction-gpu.c
 *
 * Calculate diffraction patterns by Fourier methods (GPU version)
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#ifdef HAVE_CL_CL_H
#include <CL/cl.h>
#else
#include <cl.h>
#endif

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "diffraction.h"
#include "cl-utils.h"
#include "beam-parameters.h"
#include "pattern_sim.h"


#define SINC_LUT_ELEMENTS (4096)


struct gpu_context
{
	cl_context ctx;
	cl_command_queue cq;
	cl_program prog;
	cl_kernel kern;
	cl_mem intensities;
	cl_mem flags;

	/* Array of sinc LUTs */
	cl_mem *sinc_luts;
	cl_float **sinc_lut_ptrs;
	int max_sinc_lut;  /* Number of LUTs, i.e. one greater than the maximum
	                    * index.  This equals the highest allowable "n". */
};


static void check_sinc_lut(struct gpu_context *gctx, int n)
{
	cl_int err;
	cl_image_format fmt;
	int i;

	if ( n > gctx->max_sinc_lut ) {

		gctx->sinc_luts = realloc(gctx->sinc_luts,
		                          n*sizeof(*gctx->sinc_luts));
		gctx->sinc_lut_ptrs = realloc(gctx->sinc_lut_ptrs,
		                              n*sizeof(*gctx->sinc_lut_ptrs));

		for ( i=gctx->max_sinc_lut; i<n; i++ ) {
			gctx->sinc_lut_ptrs[i] = NULL;
		}

		gctx->max_sinc_lut = n;
	}

	if ( gctx->sinc_lut_ptrs[n-1] != NULL ) return;

	/* Create a new sinc LUT */
	gctx->sinc_lut_ptrs[n-1] = malloc(SINC_LUT_ELEMENTS*sizeof(cl_float));
	gctx->sinc_lut_ptrs[n-1][0] = n;
	if ( n == 1 ) {
		for ( i=1; i<SINC_LUT_ELEMENTS; i++ ) {
			gctx->sinc_lut_ptrs[n-1][i] = 1.0;
		}
	} else {
		for ( i=1; i<SINC_LUT_ELEMENTS; i++ ) {
			double x, val;
			x = (double)i/SINC_LUT_ELEMENTS;
			val = fabs(sin(M_PI*n*x)/sin(M_PI*x));
			gctx->sinc_lut_ptrs[n-1][i] = val;
		}
	}

	fmt.image_channel_order = CL_INTENSITY;
	fmt.image_channel_data_type = CL_FLOAT;

	gctx->sinc_luts[n-1] = clCreateImage2D(gctx->ctx,
	                                CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	                                &fmt, SINC_LUT_ELEMENTS, 1, 0,
	                                gctx->sinc_lut_ptrs[n-1], &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create LUT for %i\n", n);
		return;
	}
}


static int set_arg_float(struct gpu_context *gctx, int idx, float val)
{
	cl_int err;
	err = clSetKernelArg(gctx->kern, idx, sizeof(cl_float), &val);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't set kernel argument %i: %s\n",
		      idx, clError(err));
		return 1;
	}

	return 0;
}


static int set_arg_int(struct gpu_context *gctx, int idx, int val)
{
	cl_int err;

	err = clSetKernelArg(gctx->kern, idx, sizeof(cl_int), &val);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't set kernel argument %i: %s\n",
		      idx, clError(err));
		return 1;
	}

	return 0;
}


static int set_arg_mem(struct gpu_context *gctx, int idx, cl_mem val)
{
	cl_int err;

	err = clSetKernelArg(gctx->kern, idx, sizeof(cl_mem), &val);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't set kernel argument %i: %s\n",
		      idx, clError(err));
		return 1;
	}

	return 0;
}


void get_diffraction_gpu(struct gpu_context *gctx, struct image *image,
                         int na, int nb, int nc, UnitCell *ucell)
{
	cl_int err;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	int i;
	cl_float16 cell;
	cl_int4 ncells;
	int n_inf = 0;
	int n_neg = 0;
	int n_nan = 0;

	if ( gctx == NULL ) {
		ERROR("GPU setup failed.\n");
		return;
	}

	cell_get_cartesian(ucell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	cell.s[0] = ax;  cell.s[1] = ay;  cell.s[2] = az;
	cell.s[3] = bx;  cell.s[4] = by;  cell.s[5] = bz;
	cell.s[6] = cx;  cell.s[7] = cy;  cell.s[8] = cz;

	ncells.s[0] = na;
	ncells.s[1] = nb;
	ncells.s[2] = nc;
	ncells.s[3] = 0;  /* unused */

	/* Ensure all required LUTs are available */
	check_sinc_lut(gctx, na);
	check_sinc_lut(gctx, nb);
	check_sinc_lut(gctx, nc);

	if ( set_arg_float(gctx, 2, 1.0/image->lambda) ) return;
	if ( set_arg_mem(gctx, 9, gctx->intensities) ) return;
	if ( set_arg_mem(gctx, 10, gctx->sinc_luts[na-1]) ) return;
	if ( set_arg_mem(gctx, 11, gctx->sinc_luts[nb-1]) ) return;
	if ( set_arg_mem(gctx, 12, gctx->sinc_luts[nc-1]) ) return;
	if ( set_arg_mem(gctx, 13, gctx->flags) ) return;

	/* Unit cell */
	err = clSetKernelArg(gctx->kern, 8, sizeof(cl_float16), &cell);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't set unit cell: %s\n", clError(err));
		return;
	}

	/* Allocate memory for the result */
	image->data = calloc(image->width * image->height, sizeof(float));
	image->twotheta = calloc(image->width * image->height, sizeof(double));

	/* Iterate over panels */
	for ( i=0; i<image->det->n_panels; i++ ) {

		size_t dims[2];
		size_t ldims[2] = {1, 1};
		struct panel *p;
		cl_mem tt;
		size_t tt_size;
		cl_mem diff;
		size_t diff_size;
		float *diff_ptr;
		float *tt_ptr;
		int pan_width, pan_height;
		int fs, ss;

		p = &image->det->panels[i];

		pan_width = 1 + p->max_fs - p->min_fs;
		pan_height = 1 + p->max_ss - p->min_ss;

		/* Buffer for the results of this panel */
		diff_size = pan_width * pan_height * sizeof(cl_float);
		diff = clCreateBuffer(gctx->ctx, CL_MEM_WRITE_ONLY,
	                              diff_size, NULL, &err);
		if ( err != CL_SUCCESS ) {
			ERROR("Couldn't allocate diffraction memory\n");
			return;
		}
		tt_size = pan_width * pan_height * sizeof(cl_float);
		tt = clCreateBuffer(gctx->ctx, CL_MEM_WRITE_ONLY, tt_size,
			            NULL, &err);
		if ( err != CL_SUCCESS ) {
			ERROR("Couldn't allocate twotheta memory\n");
			return;
		}

		if ( set_arg_mem(gctx, 0, diff) ) return;
		if ( set_arg_mem(gctx, 1, tt) ) return;
		if ( set_arg_int(gctx, 3, pan_width) ) return;
		if ( set_arg_float(gctx, 4, p->cnx) ) return;
		if ( set_arg_float(gctx, 5, p->cny) ) return;
		if ( set_arg_float(gctx, 6, p->res) ) return;
		if ( set_arg_float(gctx, 7, p->clen) ) return;
		if ( set_arg_float(gctx, 14, p->fsx) ) return;
		if ( set_arg_float(gctx, 15, p->fsy) ) return;
		if ( set_arg_float(gctx, 16, p->ssx) ) return;
		if ( set_arg_float(gctx, 17, p->ssy) ) return;

		dims[0] = pan_width;
		dims[1] = pan_height;

		err = clEnqueueNDRangeKernel(gctx->cq, gctx->kern, 2, NULL,
		                             dims, ldims, 0, NULL, NULL);
		if ( err != CL_SUCCESS ) {
			ERROR("Couldn't enqueue diffraction kernel: %s\n",
			      clError(err));
			return;
		}

		clFinish(gctx->cq);

		diff_ptr = clEnqueueMapBuffer(gctx->cq, diff, CL_TRUE,
		                              CL_MAP_READ, 0, diff_size,
		                              0, NULL, NULL, &err);
		if ( err != CL_SUCCESS ) {
			ERROR("Couldn't map diffraction buffer: %s\n",
			      clError(err));
			return;
		}
		tt_ptr = clEnqueueMapBuffer(gctx->cq, tt, CL_TRUE, CL_MAP_READ,
		                            0, tt_size, 0, NULL, NULL, &err);
		if ( err != CL_SUCCESS ) {
			ERROR("Couldn't map tt buffer\n");
			return;
		}

		for ( fs=0; fs<pan_width; fs++ ) {
		for ( ss=0; ss<pan_height; ss++ ) {

			float val, tt;
			int tfs, tss;

			val = diff_ptr[fs + pan_width*ss];
			if ( isinf(val) ) n_inf++;
			if ( val < 0.0 ) n_neg++;
			if ( isnan(val) ) n_nan++;
			tt = tt_ptr[fs + pan_width*ss];

			tfs = p->min_fs + fs;
			tss = p->min_ss + ss;
			image->data[tfs + image->width*tss] = val;
			image->twotheta[tfs + image->width*tss] = tt;

		}
		}

		clEnqueueUnmapMemObject(gctx->cq, diff, diff_ptr,
		                        0, NULL, NULL);
		clEnqueueUnmapMemObject(gctx->cq, tt, tt_ptr,
		                        0, NULL, NULL);

		clReleaseMemObject(diff);
		clReleaseMemObject(tt);

	}


	if ( n_neg + n_inf + n_nan ) {
		ERROR("WARNING: The GPU calculation produced %i negative"
		      " values, %i infinities and %i NaNs.\n",
		      n_neg, n_inf, n_nan);
	}

}


/* Setup the OpenCL stuff, create buffers, load the structure factor table */
struct gpu_context *setup_gpu(int no_sfac,
                              const double *intensities, unsigned char *flags,
                              const char *sym, int dev_num)
{
	struct gpu_context *gctx;
	cl_uint nplat;
	cl_platform_id platforms[8];
	cl_context_properties prop[3];
	cl_int err;
	cl_device_id dev;
	size_t intensities_size;
	float *intensities_ptr;
	size_t flags_size;
	float *flags_ptr;
	size_t maxwgsize;
	int i;
	char cflags[512] = "";

	STATUS("Setting up GPU...\n");

	err = clGetPlatformIDs(8, platforms, &nplat);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't get platform IDs: %i\n", err);
		return NULL;
	}
	if ( nplat == 0 ) {
		ERROR("Couldn't find at least one platform!\n");
		return NULL;
	}
	prop[0] = CL_CONTEXT_PLATFORM;
	prop[1] = (cl_context_properties)platforms[0];
	prop[2] = 0;

	gctx = malloc(sizeof(*gctx));
	gctx->ctx = clCreateContextFromType(prop, CL_DEVICE_TYPE_GPU,
	                                    NULL, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create OpenCL context: %i\n", err);
		free(gctx);
		return NULL;
	}

	dev = get_cl_dev(gctx->ctx, dev_num);

	gctx->cq = clCreateCommandQueue(gctx->ctx, dev, 0, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create OpenCL command queue\n");
		free(gctx);
		return NULL;
	}

	/* Create a single-precision version of the scattering factors */
	intensities_size = IDIM*IDIM*IDIM*sizeof(cl_float);
	intensities_ptr = malloc(intensities_size);
	if ( intensities != NULL ) {
		for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
			intensities_ptr[i] = intensities[i];
		}
	} else {
		for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
			intensities_ptr[i] = 1e5;
		}
		strncat(cflags, "-DFLAT_INTENSITIES ", 511-strlen(cflags));
	}
	gctx->intensities = clCreateBuffer(gctx->ctx,
	                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	                             intensities_size, intensities_ptr, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate intensities memory\n");
		free(gctx);
		return NULL;
	}
	free(intensities_ptr);

	if ( sym != NULL ) {
		if ( strcmp(sym, "1") == 0 ) {
			strncat(cflags, "-DPG1 ", 511-strlen(cflags));
		} else if ( strcmp(sym, "-1") == 0 ) {
			strncat(cflags, "-DPG1BAR ", 511-strlen(cflags));
		} else if ( strcmp(sym, "6/mmm") == 0 ) {
			strncat(cflags, "-DPG6MMM ", 511-strlen(cflags));
		} else if ( strcmp(sym, "6") == 0 ) {
			strncat(cflags, "-DPG6 ", 511-strlen(cflags));
		} else if ( strcmp(sym, "6/m") == 0 ) {
			strncat(cflags, "-DPG6M ", 511-strlen(cflags));
		} else if ( strcmp(sym, "mmm") == 0 ) {
			strncat(cflags, "-DPGMMM ", 511-strlen(cflags));
		} else {
			ERROR("Sorry!  Point group '%s' is not currently"
			      " supported on the GPU."
			      " I'm using '1' instead.\n", sym);
			strncat(cflags, "-DPG1 ", 511-strlen(cflags));
		}
	} else {
		if ( intensities != NULL ) {
			ERROR("You gave me an intensities file but no point"
			      " group.  I'm assuming '1'.\n");
			strncat(cflags, "-DPG1 ", 511-strlen(cflags));
		}
	}

	/* Create a flag array */
	flags_size = IDIM*IDIM*IDIM*sizeof(cl_float);
	flags_ptr = malloc(flags_size);
	if ( flags != NULL ) {
		for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
			flags_ptr[i] = flags[i];
		}
	} else {
		for ( i=0; i<IDIM*IDIM*IDIM; i++ ) {
			flags_ptr[i] = 1.0;
		}
	}
	gctx->flags = clCreateBuffer(gctx->ctx,
	                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	                             flags_size, flags_ptr, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate flag buffer\n");
		free(gctx);
		return NULL;
	}
	free(flags_ptr);

	gctx->prog = load_program(DATADIR"/crystfel/diffraction.cl", gctx->ctx,
	                          dev, &err, cflags);
	if ( err != CL_SUCCESS ) {
		free(gctx);
		return NULL;
	}

	gctx->kern = clCreateKernel(gctx->prog, "diffraction", &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create kernel\n");
		free(gctx);
		return NULL;
	}

	gctx->max_sinc_lut = 0;
	gctx->sinc_lut_ptrs = NULL;
	gctx->sinc_luts = NULL;

	clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_GROUP_SIZE,
	                sizeof(size_t), &maxwgsize, NULL);
	STATUS("Maximum work group size = %lli\n", (long long int)maxwgsize);

	return gctx;
}


void cleanup_gpu(struct gpu_context *gctx)
{
	int i;

	clReleaseProgram(gctx->prog);
	clReleaseMemObject(gctx->intensities);

	/* Release LUTs */
	for ( i=1; i<=gctx->max_sinc_lut; i++ ) {
		if ( gctx->sinc_lut_ptrs[i-1] != NULL ) {
			clReleaseMemObject(gctx->sinc_luts[i-1]);
			free(gctx->sinc_lut_ptrs[i-1]);
		}
	}

	free(gctx->sinc_luts);
	free(gctx->sinc_lut_ptrs);

	clReleaseCommandQueue(gctx->cq);
	clReleaseContext(gctx->ctx);
	free(gctx);
}
