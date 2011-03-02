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
#include "sfac.h"
#include "cl-utils.h"
#include "beam-parameters.h"


#define SAMPLING (4)
#define BWSAMPLING (10)

#define SINC_LUT_ELEMENTS (4096)


struct gpu_context
{
	cl_context ctx;
	cl_command_queue cq;
	cl_program prog;
	cl_kernel kern;
	cl_mem intensities;
	cl_mem flags;

	cl_mem tt;
	size_t tt_size;

	cl_mem diff;
	size_t diff_size;

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


static int sfloat(struct gpu_context *gctx, int idx, float val)
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


static int setint(struct gpu_context *gctx, int idx, int val)
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


static int setmem(struct gpu_context *gctx, int idx, cl_mem val)
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
	float klow, khigh;
	cl_event *event;
	int p;
	float *tt_ptr;
	int x, y;
	cl_float16 cell;
	float *diff_ptr;
	cl_int4 ncells;
	const int sampling = SAMPLING;
	cl_float bwstep;
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

	/* Calculate wavelength */
	klow = 1.0/(image->lambda*(1.0 + image->beam->bandwidth/2.0));
	khigh = 1.0/(image->lambda*(1.0 - image->beam->bandwidth/2.0));
	bwstep = (khigh-klow) / BWSAMPLING;

	ncells.s[0] = na;
	ncells.s[1] = nb;
	ncells.s[2] = nc;
	ncells.s[3] = 0;  /* unused */

	/* Ensure all required LUTs are available */
	check_sinc_lut(gctx, na);
	check_sinc_lut(gctx, nb);
	check_sinc_lut(gctx, nc);

	if ( setmem(gctx, 0, gctx->diff) ) return;
	if ( setmem(gctx, 1, gctx->tt) ) return;
	if ( setmem(gctx, 9, gctx->intensities) ) return;
	if ( setmem(gctx, 15, gctx->sinc_luts[na-1]) ) return;
	if ( setmem(gctx, 16, gctx->sinc_luts[nb-1]) ) return;
	if ( setmem(gctx, 17, gctx->sinc_luts[nc-1]) ) return;
	if ( setmem(gctx, 18, gctx->flags) ) return;

	/* Unit cell */
	clSetKernelArg(gctx->kern, 8, sizeof(cl_float16), &cell);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't set unit cell: %s\n", clError(err));
		return;
	}

	/* Local memory for reduction */
	clSetKernelArg(gctx->kern, 13,
	               BWSAMPLING*SAMPLING*SAMPLING*sizeof(cl_float), NULL);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't set local memory: %s\n", clError(err));
		return;
	}


	if ( sfloat(gctx, 2, klow) ) return;
	if ( setint(gctx, 3, image->width) ) return;
	if ( setint(gctx, 12, sampling) ) return;
	if ( sfloat(gctx, 14, bwstep) ) return;

	/* Iterate over panels */
	event = malloc(image->det->n_panels * sizeof(cl_event));
	for ( p=0; p<image->det->n_panels; p++ ) {

		size_t dims[3];
		size_t ldims[3] = {SAMPLING, SAMPLING, BWSAMPLING};

		/* In a future version of OpenCL, this could be done
		 * with a global work offset.  But not yet... */
		dims[0] = 1+image->det->panels[p].max_fs
		           -image->det->panels[p].min_fs;
		dims[1] = 1+image->det->panels[p].max_ss
		           -image->det->panels[p].min_ss;
		dims[0] *= SAMPLING;
		dims[1] *= SAMPLING;
		dims[2] = BWSAMPLING;

		if ( sfloat(gctx, 4, image->det->panels[p].cnx) ) return;
		if ( sfloat(gctx, 5, image->det->panels[p].cny) ) return;
		if ( sfloat(gctx, 6, image->det->panels[p].res) ) return;
		if ( sfloat(gctx, 7, image->det->panels[p].clen) ) return;
		if ( setint(gctx, 10, image->det->panels[p].min_fs) ) return;
		if ( setint(gctx, 11, image->det->panels[p].min_ss) ) return;
		if ( sfloat(gctx, 19, image->det->panels[p].fsx) ) return;
		if ( sfloat(gctx, 19, image->det->panels[p].fsy) ) return;
		if ( sfloat(gctx, 20, image->det->panels[p].ssx) ) return;
		if ( sfloat(gctx, 21, image->det->panels[p].ssy) ) return;

		err = clEnqueueNDRangeKernel(gctx->cq, gctx->kern, 3, NULL,
		                             dims, ldims, 0, NULL, &event[p]);
		if ( err != CL_SUCCESS ) {
			ERROR("Couldn't enqueue diffraction kernel: %s\n",
			      clError(err));
			return;
		}
	}

	diff_ptr = clEnqueueMapBuffer(gctx->cq, gctx->diff, CL_TRUE,
	                              CL_MAP_READ, 0, gctx->diff_size,
	                              image->det->n_panels, event, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't map diffraction buffer: %s\n", clError(err));
		return;
	}
	tt_ptr = clEnqueueMapBuffer(gctx->cq, gctx->tt, CL_TRUE, CL_MAP_READ, 0,
	                            gctx->tt_size, image->det->n_panels, event,
	                            NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't map tt buffer\n");
		return;
	}

	free(event);

	image->data = calloc(image->width * image->height, sizeof(float));
	image->twotheta = calloc(image->width * image->height, sizeof(double));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		float val, tt;

		val = diff_ptr[x + image->width*y];
		if ( isinf(val) ) n_inf++;
		if ( val < 0.0 ) n_neg++;
		if ( isnan(val) ) n_nan++;
		tt = tt_ptr[x + image->width*y];

		image->data[x + image->width*y] = val;
		image->twotheta[x + image->width*y] = tt;

	}
	}

	if ( n_neg + n_inf + n_nan ) {
		ERROR("WARNING: The GPU calculation produced %i negative"
		      " values, %i infinities and %i NaNs.\n",
		      n_neg, n_inf, n_nan);
	}

	clEnqueueUnmapMemObject(gctx->cq, gctx->diff, diff_ptr, 0, NULL, NULL);
	clEnqueueUnmapMemObject(gctx->cq, gctx->tt, tt_ptr, 0, NULL, NULL);
}


/* Setup the OpenCL stuff, create buffers, load the structure factor table */
struct gpu_context *setup_gpu(int no_sfac, struct image *image,
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

	/* Create buffer for the picture */
	gctx->diff_size = image->width*image->height*sizeof(cl_float);
	gctx->diff = clCreateBuffer(gctx->ctx, CL_MEM_WRITE_ONLY,
	                            gctx->diff_size, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate diffraction memory\n");
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
	} else {
		ERROR("Sorry!  Point group '%s' is not currently supported"
		      " on the GPU.  I'm using '1' instead.\n", sym);
		strncat(cflags, "-DPG1 ", 511-strlen(cflags));
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

	gctx->tt_size = image->width*image->height*sizeof(cl_float);
	gctx->tt = clCreateBuffer(gctx->ctx, CL_MEM_WRITE_ONLY, gctx->tt_size,
	                          NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate twotheta memory\n");
		free(gctx);
		return NULL;
	}

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
	clReleaseMemObject(gctx->diff);
	clReleaseMemObject(gctx->tt);
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
