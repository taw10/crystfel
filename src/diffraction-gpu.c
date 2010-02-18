/*
 * diffraction-gpu.c
 *
 * Calculate diffraction patterns by Fourier methods (GPU version)
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <CL/cl.h>

#include "image.h"
#include "utils.h"
#include "cell.h"
#include "diffraction.h"
#include "sfac.h"


#define SAMPLING (5)
#define BWSAMPLING (10)
#define BANDWIDTH (1.0 / 100.0)


static const char *clstrerr(cl_int err)
{
	switch ( err ) {
	case CL_SUCCESS : return "no error";
	case CL_INVALID_PLATFORM : return "invalid platform";
	case CL_INVALID_KERNEL : return "invalid kernel";
	case CL_INVALID_ARG_INDEX : return "invalid argument index";
	case CL_INVALID_ARG_VALUE : return "invalid argument value";
	case CL_INVALID_MEM_OBJECT : return "invalid memory object";
	case CL_INVALID_SAMPLER : return "invalid sampler";
	case CL_INVALID_ARG_SIZE : return "invalid argument size";
	case CL_INVALID_COMMAND_QUEUE  : return "invalid command queue";
	case CL_INVALID_CONTEXT : return "invalid context";
	case CL_INVALID_VALUE : return "invalid value";
	case CL_INVALID_EVENT_WAIT_LIST : return "invalid wait list";
	case CL_MAP_FAILURE : return "map failure";
	case CL_MEM_OBJECT_ALLOCATION_FAILURE : return "object allocation failure";
	case CL_OUT_OF_HOST_MEMORY : return "out of host memory";
	case CL_OUT_OF_RESOURCES : return "out of resources";
	case CL_INVALID_KERNEL_NAME : return "invalid kernel name";
	default :
		ERROR("Error code: %i\n", err);
		return "unknown error";
	}
}


static cl_device_id get_first_dev(cl_context ctx)
{
	cl_device_id dev;
	cl_int r;

	r = clGetContextInfo(ctx, CL_CONTEXT_DEVICES, sizeof(dev), &dev, NULL);
	if ( r != CL_SUCCESS ) {
		ERROR("Couldn't get device\n");
		return 0;
	}

	return dev;
}


static void show_build_log(cl_program prog, cl_device_id dev)
{
	cl_int r;
	char log[4096];
	size_t s;

	r = clGetProgramBuildInfo(prog, dev, CL_PROGRAM_BUILD_LOG, 4096, log,
	                          &s);

	STATUS("%s\n", log);
}


static cl_program load_program(const char *filename, cl_context ctx,
                               cl_device_id dev, cl_int *err)
{
	FILE *fh;
	cl_program prog;
	char *source;
	size_t len;
	cl_int r;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		*err = CL_INVALID_PROGRAM;
		return 0;
	}
	source = malloc(16384);
	len = fread(source, 1, 16383, fh);
	fclose(fh);
	source[len] = '\0';

	prog = clCreateProgramWithSource(ctx, 1, (const char **)&source,
	                                 NULL, err);
	if ( *err != CL_SUCCESS ) {
		ERROR("Couldn't load source\n");
		return 0;
	}

	r = clBuildProgram(prog, 0, NULL, "-Werror", NULL, NULL);
	if ( r != CL_SUCCESS ) {
		ERROR("Couldn't build program '%s'\n", filename);
		show_build_log(prog, dev);
		*err = r;
		return 0;
	}

	*err = CL_SUCCESS;
	return prog;
}


void get_diffraction_gpu(struct image *image, int na, int nb, int nc,
                         int no_sfac)
{
	cl_uint nplat;
	cl_platform_id platforms[8];
	cl_context_properties prop[3];
	cl_context ctx;
	cl_int err;
	cl_command_queue cq;
	cl_program prog;
	cl_device_id dev;
	cl_kernel kern;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double a, b, c, d;
	float kc;
	const size_t dims[2] = {1024, 1024};

	cl_mem sfacs;
	size_t sfac_size;
	float *sfac_ptr;
	cl_mem tt;
	size_t tt_size;
	float *tt_ptr;
	int x, y;
	cl_float16 cell;
	cl_mem diff;
	size_t diff_size;
	float *diff_ptr;
	int i;
	cl_float4 orientation;

	if ( image->molecule == NULL ) return;

	/* Generate structure factors if required */
	if ( !no_sfac ) {
		if ( image->molecule->reflections == NULL ) {
			get_reflections_cached(image->molecule,
			                       ph_lambda_to_en(image->lambda));
		}
	}

	cell_get_cartesian(image->molecule->cell, &ax, &ay, &az,
		                                  &bx, &by, &bz,
		                                  &cx, &cy, &cz);
	cell[0] = ax;  cell[1] = ay;  cell[2] = az;
	cell[3] = bx;  cell[4] = by;  cell[5] = bz;
	cell[6] = cx;  cell[7] = cy;  cell[8] = cz;

	cell_get_parameters(image->molecule->cell,
	                    &a, &b, &c, &d, &d, &d);
	STATUS("Particle size = %i x %i x %i (=%5.2f x %5.2f x %5.2f nm)\n",
	       na, nb, nc, na*a/1.0e-9, nb*b/1.0e-9, nc*c/1.0e-9);

	err = clGetPlatformIDs(8, platforms, &nplat);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't get platform IDs: %i\n", err);
		return;
	}
	if ( nplat == 0 ) {
		ERROR("Couldn't find at least one platform!\n");
		return;
	}
	prop[0] = CL_CONTEXT_PLATFORM;
	prop[1] = (cl_context_properties)platforms[0];
	prop[2] = 0;

	ctx = clCreateContextFromType(prop, CL_DEVICE_TYPE_GPU, NULL, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create OpenCL context: %i\n", err);
		return;
	}

	dev = get_first_dev(ctx);

	cq = clCreateCommandQueue(ctx, dev, 0, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create OpenCL command queue\n");
		return;
	}

	/* Create buffer for the picture */
	diff_size = image->width*image->height*sizeof(cl_float)*2; /* complex */
	diff = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, diff_size, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate diffraction memory\n");
		return;
	}

	/* Create a single-precision version of the scattering factors */
	sfac_size = IDIM*IDIM*sizeof(cl_float)*2; /* complex */
	sfac_ptr = malloc(IDIM*IDIM*sizeof(cl_float)*2);
	for ( i=0; i<IDIM*IDIM; i++ ) {
		sfac_ptr[2*i+0] = creal(image->molecule->reflections[i]);
		sfac_ptr[2*i+1] = cimag(image->molecule->reflections[i]);
	}
	sfacs = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
	                       sfac_size, sfac_ptr, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate sfac memory\n");
		return;
	}

	tt_size = image->width*image->height*sizeof(cl_float);
	tt = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, tt_size, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't allocate twotheta memory\n");
		return;
	}

	prog = load_program(DATADIR"/crystfel/diffraction.cl", ctx, dev, &err);
	if ( err != CL_SUCCESS ) {
		return;
	}

	kern = clCreateKernel(prog, "diffraction", &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't create kernel\n");
		return;
	}

	/* Calculate wavelength */
	kc = 1.0/image->lambda;  /* Centre value */

	/* Orientation */
	orientation[0] = image->orientation.w;
	orientation[1] = image->orientation.x;
	orientation[2] = image->orientation.y;
	orientation[3] = image->orientation.z;

	clSetKernelArg(kern, 0, sizeof(cl_mem), &diff);
	clSetKernelArg(kern, 1, sizeof(cl_mem), &tt);
	clSetKernelArg(kern, 2, sizeof(cl_float), &kc);
	clSetKernelArg(kern, 3, sizeof(cl_int), &image->width);
	clSetKernelArg(kern, 4, sizeof(cl_float), &image->det.panels[0].cx);
	clSetKernelArg(kern, 5, sizeof(cl_float), &image->det.panels[0].cy);
	clSetKernelArg(kern, 6, sizeof(cl_float), &image->det.panels[0].res);
	clSetKernelArg(kern, 7, sizeof(cl_float), &image->det.panels[0].clen);
	clSetKernelArg(kern, 8, sizeof(cl_float16), &cell);
	clSetKernelArg(kern, 9, sizeof(cl_mem), &sfacs);
	clSetKernelArg(kern, 10, sizeof(cl_float4), &orientation);

	err = clEnqueueNDRangeKernel(cq, kern, 2, NULL, dims, NULL,
	                             1, &event_q, &event_d);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't enqueue diffraction kernel\n");
		return;
	}

	diff_ptr = clEnqueueMapBuffer(cq, diff, CL_TRUE, CL_MAP_READ, 0,
	                              diff_size, 1, &event_d, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't map diffraction buffer: %s\n", clstrerr(err));
		return;
	}
	tt_ptr = clEnqueueMapBuffer(cq, tt, CL_TRUE, CL_MAP_READ, 0,
	                            tt_size, 1, &event_d, NULL, &err);
	if ( err != CL_SUCCESS ) {
		ERROR("Couldn't map tt buffer\n");
		return;
	}

	image->sfacs = calloc(image->width * image->height,
	                      sizeof(double complex));
	image->twotheta = calloc(image->width * image->height, sizeof(double));

	for ( x=0; x<image->width; x++ ) {
	for ( y=0; y<image->height; y++ ) {

		float re, im, tt;

		re = diff_ptr[2*(x + image->width*y)+0];
		im = diff_ptr[2*(x + image->width*y)+1];
		tt = tt_ptr[x + image->width*y];

		image->sfacs[x + image->width*y] = re + I*im;
		image->twotheta[x + image->width*y] = tt;

	}
	}

	clReleaseProgram(prog);
	clReleaseMemObject(sfacs);
	clReleaseCommandQueue(cq);
	clReleaseContext(ctx);

}
