/*
 * cl-utils.c
 *
 * OpenCL utility functions
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cl.h>

#include "utils.h"


const char *clError(cl_int err)
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
	case CL_INVALID_KERNEL_ARGS : return "invalid kernel arguments";
	case CL_INVALID_WORK_GROUP_SIZE : return "invalid work group size";
	case CL_IMAGE_FORMAT_NOT_SUPPORTED : return "image format not supported";
	default :
		ERROR("Error code: %i\n", err);
		return "unknown error";
	}
}


cl_device_id get_first_dev(cl_context ctx)
{
	cl_device_id dev;
	cl_int r;

	r = clGetContextInfo(ctx, CL_CONTEXT_DEVICES, sizeof(dev), &dev, NULL);
	if ( r != CL_SUCCESS ) {
		ERROR("Couldn't get device: %s\n", clError(r));
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


cl_program load_program(const char *filename, cl_context ctx,
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

	r = clBuildProgram(prog, 0, NULL,
	                   "-Werror -I"DATADIR"/crystfel/ -cl-no-signed-zeros",
	                   NULL, NULL);
	if ( r != CL_SUCCESS ) {
		ERROR("Couldn't build program '%s'\n", filename);
		show_build_log(prog, dev);
		*err = r;
		return 0;
	}

	free(source);
	*err = CL_SUCCESS;
	return prog;
}
