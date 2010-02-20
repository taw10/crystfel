/*
 * cl-utils.h
 *
 * OpenCL utility functions
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef CLUTILS_H
#define CLUTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


extern const char *clError(cl_int err);
extern cl_device_id get_first_dev(cl_context ctx);
extern cl_program load_program(const char *filename, cl_context ctx,
                        cl_device_id dev, cl_int *err);


#endif	/* CLUTILS_H */
