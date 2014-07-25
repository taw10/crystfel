/*
 * cl-utils.h
 *
 * OpenCL utility functions
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

#ifndef CLUTILS_H
#define CLUTILS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


extern const char *clError(cl_int err);
extern cl_device_id get_cl_dev(cl_context ctx, int n);
extern cl_program load_program(const char *filename, cl_context ctx,
                        cl_device_id dev, cl_int *err,
                        const char *extra_cflags, const char *insert_stuff);


#endif	/* CLUTILS_H */
