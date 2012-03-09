/*
 * scaling-report.h
 *
 * Write a nice PDF of scaling parameters
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

#ifndef SCALING_REPORT_H
#define SCALING_REPORT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "utils.h"

typedef struct _srcontext SRContext;  /* Opaque */

#ifdef HAVE_CAIRO

extern SRContext *sr_titlepage(struct image *images, int n,
                               const char *filename,
                               const char *stream_filename,
                               const char *cmdline);

extern void sr_iteration(SRContext *sr, int iteration, struct image *images,
                         int n, RefList *full);

extern void sr_finish(SRContext *sr);

#else

SRContext *sr_titlepage(struct image *images, int n, const char *filename,
                        const char *stream_filename, const char *cmdline)
{
	return NULL;
}

void sr_iteration(SRContext *sr, int iteration, struct image *images, int n,
                  RefList *full)
{
}

void sr_finish(SRContext *sr)
{
}

#endif


#endif	/* SCALING_REPORT_H */
