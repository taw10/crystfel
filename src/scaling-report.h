/*
 * scaling-report.h
 *
 * Write a nice PDF of scaling parameters
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
