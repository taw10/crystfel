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

extern SRContext *sr_header(const char *filename, const char *stream_filename,
                            const char *cmdline);

extern void sr_before(SRContext *sr, struct image *images, int n,
                      RefList *full);

extern void sr_after(SRContext *sr, struct image *images, int n,
                     RefList *full);

#else

SRContext *sr_header(const char *filename, const char *stream_filename,
                     const char *cmdline)
{
	return NULL;
}

void sr_before(SRContext *sr, struct image *images, int n, RefList *full)
{
}

void sr_after(SRContext *sr, struct image *images, int n, RefList *full)
{
}

#endif

#endif	/* SCALING_REPORT_H */
