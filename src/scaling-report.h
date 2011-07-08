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

#ifdef HAVE_CAIRO
extern void scaling_report(const char *filename, const struct image *images,
                           int n, const char *stream_filename);
#else
static inline void scaling_report(const char *filename,
                                  const struct image *images, int n,
                                  const char *stream_filename)
{
	ERROR("Not writing scaling report - no Cairo support.\n");
}
#endif

#endif	/* SCALING_REPORT_H */
