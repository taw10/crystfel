/*
 * stream.h
 *
 * Indexed stream tools
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef STREAM_H
#define STREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


struct image;

extern int count_patterns(FILE *fh);
extern int find_chunk(FILE *fh, UnitCell **cell, char **filename, double *ev);
extern void write_chunk(FILE *ofh, struct image *image, int flags);


#endif	/* STREAM_H */
