/*
 * stream.h
 *
 * Stream tools
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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

/* Possible options dictating what goes into the output stream */
enum
{
	STREAM_NONE             = 0,
	STREAM_INTEGRATED       = 1<<0,
	STREAM_PIXELS           = 1<<1,
	STREAM_PEAKS            = 1<<2,
	STREAM_PEAKS_IF_INDEXED = 1<<3,
};


extern int count_patterns(FILE *fh);

extern void write_chunk(FILE *ofh, struct image *image, int flags);

extern int parse_stream_flags(const char *a);

extern int read_chunk(FILE *fh, struct image *image);

extern int skip_some_files(FILE *fh, int n);

#endif	/* STREAM_H */
