/*
 * stream.h
 *
 * Stream tools
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

#ifndef STREAM_H
#define STREAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


struct image;
struct hdfile;

/* Possible options dictating what goes into the output stream */
enum
{
	STREAM_NONE                 = 0,
	STREAM_INTEGRATED           = 1<<0,
	STREAM_PEAKS                = 1<<2,
	STREAM_PEAKS_IF_INDEXED     = 1<<3,
	STREAM_PEAKS_IF_NOT_INDEXED = 1<<4,
};


extern int count_patterns(FILE *fh);

extern void write_stream_header(FILE *ofh, int argc, char *argv[]);

extern void write_chunk(FILE *ofh, struct image *image, struct hdfile *hdfile,
                        int flags);

extern int parse_stream_flags(const char *a);

extern int read_chunk(FILE *fh, struct image *image);

extern int skip_some_files(FILE *fh, int n);

extern int is_stream(const char *filename);

#endif	/* STREAM_H */
