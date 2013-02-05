/*
 * stream.h
 *
 * Stream tools
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
 *   2011      Andrew Aquila
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

#define CHUNK_START_MARKER "----- Begin chunk -----"
#define CHUNK_END_MARKER "----- End chunk -----"
#define PEAK_LIST_START_MARKER "Peaks from peak search"
#define PEAK_LIST_END_MARKER "End of peak list"
#define CRYSTAL_START_MARKER "--- Begin crystal"
#define CRYSTAL_END_MARKER "--- End crystal"
#define REFLECTION_START_MARKER "Reflections measured after indexing"
/* REFLECTION_END_MARKER is over in reflist-utils.h because it is also
 * used to terminate a standalone list of reflections */

typedef struct _stream Stream;

extern Stream *open_stream_for_read(const char *filename);
extern Stream *open_stream_for_write(const char *filename);
extern Stream *open_stream_fd_for_write(int fd);

extern void close_stream(Stream *st);

extern void write_chunk(Stream *st, struct image *image, struct hdfile *hdfile,
                        int include_peaks, int include_reflections);

extern int read_chunk(Stream *st, struct image *image);

extern int rewind_stream(Stream *st);

extern int is_stream(const char *filename);

#endif	/* STREAM_H */
