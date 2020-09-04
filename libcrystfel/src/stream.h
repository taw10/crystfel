/*
 * stream.h
 *
 * Stream tools
 *
 * Copyright © 2013-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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

/**
 * \file stream.h
 * Stream functions (for indexing results)
 */

struct image;

#include "datatemplate.h"
#include "cell.h"

#define STREAM_GEOM_START_MARKER "----- Begin geometry file -----"
#define STREAM_GEOM_END_MARKER "----- End geometry file -----"
#define STREAM_CELL_START_MARKER "----- Begin unit cell -----"
#define STREAM_CELL_END_MARKER "----- End unit cell -----"
#define STREAM_CHUNK_START_MARKER "----- Begin chunk -----"
#define STREAM_CHUNK_END_MARKER "----- End chunk -----"
#define STREAM_PEAK_LIST_START_MARKER "Peaks from peak search"
#define STREAM_PEAK_LIST_END_MARKER "End of peak list"
#define STREAM_CRYSTAL_START_MARKER "--- Begin crystal"
#define STREAM_CRYSTAL_END_MARKER "--- End crystal"
#define STREAM_REFLECTION_START_MARKER "Reflections measured after indexing"
#define STREAM_REFLECTION_END_MARKER "End of reflections"

/**
 * An opaque structure representing a stream being read or written
 */
typedef struct _stream Stream;

/**
 * A bitfield of things that can be read from or written to a stream.
 * Use this together with stream_{read,write}_chunk to read/write the
 * stream faster if you don't need all the information.
 *
 * General information about crystals (including unit cell parameters)
 * is always read and written.
 **/
typedef enum {

	/** Read the integrated reflections */
	STREAM_REFLECTIONS = 2,

	/** Read the peak search results */
	STREAM_PEAKS = 4,

} StreamFlags;

#ifdef __cplusplus
extern "C" {
#endif

/* Opening/closing streams */
extern Stream *stream_open_for_read(const char *filename);
extern Stream *stream_open_for_write(const char *filename,
                                     const DataTemplate *dtempl);
extern Stream *stream_open_fd_for_write(int fd,
                                        const DataTemplate *dtempl);
extern void stream_close(Stream *st);

/* Writing things to stream header */
extern void stream_write_geometry_file(Stream *st,
                                       const char *geom_filename);
extern void stream_write_target_cell(Stream *st,
                                     const UnitCell *cell);
extern void stream_write_commandline_args(Stream *st,
                                          int argc, char *argv[]);
extern void stream_write_indexing_methods(Stream *st,
                                          const char *indm_str);

/* Metadata */
extern int stream_has_old_indexers(Stream *st);
extern char *stream_audit_info(Stream *st);
extern char *stream_geometry_file(Stream *st);

/* Low-level stuff used for indexamajig sandbox */
extern int stream_get_fd(Stream *st);
extern int stream_rewind(Stream *st);

/* Random access */
typedef struct _streamindex StreamIndex;
extern StreamIndex *stream_make_index(const char *filename);
extern int stream_select_chunk(Stream *st, StreamIndex *index,
                               const char *filename,
                               const char *ev);
extern void stream_index_free(StreamIndex *index);

/* Read/write chunks */
extern struct image *stream_read_chunk(Stream *st, StreamFlags srf);
extern int stream_write_chunk(Stream *st, const struct image *image,
                              StreamFlags srf);

#ifdef __cplusplus
}
#endif

#endif	/* STREAM_H */
