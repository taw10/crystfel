/*
 * index.h
 *
 * Perform indexing (somehow)
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2017 Thomas White <taw@physics.org>
 *   2010      Richard Kirian
 *   2012      Lorenzo Galli
 *   2015      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
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

#ifndef INDEX_H
#define INDEX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * \file index.h
 * The indexing subsystem
 */


#define INDEXING_DEFAULTS_DIRAX (INDEXING_DIRAX)

#define INDEXING_DEFAULTS_ASDF (INDEXING_ASDF | INDEXING_USE_CELL_PARAMETERS)

#define INDEXING_DEFAULTS_MOSFLM (INDEXING_MOSFLM | INDEXING_USE_LATTICE_TYPE  \
                                  | INDEXING_USE_CELL_PARAMETERS)

#define INDEXING_DEFAULTS_FELIX (INDEXING_FELIX | INDEXING_USE_LATTICE_TYPE \
                                     | INDEXING_USE_CELL_PARAMETERS)

#define INDEXING_DEFAULTS_TAKETWO (INDEXING_TAKETWO \
                                   | INDEXING_USE_CELL_PARAMETERS \
                                   | INDEXING_USE_LATTICE_TYPE)

#define INDEXING_DEFAULTS_XDS (INDEXING_XDS | INDEXING_USE_LATTICE_TYPE \
                                     | INDEXING_USE_CELL_PARAMETERS)

#define INDEXING_DEFAULTS_XGANDALF (INDEXING_XGANDALF | INDEXING_USE_CELL_PARAMETERS)

#define INDEXING_DEFAULTS_PINKINDEXER (INDEXING_PINKINDEXER | INDEXING_USE_CELL_PARAMETERS)

/**
 * An enumeration of all the available indexing methods.
 **/
typedef enum {

	INDEXING_NONE   = 0,      /**< No indexing to be performed */

	INDEXING_DIRAX  = 1,      /**< Invoke DirAx program */
	INDEXING_MOSFLM = 2,      /**< Invoke MOSFLM program */
	INDEXING_FELIX = 4,       /**< Invoke Felix program */
	INDEXING_XDS = 5,         /**< Invoke XDS program (NB not nXDS) */
	INDEXING_SIMULATION = 6,  /**< Dummy value for simulated data */
	INDEXING_DEBUG = 7,       /**< Results injector for debugging */
	INDEXING_ASDF = 8,        /**< Use built-in ASDF algorithm */
	INDEXING_TAKETWO = 9,     /**< Use built-in TakeTwo algorithm */
	INDEXING_XGANDALF = 10,   /**< Use XGANDALF (via optional library) */
	INDEXING_PINKINDEXER = 11,/**< Use PinkIndexer (via optional library) */

	INDEXING_ERROR = 255,     /**< Special value for unrecognised indexing
	                           *   engine */

	/** \name Bits which can be set to modify the behaviour of the above
	 *  indexing methods */
	/**@{*/
	/** Use lattice type and centering information */
	INDEXING_USE_LATTICE_TYPE        = 2048,

	/** Use the cell parameters themselves */
	INDEXING_USE_CELL_PARAMETERS     = 4096,
	/**@}*/


} IndexingMethod;

/** This defines the bits in "IndexingMethod" which are used to represent the
 * core of the indexing method */
#define INDEXING_METHOD_MASK (0xff)

/**
 * Flags affecting how the indexing system processes the results from the
 * indexing engine
 */
typedef enum {

	/** Retry indexing if it doesn't work */
	INDEXING_RETRY = 1,

	/** Attempt to index remaining peaks to find more lattices */
	INDEXING_MULTI = 2,

	/** Refine the indexing solution */
	INDEXING_REFINE = 4,

	/* 8, 16 reserved (formerly INDEXING_CHECK_CELL_COMBINATIONS and
	 * INDEXING_CHECK_CELL_AXES respectively) */

	/** Check that the peaks agree with the indexing solution */
	INDEXING_CHECK_PEAKS = 32,

	/** Check that the unit cell agrees with the target cell */
	INDEXING_CHECK_CELL = 64,

} IndexingFlags;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * This is an opaque data structure containing information needed by the
 * indexing system.
 **/
typedef struct _indexingprivate IndexingPrivate;

/* Convert indexing methods to/from text */
extern char *indexer_str(IndexingMethod indm);
extern IndexingMethod get_indm_from_string(const char *method);
extern IndexingMethod get_indm_from_string_2(const char *method, int *err);

#include "detector.h"
#include "cell.h"
#include "image.h"
#include "taketwo.h"
#include "xgandalf.h"
#include "pinkindexer.h"
#include "felix.h"

extern IndexingPrivate *setup_indexing(const char *methods, UnitCell *cell,
                                       struct detector *det,
                                       struct beam_params *beam, float *ltl,
                                       IndexingFlags flags,
                                       struct taketwo_options *ttopts,
                                       struct xgandalf_options *xgandalf_opts,
                                       struct pinkIndexer_options *pinkIndexer_opts,
                                       struct felix_options *felix_opts);

extern const IndexingMethod *indexing_methods(IndexingPrivate *p, int *n);

extern char *detect_indexing_methods(UnitCell *cell);

extern void index_pattern(struct image *image, IndexingPrivate *ipriv);

extern void index_pattern_2(struct image *image, IndexingPrivate *ipriv,
                            int *ping);

extern void index_pattern_3(struct image *image, IndexingPrivate *ipriv,
                            int *ping, char *last_task);

extern void cleanup_indexing(IndexingPrivate *ipriv);

#ifdef __cplusplus
}
#endif

#endif	/* INDEX_H */
