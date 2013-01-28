/*
 * index.h
 *
 * Perform indexing (somehow)
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
 *   2010      Richard Kirian
 *   2012      Lorenzo Galli
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


#include "cell.h"
#include "image.h"
#include "detector.h"


/* Indexing methods */
typedef enum {
	INDEXING_NONE,
	INDEXING_DIRAX,
	INDEXING_MOSFLM,
        INDEXING_XDS,
	INDEXING_REAX,
} IndexingMethod;


/* Cell reduction methods */
enum {
	CELLR_NONE,
	CELLR_REDUCE,
	CELLR_COMPARE,
	CELLR_COMPARE_AB,
};


typedef struct _indexingprivate IndexingPrivate;

extern IndexingPrivate *indexing_private(IndexingMethod indm);

extern IndexingPrivate **prepare_indexing(IndexingMethod *indm, UnitCell *cell,
                                         const char *filename,
                                         struct detector *det,
                                         double nominal_photon_energy);

extern void map_all_peaks(struct image *image);

extern void index_pattern(struct image *image, UnitCell *cell,
                          IndexingMethod *indm, int cellr, int verbose,
                          IndexingPrivate **priv, int config_insane,
                          const float *ltl);

extern void cleanup_indexing(IndexingPrivate **priv);

extern IndexingMethod *build_indexer_list(const char *str, int *need_cell);

#endif	/* INDEX_H */
