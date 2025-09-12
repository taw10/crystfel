/*
 * mosflm.h
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010      Richard Kirian <rkirian@asu.edu>
 *   2012-2017 Thomas White <taw@physics.org>
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

#ifndef MOSFLM_H
#define MOSFLM_H

#include "index.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file mosflm.h
 * MOSFLM indexer interface
 */

extern int run_mosflm(struct image *image, void *ipriv);

extern void *mosflm_prepare(IndexingMethod indm, UnitCell *cell);
extern const char *mosflm_probe(UnitCell *cell);

extern void mosflm_cleanup(void *pp);

#ifdef __cplusplus
}
#endif

#endif	/* MOSFLM_H */
