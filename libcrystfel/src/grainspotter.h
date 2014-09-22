/*
 * grainspotter.h
 *
 * Invoke GrainSpotter for multi-crystal autoindexing
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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

#ifndef GRAINSPOTTER_H
#define GRAINSPOTTER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cell.h"

#ifdef __cplusplus
extern "C" {
#endif

extern IndexingPrivate *grainspotter_prepare(IndexingMethod *indm,
                                             UnitCell *cell,
                                             struct detector *det,
                                             float *ltl);

extern void grainspotter_cleanup(IndexingPrivate *pp);

extern int grainspotter_index(struct image *image, IndexingPrivate *p);

#ifdef __cplusplus
}
#endif

#endif	/* GRAINSPOTTER_H */
