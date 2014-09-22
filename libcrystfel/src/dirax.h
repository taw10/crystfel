/*
 * dirax.h
 *
 * Invoke the DirAx auto-indexing program
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010,2012-2014 Thomas White <taw@physics.org>
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

#ifndef DIRAX_H
#define DIRAX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "index.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int run_dirax(struct image *image, IndexingPrivate *ipriv);

extern IndexingPrivate *dirax_prepare(IndexingMethod *indm,
                                      UnitCell *cell, struct detector *det,
                                      float *ltl);

extern void dirax_cleanup(IndexingPrivate *pp);

#ifdef __cplusplus
}
#endif

#endif	/* DIRAX_H */
