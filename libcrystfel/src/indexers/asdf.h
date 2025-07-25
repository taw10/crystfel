/*
 * asdf.h
 *
 * Alexandra's Superior Direction Finder, or
 * Algorithm Similar to DirAx, FFT-based
 *
 * Copyright © 2014-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014-2015 Alexandra Tolstikova <alexandra.tolstikova@desy.de>
 *   2015-2021 Thomas White <taw@physics.org>
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

#ifndef ASDF_H
#define ASDF_H

#include <argp.h>
#include "index.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file asdf.h
 * The ASDF indexing algorithm.
 */

extern int run_asdf(struct image *image, void *ipriv);
extern int asdf_default_options(struct asdf_options **opts_ptr);

extern void *asdf_prepare(IndexingMethod indm, UnitCell *cell,
                          struct asdf_options *asdf_opts);
extern const char *asdf_probe(UnitCell *cell);

extern void asdf_cleanup(void *pp);

#ifdef __cplusplus
}
#endif

#endif	/* ASDF_H */
