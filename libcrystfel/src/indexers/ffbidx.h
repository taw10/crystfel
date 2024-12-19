/*
 * ffbidx.h
 *
 * Invoke the Fast Feedback Indexer library
 *
 * Copyright © 2023 Paul Scherrer Institute
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Filip Leonarski <filip.leonarski@psi.ch>
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

#ifndef CRYSTFEL_FFBIDX_H
#define CRYSTFEL_FFBIDX_H

#include "index.h"

#include <argp.h>

#ifdef __cplusplus
extern "C" {
#endif

extern int ffbidx_default_options(struct ffbidx_options **opts_ptr);

extern int run_ffbidx(struct image *image, void *ipriv);

extern void *ffbidx_prepare(IndexingMethod *indm, UnitCell *cell, struct ffbidx_options *opts);

extern void ffbidx_cleanup(void *pp);

extern const char *ffbidx_probe(UnitCell *cell);

#ifdef __cplusplus
}
#endif

#endif //CRYSTFEL_FFBIDX_H
