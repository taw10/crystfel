/*
 * xgandalf.h
 *
 * Interface to XGANDALF indexer
 *
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017-2018 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
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

#ifndef LIBCRYSTFEL_SRC_XGANDALF_H
#define LIBCRYSTFEL_SRC_XGANDALF_H

#include <stddef.h>
#include <argp.h>

/**
 * \file xgandalf.h
 * XGANDALF indexer interface
 */

#include "index.h"

extern int xgandalf_default_options(struct xgandalf_options **opts_ptr);

extern int run_xgandalf(struct image *image, void *ipriv);

extern void *xgandalf_prepare(IndexingMethod indm, UnitCell *cell,
                              struct xgandalf_options *xgandalf_opts);

extern void xgandalf_cleanup(void *pp);
extern const char *xgandalf_probe(UnitCell *cell);


#endif /* LIBCRYSTFEL_SRC_XGANDALF_H */
