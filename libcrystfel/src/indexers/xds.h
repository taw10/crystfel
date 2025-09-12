/*
 * xds.h
 *
 * Invoke xds for crystal autoindexing
 *
 * Copyright © 2013-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2013 Cornelius Gati
 *
 * Authors:
 *   2010-2017 Thomas White <taw@physics.org>
 *   2013      Cornelius Gati <cornelius.gati@cfel.de>
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

#ifndef XDS_H
#define XDS_H

#include "cell.h"
#include "index.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file xds.h
 * XDS indexer interface
 */

extern int run_xds(struct image *image, void *ipriv);

extern void *xds_prepare(IndexingMethod indm, UnitCell *cell);

extern const char *xds_probe(UnitCell *cell);

extern void xds_cleanup(void *pp);

#ifdef __cplusplus
}
#endif

#endif	/* XDS_H */
