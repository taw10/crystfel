/*
 * smallcell.h
 *
 * Re-implementation of graph theory indexing algorithm for small unit cells
 *  borrowed from cctbx.small_cell
 *
 * Copyright © 2024 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2024 Isabel Costello <isabel.costello@desy.de>
 *   2024 Thomas White <thomas.white@desy.de>
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

#ifndef SMALLCELL_H
#define SMALLCELL_H

#include "image.h"
#include "cell.h"
#include "index.h"

extern int smallcell_default_options(struct smallcell_options **opts_ptr);
extern void *smallcell_prepare(IndexingMethod indm,
                               struct smallcell_options *opts, UnitCell *cell);
extern int smallcell_index(struct image *image, void *mpriv);
extern void smallcell_cleanup(void *mpriv);

#endif	/* SMALLCELL_H */
