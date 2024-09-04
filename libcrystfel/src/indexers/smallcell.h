/*
 * smallcell.h
 *
 * Perform indexing from solution file
 *
 * Copyright © 2020-2021 Max-Planck-Gesellschaft
 *                       zur Förderung der Wissenschaften e.V.
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Pascal Hogan-Lamarre <pascal.hogan.lamarre@mail.utoronto.ca>
 *   2021 Thomas White <thomas.white@desy.de>
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

#include <argp.h>

#include "image.h"

extern int smallcell_default_options(struct smallcell_options **opts_ptr);
extern void *smallcell_prepare(IndexingMethod *indm,
                              struct smallcell_options *opts, UnitCell *cell);
extern int smallcell_index(struct image *image, void *mpriv);
extern void smallcell_cleanup(void *mpriv);

#endif	/* SMALLCELL_H */
