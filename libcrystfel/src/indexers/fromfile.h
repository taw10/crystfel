/*
 * fromfile.h
 *
 * Perform indexing from solution file
 *
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

#ifndef FROMFILE_H
#define FROMFILE_H

#include <argp.h>

#include "image.h"

extern int fromfile_default_options(FromFileOptions **opts_ptr);
extern void *fromfile_prepare(IndexingMethod *indm,
                              struct fromfile_options *opts);
extern int fromfile_index(struct image *image, void *mpriv, int crystal_number);
extern void fromfile_cleanup(void *mpriv);

#endif	/* FROMFILE_H */
