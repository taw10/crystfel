/*
 * image-seedee.h
 *
 * Image loading, SeeDee parts
 *
 * Copyright © 2012-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2022 Thomas White <taw@physics.org>
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

#ifndef IMAGE_SEEDEE_H
#define IMAGE_SEEDEE_H

#include "datatemplate.h"


extern int image_seedee_read(struct image *image,
                             const DataTemplate *dtempl,
                             void *data,
                             size_t data_size,
                             char *meta_data);

#endif	/* IMAGE_SEEDEE_H */
