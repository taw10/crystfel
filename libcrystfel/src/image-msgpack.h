/*
 * image-cbf.h
 *
 * Image loading, MessagePack parts
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMAGE_MSGPACK_H
#define IMAGE_MSGPACK_H

#include "datatemplate.h"

#if defined(HAVE_MSGPACK)

#include <msgpack.h>

extern struct image *image_msgpack_read(DataTemplate *dtempl,
                                        msgpack_object *obj,
                                        int no_image_data);

extern ImageFeatureList *image_msgpack_read_peaks(const DataTemplate *dtempl,
                                                  msgpack_object *obj,
                                                  int half_pixel_shift);

#else /* defined(HAVE_MSGPACK) */

static UNUSED struct image *image_msgpack_read(DataTemplate *dtempl,
                                                void *obj,
                                                int no_image_data)
{
	return NULL;
}

static UNUSED ImageFeatureList *image_msgpack_read_peaks(const DataTemplate *dtempl,
                                                         void *obj,
                                                         int half_pixel_shift)
{
	return NULL;
}

#endif /* defined(HAVE_MSGPACK) */

#endif	/* IMAGE_MSGPACK_H */
