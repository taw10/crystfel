/*
 * hdfsee-render.h
 *
 * Rendering bits for hdfsee
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2011-2015 Thomas White <taw@physics.org>
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

#ifndef HDFSEE_RENDER_H
#define HDFSEE_RENDER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_GDKPIXBUF
#include <gdk-pixbuf/gdk-pixbuf.h>

extern GdkPixbuf **render_panels(struct image *image,
                                 int binning, int scale, double boost,
                                 int *n_pixbufs);

extern GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int scale);

#endif /* HAVE_GDKPIXBUF */

extern int render_tiff_fp(struct image *image, const char *filename, int min_x,
                          int max_x, int min_y, int max_y);

extern int render_tiff_int16(struct image *image, const char *filename, double boost,
                             int min_x, int max_x, int min_y, int max_y);

#endif	/* HDFSEE_RENDER_H */
