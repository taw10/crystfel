/*
 * hdfsee-render.h
 *
 * Rendering bits for hdfsee
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
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


#ifdef HAVE_GTK

#include <gdk-pixbuf/gdk-pixbuf.h>

extern GdkPixbuf **render_panels(struct image *image,
                                 int binning, int scale, double boost,
                                 int *n_pixbufs);

extern GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int scale);

extern int render_tiff_fp(struct image *image, const char *filename);
extern int render_tiff_int16(struct image *image, const char *filename,
                             double boost);

#endif /* HAVE_GTK */


#endif	/* HDFSEE_RENDER_H */
