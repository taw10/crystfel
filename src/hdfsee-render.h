/*
 * hdfsee-render.h
 *
 * Rendering bits for hdfsee
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
