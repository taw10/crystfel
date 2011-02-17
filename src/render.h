/*
 * render.h
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef RENDER_H
#define RENDER_H


#include <stddef.h>

#include "image.h"

enum {
	SCALE_COLOUR,
	SCALE_MONO,
	SCALE_INVMONO
};

/* Colour scale lookup */
extern void render_scale(float val, float max, int scale,
                         float *rp, float *gp, float *bp);


#ifdef HAVE_GTK

#include <gdk-pixbuf/gdk-pixbuf.h>

extern GdkPixbuf *render_get_image(struct image *image,
                                   int binning, int scale, double boost);
extern GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int scale);

extern int render_png(GdkPixbuf *pixbuf, const char *filename);
extern int render_tiff_fp(struct image *image, const char *filename);
extern int render_tiff_int16(struct image *image, const char *filename,
                             double boost);

#endif /* HAVE_GTK */



#endif	/* RENDER_H */
