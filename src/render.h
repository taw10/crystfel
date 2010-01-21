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


#include <gdk-pixbuf/gdk-pixbuf.h>
#include <stddef.h>

#include "displaywindow.h"

extern GdkPixbuf *render_get_image(DisplayWindow *dw);
extern GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int monochrome);


#endif	/* RENDER_H */
