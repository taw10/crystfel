/*
 * render.h
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * (c) 2008 Thomas White <taw27@cam.ac.uk>
 *
 *  micronview - view Ditabis Micron images
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef RENDER_H
#define RENDER_H


#include <gdk-pixbuf/gdk-pixbuf.h>
#include <stddef.h>

extern GdkPixbuf *render_get_image(struct hdfile *micron, int binning,
                                   int boostint, int monochrome);
extern GdkPixbuf *render_get_colour_scale(size_t w, size_t h, int monochrome);


#endif	/* RENDER_H */
