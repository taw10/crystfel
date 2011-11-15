/*
 * render.h
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef RENDER_H
#define RENDER_H


enum {
	SCALE_COLOUR,
	SCALE_MONO,
	SCALE_INVMONO,
	SCALE_RATIO
};

/* Colour scale lookup */
extern void render_scale(double val, double max, int scale,
                         double *rp, double *gp, double *bp);


#endif	/* RENDER_H */
