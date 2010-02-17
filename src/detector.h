/*
 * detector.h
 *
 * Detector properties
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DETECTOR_H
#define DETECTOR_H

struct image;

#include "image.h"

struct panel
{
	int      min_x;   /* Smallest x value considered to be in this panel */
	int      max_x;   /* Largest x value considered to be in this panel */
	int      min_y;   /* ... and so on */
	int      max_y;
	float    cx;      /* Location of centre */
	float    cy;
	float    clen;    /* Camera length */
	float    res;     /* Resolution */
};

struct detector
{
	struct panel *panels;
	int           n_panels;
};

extern void record_image(struct image *image, int do_water, int do_poisson,
                         int do_bloom);

#endif	/* DETECTOR_H */
