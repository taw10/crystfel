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
	int      min_x;    /* Smallest x value considered to be in this panel */
	int      max_x;    /* Largest x value considered to be in this panel */
	int      min_y;    /* ... and so on */
	int      max_y;
	float    cx;       /* Location of centre */
	float    cy;
	float    clen;     /* Camera length */
	float    res;      /* Resolution */
	char     badrow;   /* 'x' or 'y' */
	int      no_index; /* Don't index peaks in this panel if non-zero */
};

struct detector
{
	struct panel *panels;
	int           n_panels;
};


/* x,y in pixels relative to central beam */
extern int map_position(struct image *image, double x, double y,
                        double *rx, double *ry, double *rz);

extern struct rvec get_q(struct image *image, unsigned int xs, unsigned int ys,
                         unsigned int sampling, float *ttp, float k);

extern double get_tt(struct image *image, unsigned int xs, unsigned int ys);

extern void record_image(struct image *image, int do_poisson);

extern struct panel *find_panel(struct detector *det, int x, int y);

extern struct detector *get_detector_geometry(const char *filename);

#endif	/* DETECTOR_H */
