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
	float    cx;       /* Location of centre in pixels */
	float    cy;
	float    clen;     /* Camera length in metres */
	float    res;      /* Resolution in pixels per metre */
	char     badrow;   /* 'x' or 'y' */
	int      no_index; /* Don't index peaks in this panel if non-zero */
	float    peak_sep; /* Characteristic peak separation */

	signed int fsx;
	signed int fsy;
	signed int ssx;
	signed int ssy;
};

struct detector
{
	struct panel *panels;
	int           n_panels;
	int           max_x;
	int           max_y;  /* Size of overall array needed, minus 1 */
};

extern struct rvec get_q(struct image *image, double xs, double ys,
                         unsigned int sampling, float *ttp, float k);

extern double get_tt(struct image *image, double xs, double ys);

extern void record_image(struct image *image, int do_poisson);

extern struct panel *find_panel(struct detector *det, int x, int y);

extern struct detector *get_detector_geometry(const char *filename);
extern void free_detector_geometry(struct detector *det);

#endif	/* DETECTOR_H */
