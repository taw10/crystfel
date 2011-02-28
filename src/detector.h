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
struct hdfile;

#include "hdf5-file.h"
#include "image.h"


struct panel
{
	int      min_fs;  /* Smallest FS value considered to be in the panel */
	int      max_fs;  /* Largest FS value considered to be in this panel */
	int      min_ss;  /* ... and so on */
	int      max_ss;
	double   cx;       /* Location of corner (min_fs,min_ss) in pixels */
	double   cy;
	double   clen;     /* Camera length in metres */
	char    *clen_from;
	double   res;      /* Resolution in pixels per metre */
	char     badrow;   /* 'x' or 'y' */
	int      no_index; /* Don't index peaks in this panel if non-zero */
	double   peak_sep; /* Characteristic peak separation */

	signed int fsx;
	signed int fsy;
	signed int ssx;
	signed int ssy;
};

struct detector
{
	struct panel *panels;
	int           n_panels;
	int           max_fs;
	int           max_ss;  /* Size of overall array needed, minus 1 */
};

extern struct rvec get_q(struct image *image, double xs, double ys,
                         unsigned int sampling, float *ttp, float k);

extern double get_tt(struct image *image, double xs, double ys);

extern void record_image(struct image *image, int do_poisson);

extern struct panel *find_panel(struct detector *det, int x, int y);

extern struct detector *get_detector_geometry(const char *filename);

extern void free_detector_geometry(struct detector *det);

extern struct detector *simple_geometry(const struct image *image);

extern void get_pixel_extents(struct detector *det,
                              double *min_x, double *min_y,
                              double *max_x, double *max_y);

extern void fill_in_values(struct detector *det, struct hdfile *f);


#endif	/* DETECTOR_H */
