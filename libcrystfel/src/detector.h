/*
 * detector.h
 *
 * Detector properties
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *   2011-2012 Richard Kirian <rkirian@asu.edu>
 *   2011      Andrew Aquila
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DETECTOR_H
#define DETECTOR_H

struct image;
struct hdfile;

#include "hdf5-file.h"
#include "image.h"

#ifdef __cplusplus
extern "C" {
#endif

struct panel
{
	char     name[1024];  /* Name for this panel */

	int      min_fs;  /* Smallest FS value considered to be in the panel */
	int      max_fs;  /* Largest FS value considered to be in this panel */
	int      min_ss;  /* ... and so on */
	int      max_ss;
	double   cnx;       /* Location of corner (min_fs,min_ss) in pixels */
	double   cny;
	double   coffset;
	double   clen;     /* Camera length in metres */
	char    *clen_from;
	double   res;      /* Resolution in pixels per metre */
	char     badrow;   /* 'x' or 'y' */
	int      no_index; /* Don't index peaks in this panel if non-zero */
	char    *rigid_group;  /* Rigid group, or -1 for none */
	double   adu_per_eV;   /* Number of ADU per eV */
	double   max_adu;  /* Treat pixel as unreliable if higher than this */

	double fsx;
	double fsy;
	double ssx;
	double ssy;

	double xfs;
	double yfs;
	double xss;
	double yss;
};


struct badregion
{
	char name[1024];
	double min_x;
	double max_x;
	double min_y;
	double max_y;
};


struct detector
{
	struct panel     *panels;
	int               n_panels;

	int               max_fs;
	int               max_ss;  /* Size of overall array needed, minus 1 */

	struct badregion *bad;
	int               n_bad;

	char              *mask;
	unsigned int       mask_bad;
	unsigned int       mask_good;

	char             **rigid_groups;
	int                num_rigid_groups;

	struct panel       defaults;
};


extern struct rvec get_q(struct image *image, double fs, double ss,
                         double *ttp, double k);

extern struct rvec get_q_for_panel(struct panel *p, double fs, double ss,
                                   double *ttp, double k);

extern double get_tt(struct image *image, double xs, double ys, int *err);

extern int in_bad_region(struct detector *det, double fs, double ss);

extern void record_image(struct image *image, int do_poisson);

extern struct panel *find_panel(struct detector *det, double fs, double ss);

extern int find_panel_number(struct detector *det, int fs, int ss);

extern struct detector *get_detector_geometry(const char *filename);

extern void free_detector_geometry(struct detector *det);

extern struct detector *simple_geometry(const struct image *image);

extern void get_pixel_extents(struct detector *det,
                              double *min_x, double *min_y,
                              double *max_x, double *max_y);

extern void fill_in_values(struct detector *det, struct hdfile *f);

extern struct detector *copy_geom(const struct detector *in);

extern int reverse_2d_mapping(double x, double y, double *pfs, double *pss,
                              struct detector *det);

extern double largest_q(struct image *image);

extern double smallest_q(struct image *image);

extern struct panel *find_panel_by_name(struct detector *det, const char *name);

#ifdef __cplusplus
}
#endif

#endif	/* DETECTOR_H */
