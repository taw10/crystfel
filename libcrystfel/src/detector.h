/*
 * detector.h
 *
 * Detector properties
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2015 Thomas White <taw@physics.org>
 *   2011-2012 Richard Kirian <rkirian@asu.edu>
 *   2014      Valerio Mariani
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

struct rigid_group;
struct rg_collection;
struct detector;
struct panel;
struct badregion;
struct detector;
struct beam_params;
struct hdfile;
struct event;

#include "hdf5-file.h"
#include "image.h"

#ifdef __cplusplus
extern "C" {
#endif


struct rigid_group
{
	char *name;
	struct panel **panels;
	int n_panels;

	/* Updates to panel position calculated during integration */
	double d_fsx;
	double d_ssx;
	double d_cnx;
	double d_fsy;
	double d_ssy;
	double d_cny;
	int have_deltas;
};


struct rg_collection
{
	char *name;
	struct rigid_group **rigid_groups;
	int n_rigid_groups;
};


struct panel
{
        char     name[1024];  /* Name for this panel */

        /* Position of panel in the data block in memory (see below) */
        int      min_fs;  /* Smallest FS value considered to be in the panel */
        int      max_fs;  /* Largest FS value considered to be in this panel */
        int      min_ss;  /* ... and so on */
        int      max_ss;

        double   cnx;       /* Location of corner (min_fs,min_ss) in pixels */
        double   cny;
        double   coffset;
        double   clen;     /* Camera length in metres */
        char    *clen_from;
        char    *mask;
        double   res;      /* Resolution in pixels per metre */
        char     badrow;   /* 'x' or 'y' */
        int      no_index; /* Don't index peaks in this panel if non-zero */
        double   adu_per_eV;   /* Number of ADU per eV */
        double   max_adu;  /* Treat pixel as unreliable if higher than this */
        char    *data;

        struct dim_structure *dim_structure;

        double fsx;
        double fsy;
        double ssx;
        double ssy;

        double xfs;
        double yfs;
        double xss;
        double yss;

        /* Position of the panel in the data block in the file.  The panels may
         * get moved around when the file is loaded (see hdf5_read2()),
         * especially if the panels come from different HDF5 elements. */
        int orig_min_fs;
        int orig_max_fs;
        int orig_min_ss;
        int orig_max_ss;

        int w;  /* Width, calculated as max_fs-min_fs+1 */
        int h;  /* Height, calculated as max_ss-min_ss+1 */
};


struct badregion
{
	char name[1024];
	int is_fsss;
	char *panel;

	double min_x;
	double max_x;
	double min_y;
	double max_y;

	/* Specified INCLUSIVELY */
	int      min_fs;
	int      max_fs;
	int      min_ss;
	int      max_ss;

};


struct detector
{
	struct panel     *panels;
	int               n_panels;

	int               max_fs;
	int               max_ss;  /* Size of overall array needed, minus 1 */

	struct badregion *bad;
	int               n_bad;

	unsigned int       mask_bad;
	unsigned int       mask_good;

	struct rigid_group **rigid_groups;
	int                n_rigid_groups;

	struct rg_collection **rigid_group_collections;
	int               	 n_rg_collections;

	/* Location of the pixel furthest away from the beam position, which
	 * will have the largest value of 2theta regardless of camera length
	 * and wavelength */
	struct panel      *furthest_out_panel;
	double             furthest_out_fs;
	double             furthest_out_ss;

	/* As above, but for the smallest 2theta */
	struct panel      *furthest_in_panel;
	double             furthest_in_fs;
	double             furthest_in_ss;

	int                path_dim;
	int                dim_dim;

	struct panel       defaults;
};


extern struct rvec get_q(struct image *image, double fs, double ss,
                         double *ttp, double k);

extern struct rvec get_q_for_panel(struct panel *p, double fs, double ss,
                                   double *ttp, double k);

extern double get_tt(struct image *image, double xs, double ys, int *err);

extern int in_bad_region(struct detector *det, double fs, double ss);

extern void record_image(struct image *image, int do_poisson, double background,
                         gsl_rng *rng, double beam_radius, double nphotons);

extern struct panel *find_panel(struct detector *det, double fs, double ss);
extern signed int find_panel_number(struct detector *det, double fs, double ss);
extern struct panel *find_orig_panel(struct detector *det,
                                     double fs, double ss);

extern struct detector *get_detector_geometry(const char *filename,
                                              struct beam_params *beam);

extern void free_detector_geometry(struct detector *det);

extern struct detector *simple_geometry(const struct image *image);

extern void get_pixel_extents(struct detector *det,
                              double *min_x, double *min_y,
                              double *max_x, double *max_y);

extern void fill_in_values(struct detector *det, struct hdfile *f,
                           struct event* ev);

extern int panel_is_in_rigid_group(const struct rigid_group *rg,
                                   struct panel *p);

extern int rigid_group_is_in_collection(struct rg_collection *c,
                                        struct rigid_group *rg);

extern struct detector *copy_geom(const struct detector *in);

extern int reverse_2d_mapping(double x, double y, double *pfs, double *pss,
                              struct detector *det);

extern double largest_q(struct image *image);

extern double smallest_q(struct image *image);

extern struct panel *find_panel_by_name(struct detector *det, const char *name);

extern int write_detector_geometry_2(const char *geometry_filename,
                                     const char *output_filename,
                                     struct detector *det,
                                     const char *additional_comment,
                                     int write_panel_coffset);

extern int write_detector_geometry(const char *geometry_filename,
                                   const char *output_filename,
                                   struct detector *det);

extern void mark_resolution_range_as_bad(struct image *image,
                                         double min, double max);


extern int single_panel_data_source (struct detector *det, const char *element);

struct rg_collection *find_rigid_group_collection_by_name(struct detector *det,
                                                          const char *name);

#ifdef __cplusplus
}
#endif

#endif	/* DETECTOR_H */
