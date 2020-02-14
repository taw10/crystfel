/*
 * detgeom.h
 *
 * Detector geometry structure
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2019 Thomas White <taw@physics.org>
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

#ifndef DETGEOM_H
#define DETGEOM_H

struct rigid_group;
struct rg_collection;
struct detector;
struct panel;
struct badregion;

#include "hdf5-file.h"
#include "image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file detgeom.h
 * Detector geometry structure and related functions.
 */


struct rigid_group
{
	char *name;
	struct panel **panels;
	int n_panels;
};


struct rg_collection
{
	char *name;
	struct rigid_group **rigid_groups;
	int n_rigid_groups;
};


/**
 * Represents one panel of a detector
 */
struct panel
{
	/** Text name for panel */
	const char *name;

	/** \name Location of corner in units of the pixel size of this panel, \
	 *  measured from the interaction point. */
	/**@{*/
	double      cnx;
	double      cny;
	double      cnz;
	/**@}*/

	/** Pixel size in metres */
	double      pixel_pitch;

	/** Readout direction (for filtering out clusters of peaks)
	  * ('x' or 'y') */
	char        badrow;

	/** Number of detector intensity units per photon (or electron, etc) */
	double      adu_per_photon;

	/** Treat pixel as unreliable if higher than this */
	double      max_adu;

	/** \name Transformation matrix from pixel coordinates to lab frame */
	/*@{*/
	double      fsx;
	double      fsy;
	double      fsz;
	double      ssx;
	double      ssy;
	double      ssz;
	/*@}*/

	/** \name Inverse of 2D part of transformation matrix */
	/*@{*/
	double      xfs;
	double      yfs;
	double      xss;
	double      yss;
	/*@}*/

	/** \name Width and height of panel */
	/*@{*/
	int         w;
	int         h;
	/*@}*/
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

	struct rigid_group **rigid_groups;
	int                  n_rigid_groups;

	struct rg_collection **rigid_group_collections;
	int                    n_rg_collections;

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
};


extern struct rvec get_q_for_panel(struct panel *p, double fs, double ss,
                                   double *ttp, double k);

extern double get_tt(struct image *image, double xs, double ys, int *err);

extern struct panel *find_orig_panel(struct detector *det,
                                     double fs, double ss);

extern signed int find_orig_panel_number(struct detector *det,
                                         double fs, double ss);

extern int panel_number(const struct detector *det, const struct panel *p);

extern struct detector *get_detector_geometry(const char *filename,
                                              struct beam_params *beam);

extern struct detector *get_detector_geometry_2(const char *filename,
                                                struct beam_params *beam,
                                                char **hdf5_peak_path);

extern struct detector *get_detector_geometry_from_string(const char *string,
                                                          struct beam_params *beam,
                                                          char **hdf5_peak_path);

extern void free_detector_geometry(struct detector *det);

extern void get_pixel_extents(struct detector *det,
                              double *min_x, double *min_y,
                              double *max_x, double *max_y);

extern int panel_is_in_rigid_group(const struct rigid_group *rg,
                                   struct panel *p);

extern int rigid_group_is_in_collection(struct rg_collection *c,
                                        struct rigid_group *rg);

extern int reverse_2d_mapping(double x, double y, struct detector *det,
                              struct panel **pp, double *pfs, double *pss);

extern double largest_q(struct image *image);

extern double smallest_q(struct image *image);

extern struct panel *find_panel_by_name(struct detector *det, const char *name);

extern int write_detector_geometry_2(const char *geometry_filename,
                                     const char *output_filename,
                                     struct detector *det,
                                     const char *additional_comment,
                                     int write_panel_coffset);

extern int write_detector_geometry_3(const char *geometry_data,
                                     const char *output_filename,
                                     struct detector *det,
                                     const char *additional_comment,
                                     int write_panel_coffset);

extern int write_detector_geometry(const char *geometry_filename,
                                   const char *output_filename,
                                   struct detector *det);

extern void mark_resolution_range_as_bad(struct image *image,
                                         double min, double max);

struct rg_collection *find_rigid_group_collection_by_name(struct detector *det,
                                                          const char *name);

#ifdef __cplusplus
}
#endif

#endif	/* DETGEOM_H */
