/*
 * detgeom.h
 *
 * Detector geometry structure
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2009-2020 Thomas White <taw@physics.org>
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

#ifndef DETGEOM_H
#define DETGEOM_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file detgeom.h
 * Detector geometry structure and related functions.
 */

#include <gsl/gsl_matrix.h>

/**
 * Represents one panel of a detector
 */
struct detgeom_panel
{
	/** Text name for panel */
	char *name;

	/** \name Location of corner in units of the pixel size of this panel, \
	 *  measured from the interaction point. */
	/**@{*/
	double      cnx;
	double      cny;
	double      cnz;
	/**@}*/

	/** Pixel size in metres */
	double      pixel_pitch;

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

	/** \name Width and height of panel */
	/*@{*/
	int         w;
	int         h;
	/*@}*/

	/** \name Leaf group containing this panel (only) */
	const struct detgeom_panel_group *group;
};


struct detgeom_panel_group
{
	char *name;
	int n_children;

	struct detgeom_panel_group *parent;
	int serial;

	/* Center of panel group, in lab coordinate system (metres)
	 * This will be the rotation center. */
	double cx;
	double cy;
	double cz;

	/* If n_children > 0, here are the child groups */
	struct detgeom_panel_group **children;

	/* If n_children == 0, this is a leaf node, so: */
	struct detgeom_panel *panel;
};


struct detgeom
{
	struct detgeom_panel *panels;
	int n_panels;

	struct detgeom_panel_group *top_group;
};


extern void detgeom_transform_coords(struct detgeom_panel *p,
                                     double fs, double ss,
                                     double wavelength,
                                     double dx, double dy,
                                     double *r);

extern void detgeom_free(struct detgeom *detgeom);

extern double detgeom_max_resolution(struct detgeom *detgeom,
                                     double wavelength);

extern void show_panel(struct detgeom_panel *p);

extern double detgeom_mean_camera_length(struct detgeom *dg);

extern struct detgeom_panel *detgeom_find_panel(struct detgeom *dg, const char *name);

extern void detgeom_show_hierarchy(const struct detgeom *dg);

extern void detgeom_translate_detector_m(struct detgeom *dg, double x, double y, double z);

extern int detgeom_group_center(const struct detgeom_panel_group *grp,
                                double *x, double *y, double *z);

extern gsl_matrix **make_panel_minvs(struct detgeom *dg);

#ifdef __cplusplus
}
#endif

#endif	/* DETGEOM_H */
