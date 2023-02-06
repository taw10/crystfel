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

    /** Bias applied to the raw detector intensity unit (or dark-value), in  ADU */
	double      adu_bias;

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
};


struct detgeom
{
	struct detgeom_panel *panels;
	int                   n_panels;
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

#ifdef __cplusplus
}
#endif

#endif	/* DETGEOM_H */
