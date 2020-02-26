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

#ifdef __cplusplus
}
#endif

#endif	/* DETGEOM_H */
