/*
 * datatemplate_priv.h
 *
 * Data template structure (private parts)
 *
 * Copyright © 2019-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019-2020 Thomas White <taw@physics.org>
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


enum adu_per_unit
{
	ADU_PER_PHOTON,
	ADU_PER_EV
};


struct dt_rigid_group
{
	char *name;
	struct panel_template **panels;
	int n_panels;
};


struct dt_rg_collection
{
	char *name;
	struct dt_rigid_group **rigid_groups;
	int n_rigid_groups;
};


/**
 * Represents one panel of a detector
 */
struct panel_template
{
	/** Text name for panel */
	char *name;

	/** \name Location of corner in units of the pixel size of this panel */
	/**@{*/
	double cnx;
	double cny;
	/**@}*/

	/** Location to get \ref cnz from, e.g. from HDF5 file */
	char *cnz_from;

	/** The offset to be applied from \ref clen */
	double cnz_offset;

	/** Location of mask data */
	char *mask;

	/** Filename for mask data */
	char *mask_file;

	/** Location of per-pixel saturation map */
	char *satmap;

	/** Filename for saturation map */
	char *satmap_file;

	/** Mark entire panel as bad if set */
	int bad;

	/** Resolution in pixels per metre  */
	double pixel_pitch;

	/** Number of detector intensity units per photon, or eV */
	double adu_scale;
	enum adu_per_unit adu_scale_unit;

	/** Treat pixel as unreliable if higher than this */
	double max_adu;

	/** Location of data in file */
	char *data;

	/** Dimension structure */
	struct dim_structure *dim_structure;

	/** \name Transformation matrix from pixel coordinates to lab frame */
	/*@{*/
	double fsx;
	double fsy;
	double fsz;
	double ssx;
	double ssy;
	double ssz;
	/*@}*/

	/** \name Rail direction */
	/*@{*/
	double rail_x;
	double rail_y;
	double rail_z;
	/*@}*/

	/* Value of clen (without coffset) at which beam is centered */
	double clen_for_centering;

	/** \name Position of the panel in the data block in the file. */
	/*@{*/
	int orig_min_fs;
	int orig_max_fs;
	int orig_min_ss;
	int orig_max_ss;
	/*@}*/


	/* FIXME: Should only be one */
	double                     adu_per_eV;
	double                     adu_per_photon;

};


struct dt_badregion
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


struct _datatemplate
{
	struct panel_template     *panels;
	int                        n_panels;

	struct dt_badregion       *bad;
	int                        n_bad;

	char                      *wavelength_from;
	double                     photon_energy_bandwidth;  /* Eww */
	double                     photon_energy_scale;      /* Eww */

	char                      *peak_info_location;

	unsigned int               mask_bad;
	unsigned int               mask_good;

	struct dt_rigid_group    **rigid_groups;
	int                        n_rigid_groups;

	struct dt_rg_collection  **rigid_group_collections;
	int                        n_rg_collections;

	int                        path_dim;
	int                        dim_dim;

	char                      *hdf5_peak_path;

	struct panel_template      defaults;
};
