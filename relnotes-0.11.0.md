Release notes for CrystFEL version 0.11.0
=========================================

Copyright © 2012-2024 Deutsches Elektronen-Synchrotron DESY,
                      a research centre of the Helmholtz Association.

See AUTHORS as well as the individual source code files for full contributor details.

CrystFEL is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CrystFEL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CrystFEL.  If not, see <http://www.gnu.org/licenses/>.


Overview
--------

This CrystFEL release adds a totally new system for detector geometry
refinement, as well as better support for compiling on Mac OS.  The old CMake
build system has been removed in favour of Meson.  The simulation tools
`pattern_sim` and `partial_sim`, which were not well maintained nor widely used,
have been removed.

These changes are detailed below.  In addition, there were many smaller fixes
and improvements.  See the ChangeLog or the Git history for a comprehensive
list of all changes.


New system for detector geometry refinement
-------------------------------------------

The new detector geometry system is based on the Millepede-II program, the same
software used for high energy physics experiments including the CMS.  It
allows us to do a full refinement including all of the inter-dependencies
between crystal orientations, cell parameters and the detector geometry.
Previously, the refined updates to the detector positions would be biased
because they were based on biased crystal parameters, which were in turn biased
because of the incorrect detector geometry.  To break this circular dependency,
we have to refine the individual parameters for each crystal at the same time
as we refine the geometry.  Normally, doing this would result in an enormous
and computationally prohibitive least squares calculation, even with
specialised sparse matrix solvers.  However, the Millepede algorithm reduces
the calculation to one with a sized based only on the number of geometrical
parameters for the detector.

To get started with the new geometry refinement, read the manual for the new
`align_detector` program by running `man align_detector` or visiting
<https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/man/align_detector.1.md>.
An update of the geometry file will be necessary to add information about the
hierarchy of detector panels, but example files are available.  The refinement
can also be performed using the CrystFEL GUI.

A side-effect of this change is that the beam center is no longer refined
after indexing for each pattern individually.  This should make the prediction
refinement more stable.

Future CrystFEL versions will look at improving the stability and precision of
this method, in particular for three-dimensional and fine-grained refinement
tasks.


Homebrew formula and improved Mac OS support
--------------------------------------------

This release adds a Homebrew formula for CrystFEL.  If you're running on a Mac,
you can now install CrystFEL with a single `brew install` command.  See
`INSTALL.md` for more details.  Eventually we hope to submit this formula to
the main Homebrew repository to make it even easier to install.

In addition, we are now regularly testing CrystFEL on Mac OS as part of our
continuous integration pipeline.


Removal of simulation tools and CMake build system
--------------------------------------------------

Due mostly to a insufficient development resources, the simulation tools
`pattern_sim` and `partial_sim` have been removed from the CrystFEL suite
starting from this version.  For more detailed discussion of the rationale for
this, see https://gitlab.desy.de/thomas.white/crystfel/-/issues/81

If you need quantitative image simulations, there are several better options
including [nanoBragg]<https://bl831.als.lbl.gov/~jamesh/nanoBragg/>.  For
relative-scale geometry-only simulations of the kind done by `partial_sim`,
a better option is expected in a future CrystFEL version - watch this space!

The CMake-based build system, which was barely maintained and lacked many
features compared to the Meson-based system, has now been removed altogether.


API changes
-----------

Removed routines:
* `data_template_get_rigid_groups`
* `r_gradient`
* `x_gradient`
* `y_gradient`
* `image_hdf5_write`

Added routines:
* `crystfel_mille_new`
* `crystfel_mille_free`
* `crystfel_mille_delete_last_record`
* `crystfel_mille_write_record`
* `mille_label`
* `mille_unlabel`
* `write_mille`
* `data_template_show_hierarchy`
* `data_template_translate_group_px`
* `data_template_translate_group_m`
* `data_template_rotate_group`
* `data_template_write_to_fh`
* `data_template_write_to_file`
* `data_template_group_info`
* `detgeom_find_panel`
* `detgeom_show_hierarchy`
* `detgeom_translate_detector_m`
* `detgeom_group_center`
* `make_panel_minvs`
* `is_cbf_file`
* `is_cbfgz_file`
* `is_hdf5_file`
* `image_data_arrays_new`
* `image_data_arrays_free`
* `image_create_dp_bad`
* `image_set_zero_data`
* `index_pattern_4`
* `r_dev`
* `fs_dev`
* `ss_dev`
* `r_gradient`
* `fs_ss_gradient`
* `read_reflections_3`
* `stream_write_data_template`
* `show_vector`
* `matrix_mult2`
* `matrix_mult3`
* `matrix_invert`
* `sq`
* `rotate2d`
* `is_dir`

Added type definitions:
* `Mille`
* `struct dg_group_info`
* `struct detgeom_panel_group`
* `ImageDataArrays`
* `struct peak_params`
* `struct reflpeak`

Removed type definitions:
* `struct rigid_group`
* `struct rg_collection`

Changed structure definitions:
* `detgeom_panel` - addition of `group`
* `detgeom` - addition of `top_group`
* `image` - addition of `ida`
* `enum gparam` - changed values, moved from geometry.h to predict-refine.h

Changed routine prototypes:
* `image_read` - addition of `ida`
* `image_read_data_block` - addition of `ida`
* `refine_prediction` - addition of `mille` and `max_mille_level`
