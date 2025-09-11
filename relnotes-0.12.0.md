Release notes for CrystFEL version 0.12.0
=========================================

Copyright © 2012-2025 Deutsches Elektronen-Synchrotron DESY,
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

With this release, CrystFEL now enforces the geometrical restrictions of the
lattice type.  There are also further improvements to detector geometry
refinement, and two new indexing algorithms.  These changes are detailed below.
In addition, there were many smaller fixes and improvements.  See the ChangeLog
or the Git history for a comprehensive list of all changes.


Enforcement of lattice type restrictions
----------------------------------------

The geometrical restrictions of the lattice type are now enforced when indexing
diffraction patterns.  For instance, if the target unit cell specifies a cubic
lattice, then indexamajig will output cells with all angles 90° and all axis
lengths equal.

Enforcement of these restrictions gives better results in almost all cases.
Nevertheless, if you need to disable the enforcement, specify the lattice type
as triclinic (instead of whatever the real lattice type is).  The axis lengths
and angles of indexing solutions will still be checked for correspondence
(within tolerances) to the target unit cell, as before.

The restrictions are applied at the refinement stage, so will not apply if you
disable the refinement with option --no-refine.


Improvements/bug fixes for detector geometry refinement
-------------------------------------------------------

Detector geometry refinement using align_detector/Millepede benefits
substantially from the enforcement of lattice type restrictions described
above.  In addition, this release fixes a bug which was leading to inaccurate
results for out-of-plane panel tilts.  Refining out-of-plane tilts and shifts,
including the overall camera length, now works reliably provided that the
lattice has monoclinic or higher symmetry.

The low-level multiprocessing code for indexamajig was overhauled in this
release.  One of the benefits of this is that only one "Mille" data file will
now be written for all worker processes in a single indexamajig run.

A new program, show_residuals, has been added to calculate the average spot
position residuals and therefore give a figure of merit for the detector
geometry.

When writing updated geometry files, align_detector and adjust_detector will
now preserve comments placed at the top of the input geometry file.  This can
be useful for keeping a record of the geometry file's history.

The maximum number of panel groups has also been increased from 256 to 512,
which fixes a problem with detectors with very large numbers of panels (e.g.
LPD-1M).


New indexing algorithms and documentation
-----------------------------------------

This release adds two new indexing algorithms, ffbidx and smallcell.

ffbidx (Fast Feedback Indexer) is a GPU-based indexing algorithm for rapid data
evaluation, and has been implemented in CrystFEL by Hans-Christian Stadler.
It's a derivative of the PyTorch-based "TORO" indexer, using CUDA instead of
PyTorch.  See https://doi.org/10.1107/S1600576724003182 for details.

smallcell is an implementation of the graph-theoretical "cctbx.small_cell"
algorithm for indexing diffraction patterns with small numbers of reflections.
The CrystFEL implementation was written by Isabel Costello.  For details about
the algorithm, see https://doi.org/10.1107/S1399004714026145.

To help you navigate amongst the many available indexing algorithms, a new
document has been added in doc/articles/indexer-choice.rst
(https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/indexer-choice.rst)

This CrystFEL release also fixes a bug when using asdf to index rhombohedral
lattices on hexagonal axes ("H-centered" cells).


API changes
-----------

Removed routines:
* `data_template_get_2d_detgeom_if_possible`
* `stream_get_fd`

Added routines:
* `powder_rings`
* `impose_bravais`
* `cell_rotate_gsl_direct`
* `data_template_get_detgeom_if_possible`
* `data_template_reset_total_movements`
* `data_template_print_total_movements`
* `crystfel_mille_new_fd`
* `image_read_header_int`
* `convert_long_int`
* `stream_get_fh`
* `crossp_norm`,
* `rotate3d`
* `adjust_vector_length`
* `print_indexing_info`
* `index_pattern_5` (without `ping` and `last_task`)

Added type definitions:
* `struct powder_ring`

Changed type definitions:
* `dg_group_info` - addition of `leaf`
* `header_cache_entry` - `val_int` is now `long long int`
* `image` - removal of `id`
* `gparam` - `GPARAM_ASX` etc removed, `GPARAM_A_STAR` etc added

Changed routine prototypes:
* `write_mille` - addition of `rvl` and `nl`
* `image_cache_header_int` - `header_val` is now `long long int`
* `refine_prediction` - addition of `target`
* `write_to_mtz` - addition of `spg_input`
* `default_method_options` and `setup_indexing` - addition of `ffbidx_opts` and `smallcell_opts`
