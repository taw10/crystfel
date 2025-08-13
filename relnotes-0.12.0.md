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

This release includes a significant overhaul of the indexamajig sandbox, as well
as adding new features.

These changes are detailed below.  In addition, there were many smaller fixes
and improvements.  See the ChangeLog or the Git history for a comprehensive
list of all changes.

A screencast to accompany this release is available at the following locations:
<https://desy.de/~twhite/crystfel/presentations.html>
<https://vimeo.com/1017833941>


Enforcement of lattice type restrictions
----------------------------------------

Equal angles, lengths, 90° angles.

Improves detector geometry refinement a lot.



Improvements/bug fixes for detector geometry refinement
-------------------------------------------------------

Significant improvement due to lattice type.

Fixes for out-of-plane tilts.  Can usually refine camera length.

One Mille file for all workers.

New show_residuals program.

Preservation of comments.

Increased number of panel groups.



New indexing algorithms and documentation
-----------------------------------------

ffbidx

smallcell

New article about choice of algorithms

H-centering fix for asdf



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
