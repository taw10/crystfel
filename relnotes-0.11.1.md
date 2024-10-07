Release notes for CrystFEL version 0.11.1
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

This release is mainly a bug-fixing release to 0.11.0.

These changes are detailed below.  In addition, there were many smaller fixes
and improvements.  See the ChangeLog or the Git history for a comprehensive
list of all changes.

A screencast to accompany this release is available at the following locations:
<https://desy.de/~twhite/crystfel/presentations.html>
<https://vimeo.com/1017833941>


Improvements/bug fixes for detector geometry refinement
-------------------------------------------------------

Improved detector geometry refinement, based on the "Millepede" algorithm
borrowed from particle physics, was introduced in the previous CrystFEL version
(0.11.0).  Version 0.11.1 fixes several problems in the geometry files written
out by align_detector (and also adjust_detector), which prevented the refined
geometry files from being used in some circumstances.  Corrupted "Mille" data
was also being written by indexamajig in some cases, which has been fixed.
A further problem with the graphical user interface (GUI), where the geometry
file was not correctly found by the geometry refinement task, was corrected.

Two new options have been added to align_detector.  `--out-of-plane-tilts`
allows panels to rotate outside the plane of the detector, and
`--camera-length` enables refinement of the overall camera length.
Previously, both of these motions were enabled with `--out-of-plane`.  In the
new version. `--out-of-plane` only enables translations out of the detector
plane, not rotations, and does not allow the overall camera length to vary.

If no panel groups are defined in the geometry file at all, CrystFEL will now
automatically create a single top-level group.  This allows for basic (level 0)
geometry refinement without any extra work.

As a convenience, CrystFEL now accepts panel group definitions appearing near
the start of the geometry file, before the panels have been defined.  The
groups must still appear in correct hierarchy order (lowest levels first).
To maintain compatability with older CrystFEL versions, the panel groups will
be written at the end of the geometry file by align/adjust_detector.


GUI improvements
----------------

A new option "Request exclusive use of compute nodes" has been added to the
Slurm interface of the CrystFEL graphical user interface.  This option should
allow for more efficient cluster usage.

In addition, the "Jump to frame" dialogue box now gives two different ways to
specify which frame to jump to.

The progress bar for indexing jobs run via Slurm has also been made more
realistic.


Julia bindings
--------------

An initial version of Julia bindings for "libcrystfel" has been introduced.
This will allow the various building blocks from CrystFEL (peak search,
indexing, spot prediction, reflection data management and so on) to be used
from a high-level programming language.   To get started, see
<https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/julia.rst>


API changes
-----------

Removed routines:
* `crystal_copy_deep`
* `crystal_get_cell_const`
* `crystal_get_reflections`
* `crystal_get_image`
* `crystal_get_image_const`
* `crystal_set_reflections`
* `crystal_set_image`

Added routines:
* `get_lattice_symmetry`
* `crystal_reliquish_cell`
* `image_add_crystal_refls`
* `free_reflistiterator`
* `name_equiv`
* `solve_inv`
* `set_mm_funcs` and `cfmalloc`/`cffree` etc.

Added type definitions:
* `struct crystal_refls`

Changed structure definitions:
* `image` - change of `crystals` to `struct crystal_refls`
* `image` - addition of `owns_peaklist`

Changed routine prototypes:
* `update_predictions` - addition of `image`
* `predict_to_res` - addition of `image`
* `calculate_partialities` - addition of `image`
* `peakfinder8` - return type is now `ImageFeatureList`
* `search_peaks` - return type is now `ImageFeatureList`
* `search_peaks_peakfinder8` - return type is now `ImageFeatureList`
* `search_peaks_peakfinder9` - return type is now `ImageFeatureList`
* `indexing_peak_check` - addition of `peaks`
* `validate_peaks` - addition of `peaks`
