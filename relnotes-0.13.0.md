Release notes for CrystFEL version 0.13.0
=========================================

Copyright © 2012-2026 Deutsches Elektronen-Synchrotron DESY,
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

This version of CrystFEL improves support for serial electron crystallography,
as well as simplifying indexing options, figure of merit calculation and data
export.  There were also many smaller fixes and enhancements in other areas.
See the ChangeLog or the Git history for a comprehensive list of all changes.


TIFF support and improvements for electron diffraction
------------------------------------------------------

CrystFEL 0.13.0 adds support for reading TIFF images, as commonly output by
electron microscope detectors.  In addition, several parts of the indexing
system have been tuned for better success rates when working on electron
diffraction patterns.  Further improvements are planned for upcoming versions.


Removal of latt/cell indexer options
------------------------------------

The options for indexing patterns have been greatly simplified by removing the
`latt` and `cell` switches.  These were previously used to control what prior
information was given to the indexing algorithm.  Now, the indexing algorithms
will always make the maximum possible use of whatever information they are
given via the target unit cell.


Automatic unit cell for figures of merit and export
---------------------------------------------------

Partialator can now calculate the average cell parameters from all crystals
contributing to the merged dataset.  On the command line, you can get the
average unit cell by using the new option `partialator --output-cell`.

These parameters will automatically be used by the GUI for the figure of merit
calculation and data export, saving you a step.  You can override this if you
prefer by selecting your own unit cell file in the usual way.

In addition, the figure of merit window in the GUI now features automatic
generation of "Table 1".


GUI project now automatically saves
-----------------------------------

The GUI will now automatically save the session at certain suitable moments,
such as when starting a new indexing or merging job.


New Homebrew tap
----------------

We now have a private [Homebrew](https://brew.sh) repository ("tap"), greatly
simplifying installation primarily on Mac OS, but also on other platforms
supported by Homebrew.


API changes
-----------

Removed routines:
* `data_template_get_2d_detgeom_if_possible`
* `stream_get_fd`
* `base_indexer_str`

Added routines:
* `write_cell_to_file`
* `fom_free`

Changes function prototypes:
* `indexer_str` now returns `const`

Changed type definitions:
* `DataSourceType` - addition of `DATA_SOURCE_TYPE_TIFF` and `DATA_SOURCE_TYPE_RAW`
* Removal of `INDEXING_DEFAULTS_*`
* `IndexingMethod` - removal of `INDEXING_USE_LATTICE_TYPE` and `INDEXING_USE_CELL_PARAMETERS`
* Removal of `INDEXING_METHOD_MASK`
* `struct pinkindexer_options` - removal of (vestigial) `customBandwidth`
