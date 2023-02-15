Release notes for CrystFEL version 0.10.2
=========================================

Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
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

This release is primarily a bug-fixing update to CrystFEL, and there have been
improvements around the entire suite.  The most important changes are detailed
below.  See the ChangeLog or the Git history for a comprehensive list of all
changes.


Graphical user interface (GUI)
------------------------------

This version adds a colour scale widget to the GUI, which can be used to
manually adjust the colours used to display the image.  The colour scale is
found directly to the right of the image.  Click and drag to move
the histogram of image values relative to the colour scale, and use the scroll
wheel to compress and expand the range.

Resolution rings can now be added to the image.  Enable them in the "View" menu.
There are also menu items for running `detector-shift` and `peakogram`, found
under the "Tools" menu.

The figure of merit calculation tool has been improved, and now displays the
results as a graph.  The graph values can be exported as a CSV file for further
processing.

This release also adds specially-designed icons to the GUI, which fixes
problems with "broken" icons on some systems.


Performance
-----------

Speed has been improved in many areas of CrystFEL.  The profiling option
(`--profile`) of `indexamajig` has been re-implemented and can help you to find
bottlenecks in processing on your data.  In addition, some new options were
added to `indexamajig`:

New option `--peakfinder8-fast` tells peakfinder8 to calculate some values in
advance, approximately halving the time taken by the peak search.

Another new option `--asdf-fast` slightly alters some thresholds within asdf,
resulting in a large speedup.  Since this changes the behaviour from before, we
decided to make it optional.

Finally, `--cell-parameters-only` completely circumvents the spot prediction and
integration calculations, which can make a very large speedup as long as you're
not interested in spot positions.  Note that this isn't the same as the
existing option ``--integration=none``, where the spot positions are predicted
but the intensity integration is skipped.  For very large unit cells, which are
often produced by unconstrained indexing (e.g. for online monitoring), the spot
prediction can take signficant time.

For advice about making CrystFEL faster, see this page:  
<https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/speed.rst>


Documentation improvements
--------------------------

Two long-form documentation articles have been added, on top of the other two
mentioned in other sections of this page.

A completely new tutorial has been added, describing how to process data using
the graphical user interface:  
<https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/tutorial.rst>

The following new article has been added, to help with the common question of
choosing the correct point group for merging data:  
<https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/pointgroup.rst>


Streaming interfaces
--------------------

ZeroMQ and ASAP::O interfaces have been added, to allow data to be streamed
directly to CrystFEL without any intermediate file storage.  For some
information on how to get started, read the following document:  
<https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/online.rst>


Refuses to overwrite stream
---------------------------

Starting from this version, `indexamajig` and other CrystFEL tools will refuse to
create a new stream if a file already exists with the same name.  This is meant
to avoid the possibility of losing valuable results by (e.g.) re-running a
processing script.  This becomes particularly important in a streaming data
situation, where the original data might no longer be available.  If you really
want to overwrite the results, you'll need to move or delete the old file
yourself.


Installation
------------

CrystFEL has become easier to install because of several developments.  First,
the Slurm API headers are no longer needed for the GUI to submit Slurm jobs.
In this version, it uses the Slurm commands directly
(`sbatch/scontrol/scancel`).

For environments based totally on CBFs or other types of files, it's now
possible to compile CrystFEL without any reference to the HDF5 library,
avoiding a lot of possible complications.

There is a new script (`scripts/install-indexers`) which helps with installing
Mosflm, DirAx and XDS.  After downloading CrystFEL, run the script with
`--help` to get started.

Finally, Docker images are now available.  These can be used with various
container tools (not just Docker).  For example, using Apptainer (the new name
for Singularity):
```
$ apptainer pull docker://gitlab.desy.de:5555/thomas.white/crystfel/crystfel:latest
$ apptainer run -B /path/to/data crystfel_latest.sif
```

More details can be found on the installation page:
<https://www.desy.de/~twhite/crystfel/install.html>


API changes
-----------

Added routines:
* `data_template_slabby_file_to_panel_coords`
* `data_template_get_2d_detgeom_if_possible`
* `data_template_get_clen_if_possible`
* `image_create_dp_bad_sat`
* `profile_{init,start,end,print_and_reset}`

Removed routines:
* `image_set_zero_mask`

Changed routine prototypes:
* `default_method_options` - addition of `asdf_opts`
* `peakfinder8` - addition of `fast_mode` and `private_data`

Changed enumerations:
* `DataSourceType` - addition of `DATA_SOURCE_TYPE_SEEDEE`
