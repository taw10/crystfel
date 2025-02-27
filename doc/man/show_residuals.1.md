% show_residuals(1)

NAME
====

show_residuals - read and display spot position residuals from Millepede data


SYNOPSIS
========

show_residuals -g _input.geom_ _millepede-files_


DESCRIPTION
===========

**show_residuals** reads the calibration data written by **indexamajig --mille**,
and displays the average Bragg peak position offsets (predicted minus observed
spot coordinates).

To interpret the Millepede data, the geometry file is required.  The geometry
file should match the one used for the **indexamajig** run that generated the
data.

The program will show the spot position offsets and the average (absolute)
excitation errors of the reflections, averaged across the detector's hierarchy
groups.  Averages will be calculated at every level of hierarchy found in the
Millepede files, limited by **indexamajig --max-mille-level**.


OPTIONS
=======

**-g** _input.geom_, **--geometry**=_input.geom_
: Specify the input geometry filename.


AUTHOR
======

This page was written by Thomas White.


REPORTING BUGS
==============

Report bugs to <taw@physics.org>, or visit <http://www.desy.de/~twhite/crystfel>.


COPYRIGHT AND DISCLAIMER
========================

Copyright © 2025 Deutsches Elektronen-Synchrotron DESY, a research centre of
the Helmholtz Association.

show_residuals, and this manual, are part of CrystFEL.

CrystFEL is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CrystFEL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CrystFEL.  If not, see <http://www.gnu.org/licenses/>.


SEE ALSO
========

**crystfel**(7), **indexamajig**(1), **align_detector**(1),
**crystfel_geometry**(5)
