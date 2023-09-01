% adjust_detector(1)

NAME
====

adjust_detector - move detector panels


SYNOPSIS
========

adjust_detector -g _input.geom_ -o _output.geom_ -p _group_ [_movement_]


DESCRIPTION
===========

**adjust_detector** moves a panel (or group of panels) in a CrystFEL detector
geometry file, and writes an updated geometry file.

To rotate a panel, use one of **--rotx**, **--roty** or **--rotz**.  The
rotation will be in degrees, clockwise when looking along the specified axis.
The center of rotation will be the centroid of the corners of the panel, or the
centroid of all the panel centers in the specified group.

To translate a panel or group, use one of **--shiftx**, **--shifty** or
**--shiftz**.  The translation will be along the positive direction of the
specified axis.  The units are pixels, unless you additionall specify
**--mm**.

You can use multiple movements together, but the results are unspecified if you
combine a rotation with any other movement.


OPTIONS
=======

**-g** _input.geom_
: Specify the input geometry filename.

**-o** _output.geom_
: Specify the output geometry filename.
: Note that the geometry file will be re-written, meaning that any formatting
: and comments will be lost.

**-p** _group_
: Specify the panel, or panel group, to move.

**--mm**
: Interpret panel shifts as millimetres, not pixels.

**--shiftx=***n*, **--shifty=***n*, **--shiftz=***n*
: Shift the chosen panel (or group) by _n_ along the given direction, in units
: of pixels (unless **--mm** is used).  When moving a group of panels, the
: pixel size of the *first* panel will be used.

**--rotx=***n*, **--roty=***n*, **--rotz=***n*
: Shift the chosen panel (or group) by _n_ degrees clockwise about the given
: axis.  The center of rotation will be the centroid of the corners of the
: panel, or the centroid of all the panel centers in the specified group.


AUTHOR
======

This page was written by Thomas White.


REPORTING BUGS
==============

Report bugs to <taw@physics.org>, or visit <http://www.desy.de/~twhite/crystfel>.


COPYRIGHT AND DISCLAIMER
========================

Copyright © 2023 Deutsches Elektronen-Synchrotron DESY, a research centre of
the Helmholtz Association.

adjust_detector, and this manual, are part of CrystFEL.

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
