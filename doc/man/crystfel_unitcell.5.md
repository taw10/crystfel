% crystfel_unitcell(5)

INTRODUCTION
============

A CrystFEL "unit cell file", usually named _something.cell_, describes the
unit cell of your crystals: the Bravais lattice type, the centering, and
optionally the cell axis lengths and angles.

Unit cell files are accepted anywhere CrystFEL asks for a unit cell, for
example **indexamajig -p**, **compare_hkl -p** and **cell_tool**.  A unit cell
file can be created by hand in a text editor, or written by **cell_tool** with
its **-o** option.

A unit cell file does not have to contain the cell parameters.  If you give
only the lattice type and centering, the file describes a _lattice_ rather
than a specific cell, which allows you to say "index these patterns using
tetragonal primitive lattices, but with any parameters".

Wherever a unit cell file can be used, a PDB file can be used instead.  See
section **PDB FILES** below.


SYNTAX
======

The first line of a unit cell file must be exactly this:

    CrystFEL unit cell file version 1.0

If the first line is anything else, the file will be treated as a PDB file.

The remaining lines have the following general form:

    parameter = value

Values which represent a physical quantity must be followed by their units,
separated by a space:

    parameter = value units

Fields can be given in any order.  Blank lines are ignored, and so are lines
beginning with a semicolon (**;**), which you can use for comments.
Unrecognised field names will be reported as errors, but will not stop the
rest of the file from being read.


LATTICE TYPE AND CENTERING
==========================

**lattice_type** = _type_
: The Bravais lattice type.  One of **triclinic**, **monoclinic**,
  **orthorhombic**, **tetragonal**, **rhombohedral**, **hexagonal** or
  **cubic**.

**centering** = _symbol_
: The centering of the lattice, as a single letter.  One of **P** (primitive),
  **A**, **B** or **C** (base-centered on the respective face), **I**
  (body-centered), **F** (face-centered), **R** (rhombohedral, on rhombohedral
  axes) or **H** (rhombohedral, on hexagonal axes).

**unique_axis** = _axis_
: The unique axis of the cell: **a**, **b** or **c**.

The unique axis **must** be given for **monoclinic**, **tetragonal** and
**hexagonal** lattice types, and the file will be rejected if you leave it
out.  It is not applicable to the other lattice types, and giving it for
**triclinic**, **orthorhombic**, **rhombohedral** or **cubic** will produce a
warning.

For a monoclinic lattice with **A**, **B** or **C** centering, the unique axis
must be different from the centering.  For example, a monoclinic C-centered
cell must have unique axis **a** or **b**, not **c**.

Note that a rhombohedral lattice can be described in two ways.  On
rhombohedral axes, use **lattice_type = rhombohedral** with **centering = R**.
On hexagonal axes, use **lattice_type = hexagonal** with **centering = H** and
**unique_axis = c**.  Structures deposited in the PDB usually use hexagonal
axes.  Running **cell_tool --uncenter** on a hexagonal **H** cell will give
you the equivalent rhombohedral **R** cell.


CELL PARAMETERS
===============

**a** = _nnn_ [**A**|**nm**], **b** = _nnn_ [**A**|**nm**], **c** = _nnn_ [**A**|**nm**]
: The lengths of the unit cell axes.  The units must be given, and must be
  **A** (Angstroms) or **nm** (nanometres).

**al** = _nnn_ [**deg**|**rad**], **be** = _nnn_ [**deg**|**rad**], **ga** = _nnn_ [**deg**|**rad**]
: The unit cell angles alpha, beta and gamma.  The units must be given, and
  must be **deg** (degrees) or **rad** (radians).

All six parameters must be given. A partially-specified cell is
treated as if none of the parameters had been given.

EXAMPLES
========

A hexagonal cell on hexagonal axes, with all parameters given:

    CrystFEL unit cell file version 1.0

    lattice_type = hexagonal
    centering = H
    unique_axis = c

    a = 66.2 A
    b = 66.2 A
    c = 150.2 A

    al = 90.0 deg
    be = 90.0 deg
    ga = 120.0 deg

A cubic body-centered cell, with the parameters given in nanometres:

    CrystFEL unit cell file version 1.0

    lattice_type = cubic
    centering = I

    a = 6.62 nm
    b = 6.62 nm
    c = 6.62 nm

    al = 90.0 deg
    be = 90.0 deg
    ga = 90.0 deg

A tetragonal primitive lattice with no parameters, for indexing patterns when
you know the lattice type but not the cell dimensions:

    CrystFEL unit cell file version 1.0

    ; Lattice type only - any cell dimensions will be accepted
    lattice_type = tetragonal
    centering = P
    unique_axis = c

For more examples, look in the **examples** folder, which can be found online
at https://gitlab.desy.de/thomas.white/crystfel/-/tree/master/doc/examples


PDB FILES
=========

Anywhere a unit cell file is accepted, you can give a PDB file instead.  The
cell parameters will be taken from the **CRYST1** record, including the space
group, from which the lattice type and centering will be deduced.

Reading the cell from a PDB file is convenient, but the unit cell file format
is preferred because it is explicit: it does not require the space group to be
interpreted, and it allows a lattice to be specified without any parameters.


REPORTING BUGS
==============

Report bugs to <taw@physics.org>, or visit <http://www.desy.de/~twhite/crystfel>.


COPYRIGHT AND DISCLAIMER
========================

Copyright © 2026 Deutsches Elektronen-Synchrotron DESY, a research centre of
the Helmholtz Association.

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

**crystfel**(7), **crystfel_geometry**(5), **indexamajig**(1),
**cell_tool**(1), **cell_explorer**(1)
