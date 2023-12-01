==========================
Using libCrystFEL in Julia
==========================

CrystFEL's shared library component, *libcrystfel*, can be accessed from
`Julia <https://www.julialang.org/>`_.  This way, you can carry out almost any
process available in the CrystFEL tools, and much more besides.

The Julia package for CrystFEL is found in the ``julia`` directory of the
CrystFEL source code package.  The easiest way to get started is to use
``Pkg.develop``.  Start Julia, press ``]``, then run the following command,
substituting the correct location of the CrystFEL package on your system::

    (@v1.9) pkg> dev /home/twhite/crystfel/julia/CrystFEL

Afterwards, import the CrystFEL package as follows::

    julia> using CrystFEL

It's imperative that the same version of the Julia package is used as CrystFEL
itself, otherwise you are likely to experience spontaneous crashes of the
whole Julia session.

You only need to perform the ``pkg dev`` process once, not every session.  Should
you ever need to, you can remove the CrystFEL package in the usual way::

    (@v1.9) pkg> rm CrystFEL

The Julia package for CrystFEL is a fairly thin wrapper around the
`libcrystfel C API <https://www.desy.de/~twhite/crystfel/reference/index.html>`_,
so you can read the API documentation for some idea of the possibilities.
Below is a brief overview.


Unit cells
==========

Create a ``UnitCell`` object as follows::

    cell = UnitCell(MonoclinicLattice, PrimitiveCell, UniqueAxisB, 123, 45, 80, 90, 97, 90)

The arguments are, in order:

* Lattice type, one of ``TriclinicLattice``, ``MonoclinicLattice``,
  ``OrthorhombicLattice``, ``TetragonalLattice``, ``HexagonalLattice``,
  ``RhombohedralLattice`` or ``CubicLattice``.
* Centering, one of ``PrimitiveCell`` (P), ``ACenteredCell``, ``BCenteredCell``,
  ``CCenteredCell``, ``BodyCenteredCell`` (I), ``FaceCenteredCell`` (F),
  ``RhombohedralCell`` (R) or ``RhombohedralCellOnHexagonalAxes`` ("H").
* Unique axis, one of ``UniqueAxisA``, ``UniqueAxisB``, ``UniqueAxisC``,
  ``NoUniqueAxis`` (for lattice types that do not have a unique axis, e.g. cubic),
  or ``UnknownUniqueAxis`` (which should be avoided, but is sometimes necessary).
* Cell axis lengths a,b,c, in Angstroms.
* Cell angles α,β,γ, in degrees.

For many cases, you don't need to specify every parameter.  Where possible, the
unique axis will be determined from the cell parameters, and can be omitted::

    cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)

Cell axis lengths, angles and centering types can be omitted if they are fixed
by the lattice type.  For example::

    julia> UnitCell(CubicLattice, FaceCenteredCell, 40)
    UnitCell(CubicLattice, FaceCenteredCell, NoUniqueAxis,
         40.0 Å, 40.0 Å, 40.0 Å, 90.0°, 90.0°, 90.0°)

or::

    julia> UnitCell(RhombohedralLattice, 23, 75)
    UnitCell(RhombohedralLattice, RhombohedralCell, NoUniqueAxis,
         23.0 Å, 23.0 Å, 23.0 Å, 75.0°, 75.0°, 75.0°)


Reflection lists
================

In Julia, a distinction is made between merged and unmerged reflections.  No
such distinction exists in CrystFEL's C API.  Merged reflections have
multiplicities (number of contributing measurements), whereas unmerged
reflections have detector locations, background levels and parameters related
to diffraction geometry such as excitation errors and Lorentz factors.  In
reality, both types of reflection have all fields, just as for the C API.  The
distinction controls how the objects are printed, and may help in writing
clearer programs.

Load a reflection list from a data file (".hkl file") using
``loadreflist``::

    julia> q = loadreflist("example.hkl")
    Merged reflection list in point group mmm
       h    k    l  intensity  σ(intens) nmeas
       0    0    5      23.45     124.51    11
       0    0    6   32302.69   10091.34     8
       0    0    7     -87.37     167.23     5
       0    0    8    8051.75    3828.24     6
       0    0    9      94.13     128.59     3
       0    0   10    3703.07    1118.85     5
       0    0   11       4.81      31.46     5
       0    0   12   27287.94   14143.77     6
       0    0   13     148.50      32.46     2
       0    0   14     818.33     447.23     3
       ⋮    ⋮    ⋮          ⋮          ⋮     ⋮

This produces a merged reflection list::

    julia> typeof(q)
    RefList{MergedReflection}

Unmerged reflection lists (``RefList{UnmergedReflection}``) come from spot
prediction - see below.

Reflection lists can be iterated over.  For example::

    julia> for refl in q
               println(refl)
           end
    MergedReflection((0, 0, 5), intensity=23.450000762939453, σ(intensity)=124.51000213623047, nmeasurements=11)
    MergedReflection((0, 0, 6), intensity=32302.689453125, σ(intensity)=10091.33984375, nmeasurements=8)
    MergedReflection((0, 0, 7), intensity=-87.37000274658203, σ(intensity)=167.22999572753906, nmeasurements=5)
    MergedReflection((0, 0, 8), intensity=8051.75, σ(intensity)=3828.239990234375, nmeasurements=6)
    MergedReflection((0, 0, 9), intensity=94.12999725341797, σ(intensity)=128.58999633789062, nmeasurements=3)
    MergedReflection((0, 0, 10), intensity=3703.070068359375, σ(intensity)=1118.8499755859375, nmeasurements=5)
    MergedReflection((0, 0, 11), intensity=4.809999942779541, σ(intensity)=31.459999084472656, nmeasurements=5)
    MergedReflection((0, 0, 12), intensity=27287.939453125, σ(intensity)=14143.76953125, nmeasurements=6)
    MergedReflection((0, 0, 13), intensity=148.5, σ(intensity)=32.459999084472656, nmeasurements=2)
    MergedReflection((0, 0, 14), intensity=818.3300170898438, σ(intensity)=447.2300109863281, nmeasurements=3)
    ...

You can subscript a RefList using Miller indices::

    julia> q[1, 13, 43]
    MergedReflection((1, 13, 43), intensity=-28.18000030517578, σ(intensity)=7.230000019073486, nmeasurements=6)

Linear indexing is **not** supported, so you **can't** do things like
``q[10:end]``.


Symmetry
========

Symmetry operations are represented by ``SymOp`` objects, which are contained
within ``SymOpList`` objects.  A point group is therefore represented by a
``SymOpList``.

Create a point group from the Herman-Mauguin symbol as follows::

    julia> s = SymOpList("2/m")
    4-element SymOpList ("2/m")
    -h,-k,l
    h,k,-l
    hkl
    -h,-k,-l

The list can be subscripted linearly::

    julia> s[1]
    SymOp("-h,-k,l")

Note that not all ``SymOpList`` objects represent a point group.  One
counterexample is lists of indexing ambiguity operations.


Peak lists
==========


Crystals
========


Indexing
========


Prediction
==========

