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

A ``UnitCell`` also keeps a record of the orientation of the lattice.  You can
get the real and reciprocal space matrices of lattice vectors with
``cell.directcartesian`` and ``cell.reciprocalcartesian``, respectively::

    julia> uc.directcartesian
    3×3 Matrix{Float64}:
     1.37e-8  8.38883e-25  8.38883e-25
     0.0      1.37e-8      8.38883e-25
     0.0      0.0          1.37e-8


Reflection lists
================

A ``RefList`` is a container for reflection data. In Julia, a distinction is
made between merged and unmerged reflections. Merged reflections have
multiplicities (number of contributing measurements), whereas unmerged
reflections have detector locations, background levels and parameters related
to diffraction geometry such as excitation errors and
Lorentz factors.

No such distinction exists in CrystFEL's C API, and in "reality", both types of
reflection have all fields.  The distinction controls how the objects are
printed, and may help in writing clearer programs.

Load a reflection list from a data file (".hkl file") using ``loadreflist``::

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
``SymOpList``, but note that not all ``SymOpList`` objects represent a symmetry
group (in the sense of group theory).  One counterexample is lists of indexing
ambiguity operations.

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


Images and DataTemplates
========================

A ``DataTemplate`` represents the contents of a CrystFEL geometry file, which
describes the layout of information in the data, the physical positions of
parts of the detector, and the values of various items of metadata (or
information about where to get those values).  Create a ``DataTemplate`` by
loading a geometry file::

    dtempl = loaddatatemplate("/path/to/my.geom")

An ``Image`` is an overall container structure representing one frame of a
serial crystallography dataset.  Create one by loading an image from file::

    image = Image(dtempl, "/path/to/mydata.cxi", "//32")

You can use any kind of file supported by CrystFEL here.  In the example,
``//32`` is the frame ID - leave it out if there is only one frame per file.

If you're simulating data, you can create an empty image like this::

    image = Image(dtempl)

However, several caveats apply to doing this.  The ``DataTemplate`` must not
say that any metadata values (e.g. the wavelength) should be taken from file
headers, because there is no file in this case.  An error will be thrown if
there is any problem.


Peak lists
==========

A ``PeakList`` represents a list of positions on the detector surface.  Create
it and add peaks like this::

    peaklist = PeakList()
    push!(peaklist, 10.0, 20.0, 1, 2000.0)

The arguments to ``push!(::PeakList, ...)`` are, in order, the fast scan
coordinate, slow scan coordinate (both relative to the panel corner), panel
number (indexed from zero) and the spot intensity in detector units.

You can assign your peaklist to an ``Image`` by setting ``image.peaklist``.
Note that any ``PeakList`` can only be assigned to a single ``Image``.  An
error will be thrown if you try to add the same ``PeakList`` again (even to the
same ``Image``).  If necessary, you can make a copy using ``deepcopy``.


Crystals
========

A ``Crystal`` is made up of a ``UnitCell`` (which includes the orientation of
the crystal, as mentioned above) and a few other parameters::

    julia> cr = Crystal(uc)
    CrystFEL.Crystal(0x000000001ad0ada0):

    UnitCell(CubicLattice, FaceCenteredCell, NoUniqueAxis,
             137.000 Å, 137.000 Å, 137.000 Å, 90.000°, 90.000°, 90.000°)

     Linear scale factor: 1.0
      Debye-Walle factor: 0.0
               Mosaicity: 0.0
          Profile radius: 0.002 nm⁻¹
        Resolution limit: Inf
                    Flag: 0


Indexing
========

Create an indexing engine like this::

    indexer = Indexer("asdf", dtempl, cell)

The first argument is handled exactly as by ``indexamajig --indexing=``.  The
second argument is the detector geometry (``DataTemplate``), and you also need
to provide the target unit cell.

Run the indexing engine on an image like this::

    index(image, indexer)


Prediction
==========

Given a ``Crystal`` and and ``Image``, the function
``predictreflections(cr, image)`` will return a ``RefList{UnmergedReflection}``
containing the predicted reflections.

You can subsequently calculate partialities with::

    calculatepartialities!(reflist, cr, image, model=XSphereModel)

