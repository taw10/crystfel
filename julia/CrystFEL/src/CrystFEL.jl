"""
    CrystFEL

Julia bindings for CrystFEL data structures and routines

## Quick start
```julia
  using CrystFEL
  ...
```
"""
module CrystFEL

libcrystfel = "libcrystfel.so"

include("cell.jl")
using .UnitCells
export UnitCell, LatticeType, CenteringType, UniqueAxis
export TriclinicLattice, MonoclinicLattice, OrthorhombicLattice
export TetragonalLattice, HexagonalLattice, RhombohedralLattice, CubicLattice
export PrimitiveCell, ACenteredCell, BCenteredCell, CCenteredCell
export BodyCenteredCell, FaceCenteredCell, RhombohedralCell, RhombohedralCellOnHexagonalAxes
export NoUniqueAxis, UnknownUniqueAxis, UniqueAxisA, UniqueAxisB, UniqueAxisC
export rotatecell

include("detgeom.jl")
using .DetGeoms
export Panel, DetGeom

include("symmetry.jl")
using .Symmetry
export SymOpList

include("datatemplates.jl")
using .DataTemplates
export DataTemplate, loaddatatemplate, wavelength, cameralength

include("peaklist.jl")
using .PeakLists
export PeakList

include("image.jl")
using .Images
export Image

include("reflists.jl")
using .RefLists
export RefList, loadreflist
export Reflection, MergedReflection, UnmergedReflection

include("crystal.jl")
using .Crystals
export Crystal, InternalCrystal

include("diffcalc.jl")
using .DiffractionCalculations
export predictreflections

include("indexing.jl")
using .Indexing
export Indexer, index

include("streams.jl")
using .Streams
export Stream

include("millepede.jl")
using .Millepede
export Mille

end  # of module
