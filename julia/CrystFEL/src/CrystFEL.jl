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

include("detgeom.jl")
using .DetGeoms
export Panel, DetGeom

include("symmetry.jl")
using .Symmetry
export SymOpList

include("datatemplates.jl")
using .DataTemplates
export DataTemplate, loaddatatemplate

include("image.jl")
using .Images
export Image

include("reflists.jl")
using .RefLists
export RefList, loadreflist
export Reflection, MergedReflection, UnmergedReflection

end  # of module