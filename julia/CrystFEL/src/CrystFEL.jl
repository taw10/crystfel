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

# Configure libcrystfel to use Julia's memory management.  This is needed so
# that the Julia GC knows about the memory we allocate via libcrystfel
# routines.  Otherwise, potentially very large objects will be kept hanging
# around in memory because Julia thinks it's using a very small amount of
# memory, and rarely runs the GC.  In the case of image structures, the
# difference between apparent and true memory use can be a factor of a million!
function __init__()
    @ccall libcrystfel.set_mm_funcs(cglobal(:jl_malloc)::Ptr{Cvoid},
                                    cglobal(:jl_free)::Ptr{Cvoid},
                                    cglobal(:jl_calloc)::Ptr{Cvoid},
                                    cglobal(:jl_realloc)::Ptr{Cvoid})::Cint
end

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
export SymOpList, asymmetricindices

include("datatemplates.jl")
using .DataTemplates
export DataTemplate, loaddatatemplate, wavelength, cameralength, translategroup!

include("peaklist.jl")
using .PeakLists
export PeakList

include("reflists.jl")
using .RefLists
export RefList, loadreflist
export Reflection, MergedReflection, UnmergedReflection

include("crystal.jl")
using .Crystals
export Crystal, InternalCrystal

include("image.jl")
using .Images
export Image

include("diffcalc.jl")
using .DiffractionCalculations
export predictreflections, calculatepartialities!
export PartialityModel, UnityModel, XSphereModel, OffsetModel, RandomModel, GeneralGaussianModel

include("indexing.jl")
using .Indexing
export Indexer, index

include("streams.jl")
using .Streams
export Stream, chunkwrite, chunkread

include("millepede.jl")
using .Millepede
export Mille

include("peaksearch.jl")
using .PeakSearch
export zaefpeaks, peakfinder8, peakfinder9

end  # of module
