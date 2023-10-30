module UnitCells

import ..CrystFEL: libcrystfel
export UnitCell, LatticeType
export TriclinicLattice, MonoclinicLattice, OrthorhombicLattice
export TetragonalLattice, HexagonalLattice, RhombohedralLattice, CubicLattice


# Represents the real C-side (opaque) structure.
mutable struct InternalUnitCell end

# The Julia-side structure, needed to house the pointer to the C structure
# Without this, we would only ever have a Ptr{DataTemplate}, not a DataTemplate.
mutable struct UnitCell
    internalptr::Ptr{InternalUnitCell}
end

@enum LatticeType begin
    TriclinicLattice
    MonoclinicLattice
    OrthorhombicLattice
    TetragonalLattice
    HexagonalLattice
    RhombohedralLattice
    CubicLattice
end

"""
    UnitCell(a, b, c, α, β, γ, centering='P', latticetype=TriclinicLattice)

Creates a CrystFEL UnitCell, in an undefined orientation, from the given
parameters.  The angles (α, β, γ) should be in *degrees* - note that this is
different to the equivalent function in CrystFEL's C API.

Corresponds to CrystFEL C API function `cell_new_from_parameters()`.
"""
function UnitCell(a, b, c, α, β, γ,
                  centering='P',
                  latticetype=TriclinicLattice,
                  uniqueaxis='*')

    out = ccall((:cell_new_from_parameters, libcrystfel),
                Ptr{InternalUnitCell},
                (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),
                a, b, c, deg2rad(α), deg2rad(β), deg2rad(γ))
    if out == C_NULL
        throw(ArgumentError("Failed to create unit cell"))
    end

    ccall((:cell_set_centering, libcrystfel),
          Cvoid, (Ptr{InternalUnitCell},Cchar),
          out, centering)

    ccall((:cell_set_unique_axis, libcrystfel),
          Cvoid, (Ptr{InternalUnitCell},Cchar),
          out, uniqueaxis)

    cell = UnitCell(out)

    finalizer(cell) do x
        ccall((:cell_free, libcrystfel),
              Cvoid, (Ptr{InternalUnitCell},), x.internalptr)
    end

    return cell
end


function Base.show(io::IO, uc::UnitCell)
    write(io, "UnitCell()")
end


end   # of module
