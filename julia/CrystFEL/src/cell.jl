module UnitCells

import ..CrystFEL: libcrystfel
export UnitCell, LatticeType, CenteringType, UniqueAxis
export TriclinicLattice, MonoclinicLattice, OrthorhombicLattice
export TetragonalLattice, HexagonalLattice, RhombohedralLattice, CubicLattice
export PrimitiveCell, ACenteredCell, BCenteredCell, CCenteredCell
export BodyCenteredCell, FaceCenteredCell, RhombohedralCell, RhombohedralCellOnHexagonalAxes
export NoUniqueAxis, UnknownUniqueAxis, UniqueAxisA, UniqueAxisB, UniqueAxisC


# Represents the real C-side (opaque) structure.
mutable struct InternalUnitCell end

# The Julia-side structure, needed to house the pointer to the C structure
# Without this, we would only ever have a Ptr{DataTemplate}, not a DataTemplate.
mutable struct UnitCell
    internalptr::Ptr{InternalUnitCell}
end


"""
Enumeration of the seven Bravais lattice types: `TriclinicLattice`,
`MonoclinicLattice`, `OrthorhombicLattice`, `TetragonalLattice`,
`RhombohedralLattice`, `HexagonalLattice`, `CubicLattice`.
"""
@enum LatticeType begin
    TriclinicLattice
    MonoclinicLattice
    OrthorhombicLattice
    TetragonalLattice
    RhombohedralLattice
    HexagonalLattice
    CubicLattice
end


"""
Enumeration of unit cell centering possibilities: `PrimitiveCell`,
`ACenteredCell`, `BCenteredCell`, `CCenteredCell`, `BodyCenteredCell`
(I-centering), `FaceCenteredCell` (F-centering), and `RhombohedralCell`
(R-centering, primitive rhombohedral cell).  `RhombohedralCellOnHexagonalAxes`
indicates "H-centering" as used by the protein data bank, which is different
from the "triple hexagonal cell" described in the International Tables.
"""
@enum CenteringType begin
    PrimitiveCell = Int('P')
    ACenteredCell = Int('A')
    BCenteredCell = Int('B')
    CCenteredCell = Int('C')
    BodyCenteredCell = Int('I')
    FaceCenteredCell = Int('F')
    RhombohedralCell = Int('R')
    RhombohedralCellOnHexagonalAxes = Int('H')
end


"""
Enumeration of unique axis possibilities.  The possibilities are `UniqueAxisA`,
`UniqueAxisB` and `UniqueAxisC`.  Alternatively, use `NoUniqueAxis` if the type
of unit cell does not have a unique axis (triclinic, orthorhombic, cubic or
rhombohedral).  `UnknownUniqueAxis` means that the unique axis is not known.
"""
@enum UniqueAxis begin
    NoUniqueAxis = Int('*')
    UnknownUniqueAxis = Int('?')
    UniqueAxisA = Int('a')
    UniqueAxisB = Int('b')
    UniqueAxisC = Int('c')
end


"""
    UnitCell(latticetype, centering, uniqueaxis, a, b, c, α, β, γ)

Creates a CrystFEL UnitCell, in an undefined orientation, from the given
parameters.  The axis lengths (a, b, c) should be in *Ångstroms*, and the
angles (α, β, γ) should be in *degrees* - note that this is different to the
equivalent function in CrystFEL's C API.

See the documentation for `LatticeType`, `CenteringType` and `UniqueAxis` for
possible values.  You can also use the characters `'a'`, `'b'` and `'c'` for
`uniqueaxis`, or `'P'`, `'A'`, `'B'`, `'C'`, `'I'`, `'F'`, `'R'` and `'H`' for
`centering.

Corresponds to CrystFEL C API function `cell_new_from_parameters` with follow-up
calls to `cell_set_centering`, `cell_set_lattice_type` and `cell_set_unique_axis`.
"""
function UnitCell(latticetype, centering, uniqueaxis, a, b, c, α, β, γ)

    out = ccall((:cell_new_from_parameters, libcrystfel),
                Ptr{InternalUnitCell},
                (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),
                a/1e10, b/1e10, c/1e10, deg2rad(α), deg2rad(β), deg2rad(γ))
    if out == C_NULL
        throw(ArgumentError("Failed to create unit cell"))
    end

    ccall((:cell_set_centering, libcrystfel),
          Cvoid, (Ptr{InternalUnitCell},Cchar),
          out, centering)

    ccall((:cell_set_unique_axis, libcrystfel),
          Cvoid, (Ptr{InternalUnitCell},Cchar),
          out, uniqueaxis)

    ccall((:cell_set_lattice_type, libcrystfel),
          Cvoid, (Ptr{InternalUnitCell},Cint),
          out, latticetype)

    cell = UnitCell(out)

    finalizer(cell) do x
        ccall((:cell_free, libcrystfel),
              Cvoid, (Ptr{InternalUnitCell},), x.internalptr)
    end

    return cell
end

UnitCell(a, b, c, α, β, γ) = UnitCell(TriclinicLattice, PrimitiveCell, UnknownUniqueAxis,
                                      a, b, c, α, β, γ)


function getlatticetype(cell)
    lt = ccall((:cell_get_lattice_type, libcrystfel),
               Cint, (Ptr{InternalUnitCell},), cell.internalptr)
    cen = ccall((:cell_get_centering, libcrystfel),
                Cchar, (Ptr{InternalUnitCell},), cell.internalptr)
    ua = ccall((:cell_get_unique_axis, libcrystfel),
               Cchar, (Ptr{InternalUnitCell},), cell.internalptr)
    return LatticeType(lt),CenteringType(cen),UniqueAxis(ua)
end


function getcellparams(cell)
    let a=Ref{Cdouble}(0),
        b=Ref{Cdouble}(0),
        c=Ref{Cdouble}(0),
        α=Ref{Cdouble}(0),
        β=Ref{Cdouble}(0),
        γ=Ref{Cdouble}(0)
        ccall((:cell_get_parameters, libcrystfel),
              Cvoid, (Ptr{InternalUnitCell},
                      Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},
                      Ref{Cdouble},Ref{Cdouble},Ref{Cdouble}),
              cell.internalptr, a, b, c, α, β, γ)
        return (a=a[], b=b[], c=c[], α=α[], β=β[], γ=γ[])
    end
end


function Base.show(io::IO, uc::UnitCell)
    write(io, "UnitCell(")
    lt,cen,ua = getlatticetype(uc)
    show(io, lt); write(io, ", ")
    show(io, cen); write(io, ", ")
    show(io, ua); write(io, ",\n         ")
    let cp = getcellparams(uc)
        show(io, cp.a*1e10); write(io, " Å, ")
        show(io, cp.b*1e10); write(io, " Å, ")
        show(io, cp.c*1e10); write(io, " Å, ")
        show(io, rad2deg(cp.α)); write(io, "°, ")
        show(io, rad2deg(cp.β)); write(io, "°, ")
        show(io, rad2deg(cp.γ)); write(io, "°")
    end
    write(io, ")")
end


end   # of module
