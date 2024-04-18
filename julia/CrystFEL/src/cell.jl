module UnitCells

using Random
using Printf

import ..CrystFEL: libcrystfel
export UnitCell, LatticeType, CenteringType, UniqueAxis
export TriclinicLattice, MonoclinicLattice, OrthorhombicLattice
export TetragonalLattice, HexagonalLattice, RhombohedralLattice, CubicLattice
export PrimitiveCell, ACenteredCell, BCenteredCell, CCenteredCell
export BodyCenteredCell, FaceCenteredCell, RhombohedralCell, RhombohedralCellOnHexagonalAxes
export NoUniqueAxis, UnknownUniqueAxis, UniqueAxisA, UniqueAxisB, UniqueAxisC
export rotatecell


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

    ninety(a) = isapprox(a, 90, atol=0.1)

    if latticetype == OrthorhombicLattice
        if !(ninety(α) && ninety(β) && ninety(γ))
            throw(ArgumentError("All angles must be 90°"))
        end

    elseif latticetype == TetragonalLattice
        if !(ninety(α) && ninety(β) && ninety(γ))
            throw(ArgumentError("All angles must be 90°"))
        end

    elseif latticetype == CubicLattice
        if !(ninety(α) && ninety(β) && ninety(γ))
            throw(ArgumentError("All angles must be 90°"))
        end

    elseif latticetype == HexagonalLattice
        if uniqueaxis == UniqueAxisA
            if !isapprox(b, c, rtol=0.01)
                throw(ArgumentError("b and c lengths should be equal"))
            end
        elseif uniqueaxis == UniqueAxisB
            if !isapprox(a, c, rtol=0.01)
                throw(ArgumentError("a and c lengths should be equal"))
            end
        elseif uniqueaxis == UniqueAxisC
            if !isapprox(a, b, rtol=0.01)
                throw(ArgumentError("a and b lengths should be equal"))
            end
        else
            throw(ArgumentError("Hexagonal cell requires a unique axis"))
        end
    end

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



"""
    UnitCell(latticetype, centering, a, b, c, α, β, γ)

A convenience constructor which attempts to determine the unique axis
automatically from the cell parameters.  If the unique axis is not obvious,
an `ArgumentError` will be thrown.
"""
function UnitCell(latticetype, centering, a, b, c, α, β, γ)

    notninety(a) = !isapprox(a, 90, atol=0.5)
    ninety(a) = isapprox(a, 90, atol=0.1)
    onetwenty(a) = isapprox(a, 120, atol=0.1)

    if latticetype == TriclinicLattice
        ua = NoUniqueAxis
    elseif latticetype == OrthorhombicLattice
        ua = NoUniqueAxis
    elseif latticetype == RhombohedralLattice
        ua = NoUniqueAxis
    elseif latticetype == CubicLattice
        ua = NoUniqueAxis

    elseif latticetype == MonoclinicLattice
        if notninety(α) && ninety(β) && ninety(γ)
            ua = UniqueAxisA
        elseif ninety(α) && notninety(β) && ninety(γ)
            ua = UniqueAxisB
        elseif ninety(α) && ninety(β) && notninety(γ)
            ua = UniqueAxisC
        else
            throw(ArgumentError("Can't determine unique axis"))
        end

    elseif latticetype == TetragonalLattice
        if isapprox(b, c, rtol=0.01) && !isapprox(a, b, rtol=0.05)
            ua = UniqueAxisA
        elseif isapprox(a, c, rtol=0.01) && !isapprox(a, b, rtol=0.05)
            ua = UniqueAxisB
        elseif isapprox(a, b, rtol=0.01) && !isapprox(c, b, rtol=0.05)
            ua = UniqueAxisC
        else
            throw(ArgumentError("Can't determine unique axis"))
        end

    elseif latticetype == HexagonalLattice
        if onetwenty(α) && ninety(β) && ninety(γ)
            ua = UniqueAxisA
        elseif ninety(α) && onetwenty(β) && ninety(γ)
            ua = UniqueAxisB
        elseif ninety(α) && ninety(β) && onetwenty(γ)
            ua = UniqueAxisC
        else
            throw(ArgumentError("Can't determine unique axis"))
        end
    end

    UnitCell(latticetype, centering, ua, a, b, c, α, β, γ)

end


"""
    UnitCell(latticetype, centering, a, b, c)

Construct a `UnitCell` for an `OrthorhombicLattice`, `TetragonalLattice` or
`CubicLattice`.
"""
function UnitCell(latticetype::LatticeType, centering::CenteringType, a::Real, b::Real, c::Real)
    if latticetype in (OrthorhombicLattice, TetragonalLattice, CubicLattice)
        UnitCell(latticetype, centering, a, b, c, 90, 90, 90)
    else
        throw(ArgumentError("More parameters needed for this type of lattice"))
    end
end


"""
    UnitCell(CubicLattice, centering, a)

Construct a `UnitCell` for a `CubicLattice`.
"""
function UnitCell(latticetype::LatticeType, centering::CenteringType, a::Real)
    if latticetype == CubicLattice
        UnitCell(latticetype, centering, a, a, a, 90, 90, 90)
    else
        throw(ArgumentError("More parameters needed for this type of lattice"))
    end
end


"""
    UnitCell(RhombohedralLattice, a, α)

Construct a `UnitCell` for a `RhombohedralLattice`.
"""
function UnitCell(latticetype::LatticeType, a::Real, α::Real)
    if latticetype == RhombohedralLattice
        UnitCell(latticetype, RhombohedralCell, a, a, a, α, α, α)
    else
        throw(ArgumentError("More parameters needed for this type of lattice"))
    end
end


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


# Returns the direct-space basis vectors as a Julia matrix
# See matrix-notation.pdf for information.  This returns an "M-matrix".
function directcartesianmatrix(uc)
    ax = Ref{Cdouble}(0)
    ay = Ref{Cdouble}(0)
    az = Ref{Cdouble}(0)
    bx = Ref{Cdouble}(0)
    by = Ref{Cdouble}(0)
    bz = Ref{Cdouble}(0)
    cx = Ref{Cdouble}(0)
    cy = Ref{Cdouble}(0)
    cz = Ref{Cdouble}(0)
    out = @ccall libcrystfel.cell_get_cartesian(uc.internalptr::Ptr{InternalUnitCell},
                                                ax::Ref{Cdouble}, ay::Ref{Cdouble}, az::Ref{Cdouble},
                                                bx::Ref{Cdouble}, by::Ref{Cdouble}, bz::Ref{Cdouble},
                                                cx::Ref{Cdouble}, cy::Ref{Cdouble}, cz::Ref{Cdouble})::Cint
    if out != 0
        throw(ErrorException("Failed to convert cell parameters"))
    end
    return [ax[] bx[] cx[]; ay[] by[] cy[]; az[] bz[] cz[]]
end


# Returns the reciprocal-space basis vectors as a Julia matrix
# See matrix-notation.pdf for information.  This returns an "R-matrix".
function reciprocalcartesianmatrix(uc)
    ax = Ref{Cdouble}(0)
    ay = Ref{Cdouble}(0)
    az = Ref{Cdouble}(0)
    bx = Ref{Cdouble}(0)
    by = Ref{Cdouble}(0)
    bz = Ref{Cdouble}(0)
    cx = Ref{Cdouble}(0)
    cy = Ref{Cdouble}(0)
    cz = Ref{Cdouble}(0)
    out = @ccall libcrystfel.cell_get_reciprocal(uc.internalptr::Ptr{InternalUnitCell},
                                                 ax::Ref{Cdouble}, ay::Ref{Cdouble}, az::Ref{Cdouble},
                                                 bx::Ref{Cdouble}, by::Ref{Cdouble}, bz::Ref{Cdouble},
                                                 cx::Ref{Cdouble}, cy::Ref{Cdouble}, cz::Ref{Cdouble})::Cint
    if out != 0
        throw(ErrorException("Failed to convert cell parameters"))
    end
    return [ax[] ay[] az[]; bx[] by[] bz[]; cx[] cy[] cz[]]
end


function Base.propertynames(uc::UnitCell; private=false)
    (:a, :b, :c, :α, :β, :γ, :latticetype, :cellparams,
     :directcartesian, :reciprocalcartesian,
     :internalptr)
end


function Base.getproperty(uc::UnitCell, name::Symbol)
    if name === :internalptr
        getfield(uc, :internalptr)
    elseif name === :cellparams
        return getcellparams(uc)
    elseif name === :latticetype
        return getlatticetype(uc)
    elseif name === :a
        return getcellparams(uc).a
    elseif name === :b
        return getcellparams(uc).b
    elseif name === :c
        return getcellparams(uc).c
    elseif name === :α
        return getcellparams(uc).α
    elseif name === :β
        return getcellparams(uc).β
    elseif name === :γ
        return getcellparams(uc).γ
    elseif name === :directcartesian
        return directcartesianmatrix(uc)
    elseif name === :reciprocalcartesian
        return reciprocalcartesianmatrix(uc)
    end
end


function Base.show(io::IO, uc::UnitCell)
    write(io, "UnitCell(")
    lt,cen,ua = uc.latticetype
    show(io, lt); write(io, ", ")
    show(io, cen); write(io, ", ")
    show(io, ua); write(io, ",\n         ")
    @printf(io, "%.3f Å, %.3f Å, %.3f Å, %.3f°, %.3f°, %.3f°",
            uc.a*1e10, uc.b*1e10, uc.c*1e10,
            rad2deg(uc.α), rad2deg(uc.β), rad2deg(uc.γ))
    write(io, ")")
end


# This type is for talking to libcrystfel, and is named to avoid conflicting
# with other quaternion libraries.
mutable struct CrystFELQuaternion
    w::Cdouble
    x::Cdouble
    y::Cdouble
    z::Cdouble
end

function randomquat()
    r = ()->2.0*rand(Float64)-1.0
    q = [r(), r(), r(), r()]
    q ./ √(sum(q.^2))
end


"""
    rotatecell(uc::UnitCell, quaternion)

Rotate a unit cell according to a quaternion (represented as a vector of 4 floats).
"""
function rotatecell(uc, quat)
    q = CrystFELQuaternion(quat...)
    out = @ccall libcrystfel.cell_rotate(uc.internalptr::Ptr{InternalUnitCell},
                                         q::CrystFELQuaternion)::Ptr{InternalUnitCell}
    if out == C_NULL
        throw(ErrorException("Failed to rotate unit cell"))
    end
    UnitCell(out)
end


"""
    rotatecell(uc::UnitCell)

Rotate a unit cell at random in three dimensions.  Use this routine for
simulating serial crystallography datasets.

Equivalent to CrystFEL routine `cell_rotate(uc, random_quaternion(<rng>))`.
"""
rotatecell(uc) = rotatecell(uc, randomquat())


end   # of module
