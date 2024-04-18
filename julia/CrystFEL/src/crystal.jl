module Crystals

using Printf

import ..CrystFEL: libcrystfel
import ..CrystFEL.RefLists: RefList, InternalRefList, UnmergedReflection
import ..CrystFEL.UnitCells: UnitCell, InternalUnitCell
export Crystal, InternalCrystal

# Represents the real C-side (opaque) structure.
mutable struct InternalCrystal end

mutable struct Crystal
    internalptr::Ptr{InternalCrystal}
    cell
end


function Crystal(cell::UnitCell; profileradius=2e6, mosaicity=0)

    out = ccall((:crystal_new, libcrystfel),
                Ptr{InternalCrystal}, ())

    if out == C_NULL
        throw(ArgumentError("Failed to create crystal"))
    end

    # We make a copy of the cell, to avoid memory model shenanigans
    uccopy = ccall((:cell_new_from_cell, libcrystfel),
                   Ptr{InternalUnitCell}, (Ptr{InternalUnitCell},),
                   cell.internalptr)

    ccall((:crystal_set_cell, libcrystfel),
          Cvoid, (Ptr{InternalCrystal},Ptr{InternalUnitCell}),
          out, uccopy)

    ccall((:crystal_set_profile_radius, libcrystfel),
          Cvoid, (Ptr{InternalCrystal},Cdouble),
          out, profileradius)

    ccall((:crystal_set_mosaicity, libcrystfel),
          Cvoid, (Ptr{InternalCrystal},Cdouble),
          out, mosaicity)

    cr = Crystal(out, nothing)

    finalizer(cr) do x
        ccall((:crystal_free, libcrystfel), Cvoid, (Ptr{InternalCrystal},),
              x.internalptr)
    end

    return cr
end


function Base.setproperty!(cr::Crystal, name::Symbol, val)
    if name === :internalptr
        setfield!(cr, :internalptr, val)
    end
end


function getcell(cr)
    out = @ccall libcrystfel.crystal_relinquish_cell(cr.internalptr::Ptr{InternalCrystal})::Ptr{InternalUnitCell}
    if getfield(cr, :cell) === nothing || getfield(cr, :cell).internalptr != out
        if out != C_NULL
            setfield!(cr, :cell, UnitCell(out))
        else
            setfield!(cr, :cell, nothing)
        end
    end
    return getfield(cr, :cell)
end


function Base.getproperty(cr::Crystal, name::Symbol)
    if name === :internalptr
        getfield(cr, :internalptr)
    elseif name === :cell
        return getcell(cr)
    elseif name === :Bfac
        return @ccall libcrystfel.crystal_get_Bfac(cr.internalptr::Ptr{InternalCrystal})::Cdouble
    elseif name === :osf
        return @ccall libcrystfel.crystal_get_osf(cr.internalptr::Ptr{InternalCrystal})::Cdouble
    elseif name === :mos
        return @ccall libcrystfel.crystal_get_mosaicity(cr.internalptr::Ptr{InternalCrystal})::Cdouble
    elseif name === :r
        return @ccall libcrystfel.crystal_get_profile_radius(cr.internalptr::Ptr{InternalCrystal})::Cdouble
    elseif name === :resolution
        return @ccall libcrystfel.crystal_get_resolution_limit(cr.internalptr::Ptr{InternalCrystal})::Cdouble
    elseif name === :flag
        return @ccall libcrystfel.crystal_get_user_flag(cr.internalptr::Ptr{InternalCrystal})::Cint
    else
        throw(ErrorException("Type Crystal has no field "*String(name)))
    end
end


function Base.show(io::IO, mime::MIME"text/plain", cr::Crystal)
    @printf(io, "CrystFEL.Crystal(%p):\n\n", cr.internalptr)
    if cr.cell !== nothing
        show(io, cr.cell)
        write(io, "\n\n")
    else
        write(io, "Unit cell parameters not set\n\n")
    end
    println(io, " Linear scale factor: ", cr.osf)
    println(io, "  Debye-Walle factor: ", cr.Bfac)
    println(io, "           Mosaicity: ", cr.mos)
    println(io, "      Profile radius: ", cr.r/1e9, " nm⁻¹")
    println(io, "    Resolution limit: ", cr.resolution)
    println(io, "                Flag: ", cr.flag)
end


function Base.show(io::IO, cr::Crystal)
    @printf(io, "CrystFEL.Crystal(%p)", cr.internalptr)
end


end   # of module
