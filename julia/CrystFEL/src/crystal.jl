module Crystals

import ..CrystFEL: libcrystfel
import ..CrystFEL.RefLists: RefList, InternalRefList, UnmergedReflection
import ..CrystFEL.UnitCells: UnitCell, InternalUnitCell
export Crystal, InternalCrystal

# Represents the real C-side (opaque) structure.
mutable struct InternalCrystal end

mutable struct Crystal
    internalptr::Ptr{InternalCrystal}
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

    cr = Crystal(out)

    finalizer(cr) do x
        ccall((:crystal_free, libcrystfel), Cvoid, (Ptr{InternalCrystal},),
              x.internalptr)
    end

    return cr
end


function Base.setproperty!(cr::Crystal, name::Symbol, val)
    if name === :internalptr
        setfield!(cr, :internalptr, val)
    else
        if name === :reflections
            if val isa RefList{UnmergedReflection}
                ccall((:crystal_set_reflections, libcrystfel),
                      Cvoid, (Ptr{InternalCrystal},Ptr{InternalRefList}),
                      cr.internalptr, val.internalptr)
            else
                throw(ArgumentError("Must be a RefList{UnmergedReflection}"))
            end
        end
    end
end


end   # of module
