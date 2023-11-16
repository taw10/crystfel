module Crystals

import ..CrystFEL: libcrystfel
import ..CrystFEL.RefLists: RefList, UnmergedReflection
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

    ccall((:crystal_set_cell, libcrystfel),
          Cvoid, (Ptr{InternalCrystal},Ptr{InternalUnitCell}),
          out, cell.internalptr)

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


end   # of module
