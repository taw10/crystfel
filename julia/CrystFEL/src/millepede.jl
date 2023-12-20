module Millepede

import ..CrystFEL: libcrystfel
export Mille

mutable struct InternalMille end

mutable struct Mille
    internalptr::Ptr{InternalMille}
end


"""
    Mille(filename)

Creates a new "Mille" object, which can be passed to `CrystFEL.Indexing.index()`
to accumulate detector geometry alignment information in the specified file.

When you've finished adding data, call `close()` on the Mille object.  This will
be done automatically when the object is freed by the garbage collector, but
that might not happen until Julia exits.

Corresponds to CrystFEL C API routine `crystfel_mille_new`.
"""
function Mille(filename::AbstractString)
    out = @ccall libcrystfel.crystfel_mille_new(filename::Cstring)::Ptr{InternalMille}
    if out == C_NULL
        throw(ArgumentError("Failed to create Millepede object"))
    end
    finalizer(close, Mille(out))
end


function Base.close(x::Mille)
    if x.internalptr != C_NULL
        @ccall libcrystfel.crystfel_mille_free(x.internalptr::Ptr{InternalMille})::Cvoid
        x.internalptr = C_NULL
    end
end


end  # of module
