module RefLists

import ..CrystFEL: libcrystfel
import ..CrystFEL.Symmetry: SymOpList, InternalSymOpList, symmetry_name
export RefList, loadreflist


# The internal libcrystfel structures, not exposed directly
# We only ever have e.g. a Ptr{InternalRefList}, never a real InternalRefList
mutable struct InternalRefList end
mutable struct InternalReflection end
mutable struct InternalRefListIterator end

# The Julian exposed types
mutable struct RefList
    internalptr::Ptr{InternalRefList}
    symmetry::SymOpList
end

mutable struct Reflection
    internalptr::Ptr{InternalReflection}
end

mutable struct RefListIterator
    lastrefl::Ptr{InternalReflection}
    internalptr::Ptr{InternalRefListIterator}
end


function Base.iterate(a::RefList)

    ir = Ref{Ptr{InternalRefListIterator}}()
    refl = ccall((:first_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalRefList},Ref{Ptr{InternalRefListIterator}}),
                 a.internalptr, ir)

    if refl == C_NULL
        throw(ArgumentError("Failed to find first reflection in list"))
    end

    return Reflection(refl),RefListIterator(refl, ir[])

end


function Base.iterate(a::RefList, iter)

    refl = ccall((:next_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalReflection},Ptr{InternalRefListIterator}),
                 iter.lastrefl, iter.internalptr)

    if refl == C_NULL
        return nothing
    end

    return Reflection(refl),RefListIterator(refl, iter.internalptr)

end


function loadreflist(filename::AbstractString)

    psym = Ref{Cstring}()
    out = ccall((:read_reflections_2, libcrystfel),
                Ptr{InternalRefList}, (Cstring,Ref{Cstring}),
                filename, psym)
    if out == C_NULL
        throw(ArgumentError("Failed to load reflection list"))
    end

    symmetryname = unsafe_string(psym[])
    return RefList(out, SymOpList(symmetryname))

end


function Base.show(io::IO, ::MIME"text/plain", reflist::RefList)
    write(io, "Reflection list in point group ", symmetry_name(reflist.symmetry))
end


end  # of module
