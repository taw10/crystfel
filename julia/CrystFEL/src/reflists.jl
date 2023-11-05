module RefLists

import ..CrystFEL: libcrystfel
import ..CrystFEL.Symmetry: SymOpList, InternalSymOpList, symmetry_name
export RefList, loadreflist, reflections


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
    reflist::RefList
    lastrefl::Ptr{InternalReflection}
    internalptr::Ptr{InternalRefListIterator}
end


function reflections(reflist::RefList)
    iter = RefListIterator(reflist, C_NULL, C_NULL)
    finalizer(iter) do x
        @async println("Finalising iterator: ", x)
    end
    return iter
end


function Base.iterate(iter::RefListIterator)

    rli = Ref{Ptr{InternalRefListIterator}}(C_NULL)
    refl = ccall((:first_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalRefList},Ref{Ptr{InternalRefListIterator}}),
                 iter.reflist.internalptr, rli)

    if refl == C_NULL
        throw(ArgumentError("Failed to find first reflection in list"))
    end

    iter.lastrefl = refl
    iter.internalptr = rli[]

    return Reflection(refl),iter

end


function Base.iterate(iter::RefListIterator, _)

    refl = ccall((:next_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalReflection},Ptr{InternalRefListIterator}),
                 iter.lastrefl, iter.internalptr)

    if refl == C_NULL
        return nothing
    end

    iter.lastrefl = refl

    return Reflection(refl),iter

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


function intensity(refl::Reflection)
    ccall((:get_intensity, libcrystfel),
          Cdouble, (Ptr{InternalReflection},), refl.internalptr)
end

function Base.show(io::IO, refl::Reflection)
    write(io, "Reflection(intensity=")
    show(io, intensity(refl))
    write(io, ")")
end

end  # of module
