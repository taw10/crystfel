module RefLists

import ..CrystFEL: libcrystfel
import ..CrystFEL.Symmetry: SymOpList, InternalSymOpList, symmetry_name
export RefList, loadreflist


# The internal libcrystfel structure, not exposed directly
mutable struct InternalRefList end


# The Julian exposed type
mutable struct RefList
    internalptr::Ptr{InternalRefList}
    symmetry::SymOpList
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
