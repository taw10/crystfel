module RefLists

import ..CrystFEL: libcrystfel
import ..CrystFEL.Symmetry: SymOpList, InternalSymOpList

export RefList, loadreflist

mutable struct InternalRefList end

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

end  # of module
