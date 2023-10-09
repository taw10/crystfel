module RefLists

using ..Symmetry

export RefList, Reflection, loadreflist, savereflections

mutable struct InternalRefList end

mutable struct RefList
    internalptr::Ptr{InternalRefList}
    symmetry::SymOpList
end

function loadreflist(filename::AbstractString)

    psym = Ptr{InternalSymOpList}()
    out = ccall((:read_reflections_2, :libcrystfel),
                Ptr{InternalRefList}, (Cstring,Ptr{InternalSymOpList}),
                filename, psym)
    if out == C_NULL
        throw(OutOfMemoryError())
    end

    return RefList(out, SymOpList(psym))

end

end  # of module
