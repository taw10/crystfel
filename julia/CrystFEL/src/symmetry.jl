module Symmetry

export SymOpList, InternalSymOpList, InternalIntegerMatrix

# Types for pointers returned by libcrystfel
mutable struct InternalSymOpList end
mutable struct InternalIntegerMatrix end

# Exposed types
mutable struct SymOpList
    internalptr::Ptr{InternalSymOpList}
end

mutable struct SymOp
    internalptr::Ptr{InternalIntegerMatrix}
end


function SymOpList(pointgroup::AbstractString)

    out = ccall((:get_pointgroup, :libcrystfel),
                Ptr{InternalSymOpList}, (Cstring,), pointgroup)
    if out == C_NULL
        throw(OutOfMemoryError())
    end

    sym = SymOpList(out)

    finalizer(sym) do x
        ccall((:free_symoplist, :libcrystfel),
              Cvoid, (Ptr{InternalSymOpList},), x.internalptr)
    end

    return sym

end


function Base.getindex(sym::SymOpList, i::Int)

    if i > length(sym)
        throw(BoundsError())
    end

    out = ccall((:get_symop, :libcrystfel),
                Ptr{InternalIntegerMatrix},
                (Ptr{InternalSymOpList},Ptr{Cvoid},Cint),
                sym.internalptr,C_NULL,i-1)

    if out == C_NULL
        throw(OutOfMemoryError())
    end

    return SymOp(out)

end


function Base.length(sym::SymOpList)
    return ccall((:num_equivs, :libcrystfel),
                 Cint, (Ptr{InternalSymOpList},Ptr{Cvoid}),
                 sym.internalptr, C_NULL)
end


function hkl_op(op::SymOp)
    s = ccall((:name_equiv, :libcrystfel),
               Cstring,
               (Ptr{InternalIntegerMatrix},),
               op.internalptr)
    return unsafe_string(s)
end


function symmetry_name(sym::SymOpList)
    s = ccall((:symmetry_name, :libcrystfel),
               Cstring,
               (Ptr{InternalSymOpList},),
               sym.internalptr)
    return unsafe_string(s)
end


function Base.iterate(sym::SymOpList)
    return (sym[1], 2)
end


function Base.iterate(sym::SymOpList, i)
    if i > length(sym)
        return nothing
    else
        return (sym[i], i+1)
    end
end


function Base.show(io::IO, sym::SymOpList)
    println(io, length(sym), "-element SymOpList (\"", symmetry_name(sym), "\")")
    for op in sym
        println(io, hkl_op(op))
    end
end


function Base.show(io::IO, op::SymOp)
    write(io, "SymOp(")
    write(io, "\"")
    write(io, hkl_op(op))
    write(io, "\")")
end

end  # of module
