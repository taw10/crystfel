module Indexing

import ..CrystFEL: libcrystfel
import ..CrystFEL.UnitCells: UnitCell, InternalUnitCell
export Indexer

mutable struct IndexingPriv end

mutable struct Indexer
    indexingpriv::Ptr{IndexingPriv}
end


function Indexer(methods, cell; tolerances=(0.05,0.05,0.05,1.5,1.5,1.5),
        retry=true, multilattice=false, refine=true, peakcheck=true, cellcheck=true,
        wavelength_estimate=0.0, clen_estimate=0.0, n_threads=1)

    taketwoopts = Ref{Ptr{Cvoid}}(C_NULL)
    xgandalfopts = Ref{Ptr{Cvoid}}(C_NULL)
    pinkindexeropts = Ref{Ptr{Cvoid}}(C_NULL)
    felixopts = Ref{Ptr{Cvoid}}(C_NULL)
    fromfileopts = Ref{Ptr{Cvoid}}(C_NULL)
    asdfopts = Ref{Ptr{Cvoid}}(C_NULL)

    @ccall libcrystfel.default_method_options(taketwoopts::Ref{Ptr{Cvoid}},
                                              xgandalfopts::Ref{Ptr{Cvoid}},
                                              pinkindexeropts::Ref{Ptr{Cvoid}},
                                              felixopts::Ref{Ptr{Cvoid}},
                                              fromfileopts::Ref{Ptr{Cvoid}},
                                              asdfopts::Ref{Ptr{Cvoid}})::Cvoid

    out = @ccall libcrystfel.setup_indexing(methods::Cstring,
                                            cell.internalptr::Ptr{InternalUnitCell},
                                            tolerances::Ref{NTuple{6,Cdouble}},
                                            0::Cint,
                                            wavelength_estimate::Cdouble,
                                            clen_estimate::Cdouble,
                                            n_threads::Cint,
                                            taketwoopts[]::Ptr{Cvoid},
                                            xgandalfopts[]::Ptr{Cvoid},
                                            pinkindexeropts[]::Ptr{Cvoid},
                                            felixopts[]::Ptr{Cvoid},
                                            fromfileopts[]::Ptr{Cvoid},
                                            asdfopts[]::Ptr{Cvoid})::Ptr{IndexingPriv}

    if out == C_NULL
        throw(ErrorException("Indexing setup failed"))
    end

    indexer = Indexer(out)

    finalizer(indexer) do x
        @ccall libcrystfel.cleanup_indexing(x.indexingpriv::Ptr{IndexingPriv})::Cvoid
    end

    return indexer

end


end   # of module
