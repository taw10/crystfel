module Indexing

import ..CrystFEL: libcrystfel
import ..CrystFEL.UnitCells: UnitCell, InternalUnitCell
import ..CrystFEL.Images: Image, InternalImage
import ..CrystFEL.DataTemplates: wavelength, cameralength
export Indexer, index

mutable struct IndexingPriv end

mutable struct Indexer
    indexingpriv::Ptr{IndexingPriv}
end


function indexflags(retry, multilattice, refine, peakcheck, cellcheck)
    flags = 0
    if retry
        flags |= 1
    end
    if multilattice
        flags |= 2
    end
    if refine
        flags |= 4
    end
    if peakcheck
        flags |= 32
    end
    if cellcheck
        flags |= 64
    end
    return flags
end


function Indexer(methods, dtempl, cell; tolerances=(0.05,0.05,0.05,1.5,1.5,1.5),
        retry=true, multilattice=false, refine=true, peakcheck=true, cellcheck=true,
        wavelength_estimate=nothing, clen_estimate=missing, n_threads=1)

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

    flags = indexflags(retry, multilattice, refine, peakcheck, cellcheck)

    let wlfromdtempl = wavelength(dtempl)
        if !isnothing(wlfromdtempl)
            wavelength_estimate = wlfromdtempl
        else
            if isnothing(wavelength_estimate)
                throw(ArgumentError("Wavelength cannot be determined from data template.  "*
                                    "Use Indexer(wavelength_estimate=...)"))
            end
        end
    end

    let clenfromdtempl = cameralength(dtempl)
        if !isnothing(clenfromdtempl)
            clen_estimate = clenfromdtempl
        else
            if isnothing(clen_estimate)
                throw(ArgumentError("Camera length cannot be determined from data template.  "*
                                    "Use Indexer(clen_estimate=...)"))
            end
        end
    end

    out = @ccall libcrystfel.setup_indexing(methods::Cstring,
                                            cell.internalptr::Ptr{InternalUnitCell},
                                            tolerances::Ref{NTuple{6,Cdouble}},
                                            flags::Cint,
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


function index(image::Image, idxr::Indexer)
    @ccall libcrystfel.index_pattern_4(image.internalptr::Ptr{InternalImage},
                                       idxr.indexingpriv::Ptr{IndexingPriv},
                                       C_NULL::Ptr{Cvoid},
                                       C_NULL::Ptr{Cvoid},
                                       C_NULL::Ptr{Cvoid})::Cvoid
end


end   # of module
