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


"""
    Indexer(methods, dtempl, cell;
            tolerances=(0.05,0.05,0.05,1.5,1.5,1.5),
            retry=true, multilattice=false, refine=true,
            peakcheck=true, cellcheck=true,
            wavelength_estimate=nothing, clen_estimate=missing,
            n_threads=1)

Creates a new CrystFEL indexing engine, which you can later apply to CrystFEL
`Image` structures using `index()`.
"""
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

    tols = Vector{Cfloat}(undef, 6)
    tols[1] = tolerances[1]
    tols[2] = tolerances[2]
    tols[3] = tolerances[3]
    tols[4] = deg2rad(tolerances[4])
    tols[5] = deg2rad(tolerances[5])
    tols[6] = deg2rad(tolerances[6])

    out = @ccall libcrystfel.setup_indexing(methods::Cstring,
                                            cell.internalptr::Ptr{InternalUnitCell},
                                            tols::Ptr{Cfloat},
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


"""
    index(image::Image, indexer::Indexer; mille=nothing)

Index `image` using `indexer`.

If `mille` is a valid `Mille` object, detector geometry alignment data will
be written.
"""
function index(image::Image, idxr::Indexer; mille=nothing)

    if mille === nothing
        imille = C_NULL
    else
        imille = mille.internalptr
    end

    @ccall libcrystfel.index_pattern_4(image.internalptr::Ptr{InternalImage},
                                       idxr.indexingpriv::Ptr{IndexingPriv},
                                       C_NULL::Ptr{Cvoid},
                                       C_NULL::Ptr{Cvoid},
                                       imille::Ptr{Cvoid},
                                       99::Cint)::Cvoid
end


end   # of module
