module Images

import ..CrystFEL: libcrystfel
import ..CrystFEL.DataTemplates: DataTemplate, InternalDataTemplate
import ..CrystFEL.DetGeoms: DetGeom
export Image

const HEADER_CACHE_SIZE = 128

mutable struct InternalImage
    dp::Ptr{Ptr{Cfloat}}
    bad::Ptr{Ptr{Cint}}
    sat::Ptr{Ptr{Cfloat}}
    hit::Cint
    crystals::Ptr{Ptr{Cvoid}}
    n_crystals::Cint
    indexed_by::Cint
    n_indexing_tries::Cint
    detgeom::Ptr{DetGeom}
    data_source_type::Cint
    filename::Cstring
    ev::Cstring
    data_block::Ptr{Cvoid}
    data_block_size::Csize_t
    meta_data::Cstring
    header_cache::NTuple{HEADER_CACHE_SIZE, Ptr{Cvoid}}
    n_cached_headers::Cint
    id::Cint
    serial::Cint
    spectrum::Ptr{Cvoid}
    lambda::Cdouble
    div::Cdouble
    bw::Cdouble
    peak_resolution::Cdouble
    features::Ptr{Cvoid}
    ida::Ptr{Cvoid}
end

mutable struct Image
    internalptr::Ptr{InternalImage}
end

function Base.getproperty(image::Image, name::Symbol)
    if name === :internalptr
        getfield(image, :internalptr)
    else
        idata = unsafe_load(image.internalptr)
        getproperty(idata, name)
    end
end

function Base.propertynames(image::Image)
    (:lambda, :internalptr)
end

function Image(dtempl::DataTemplate)

    out = ccall((:image_create_for_simulation, libcrystfel),
                Ptr{Image}, (Ref{InternalDataTemplate},), dtempl.internalptr)
    if out == C_NULL
        throw(OutOfMemoryError())
    end

    image = Image(out)

    finalizer(image) do x
        ccall((:image_free, libcrystfel), Cvoid, (Ptr{InternalImage},), x.internalptr)
    end

    return image
end

end  # of module
