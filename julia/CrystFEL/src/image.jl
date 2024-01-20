module Images

using Printf

import ..CrystFEL: libcrystfel
import ..CrystFEL.DataTemplates: DataTemplate, InternalDataTemplate
import ..CrystFEL.DetGeoms: DetGeom
import ..CrystFEL.PeakLists: PeakList, InternalPeakList
import ..CrystFEL.Crystals: Crystal, InternalCrystal
import ..CrystFEL.RefLists: RefList, InternalRefList, UnmergedReflection
import ..CrystFEL.Symmetry: SymOpList
export Image

const HEADER_CACHE_SIZE = 128

mutable struct CrystalRefListPair
    crystal::Ptr{InternalCrystal}
    reflist::Ptr{InternalRefList}
    owns_crystal::Cint
    owns_reflist::Cint
end

mutable struct InternalImage
    dp::Ptr{Ptr{Cfloat}}
    bad::Ptr{Ptr{Cint}}
    sat::Ptr{Ptr{Cfloat}}
    hit::Cint
    crystals::Ptr{CrystalRefListPair}
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
    peaklist::Ptr{InternalPeakList}
    ida::Ptr{Cvoid}
    owns_peaklist::Cint
end


mutable struct Image
    internalptr::Ptr{InternalImage}
    peaklist::Union{Nothing,PeakList}
    crystals
    reflists
end


function makecrystallist(image, listptr, n)

    crystals = []

    if listptr == C_NULL
        return crystals
    end

    for i in 1:n
        pairptr = unsafe_load(listptr, i)

        # Re-use old Crystal if possible
        n = findfirst(getfield(image, :crystals)) do x
            x.internalptr == pairptr.crystal
        end

        if n !== nothing
            cr = getfield(image, :crystals)[n]
        else
            cr = Crystal(pairptr.crystal, nothing)
        end

        if pairptr.reflist == C_NULL
            reflist = nothing
        else
            reflist = RefList{UnmergedReflection}(pairptr.reflist, SymOpList("1"))
            pairptr.owns_reflist = 0
        end
        push!(crystals, (crystal=cr, reflections=reflist))
        pairptr.owns_crystal = 0
        unsafe_store!(listptr, pairptr, i)
        # We are now responsible for freeing the Crystal and RefList
    end

    image.crystals = map(x->x.crystal, crystals)
    image.reflists = map(x->x.reflections, crystals)
    return crystals

end


function getpeaklist(image)
    idata = unsafe_load(image.internalptr)
    if (getfield(image, :peaklist) === nothing) ||
        (idata.peaklist != getfield(image, :peaklist).internalptr)
        if idata.peaklist != C_NULL
            setfield!(image, :peaklist, PeakList(idata.peaklist))
            # From now on, Julia is completely responsible for freeing the peaklist
            idata.owns_peaklist = 0
            unsafe_store!(image.internalptr, idata)
        else
            setfield!(image, :peaklist, nothing)
        end
    end
    return getfield(image, :peaklist)
end


function Base.getproperty(image::Image, name::Symbol)
    if name === :internalptr
        getfield(image, :internalptr)
    elseif name === :peaklist
        getpeaklist(image)
    else
        idata = unsafe_load(image.internalptr)

        if name === :crystals
            return makecrystallist(image,
                                   getfield(idata, :crystals),
                                   getfield(idata, :n_crystals))

        else
            getfield(idata, name)
        end
    end
end


strdup(str) = @ccall strdup(str::Cstring)::Cstring


function assert_type(val, type)
    if !(val isa type)
        throw(ArgumentError("Must be a "*string(type)*" (have "*string(typeof(val))*" instead)"))
    end
end


function set_peaklist(image, new_pl)

    assert_type(new_pl, PeakList)

    idata = unsafe_load(image.internalptr)
    if (idata.owns_peaklist == 0) && (idata.peaklist != C_NULL)
        @ccall libcrystfel.image_feature_list_free(idata.peaklist::Ptr{InternalPeakList})::Cvoid
    end
    idata.peaklist = new_pl.internalptr
    idata.owns_peaklist = 0
    unsafe_store!(image.internalptr, idata)

end


function Base.setproperty!(image::Image, name::Symbol, val)
    if name === :internalptr
        setfield!(image, :internalptr, val)
    else

        if name === :peaklist
            return set_peaklist(image, val)

        elseif name === :filename
            assert_type(val, AbstractString)
            val = strdup(val)

        elseif name === :ev
            assert_type(val, AbstractString)
            val = strdup(val)

        elseif name === :crystals
            return setfield!(image, :crystals, val)

        elseif name === :reflists
            return setfield!(image, :reflists, val)

        end

        idata = unsafe_load(image.internalptr)
        setproperty!(idata, name, val)
        unsafe_store!(image.internalptr, idata)

    end
end


function Base.propertynames(image::Image; private=false)
    if private
        fieldnames(InternalImage)
    else
        tuple(fieldnames(InternalImage)..., :internalptr)
    end
end


function Base.push!(image::Image, cr::Crystal)
    ccall((:image_add_crystal, libcrystfel),
          Cvoid, (Ptr{InternalImage},Ptr{InternalCrystal}),
          image.internalptr, cr.internalptr)
end


"""
    Image(dtempl::DataTemplate)

Creates a CrystFEL image structure, not linked to any file or data block,
i.e. for simulation purposes.  This will fail if `dtempl` contains any
references to metadata fields, e.g. `photon_energy = /LCLS/photon_energy eV`.

Corresponds to CrystFEL C API function `image_create_for_simulation()`.
"""
function Image(dtempl::DataTemplate)

    out = ccall((:image_create_for_simulation, libcrystfel),
                Ptr{Image}, (Ref{InternalDataTemplate},), dtempl.internalptr)
    if out == C_NULL
        throw(ArgumentError("Failed to create image"))
    end

    image = Image(out, nothing, [], [])

    finalizer(image) do x
        ccall((:image_free, libcrystfel), Cvoid, (Ptr{InternalImage},), x.internalptr)
    end

    return image
end


"""
    Image(dtempl::DataTemplate, filename::AbstractString, event::AbstractString,
          no_image_data=false, no_mask_data=false)

Loads an image from the filesystem.

Corresponds to CrystFEL C API function `image_read()`.
"""
function Image(dtempl::DataTemplate,
               filename::AbstractString,
               event::AbstractString="//",
               no_image_data=false,
               no_mask_data=false)

    out = @ccall libcrystfel.image_read(dtempl.internalptr::Ptr{InternalDataTemplate},
                                        filename::Cstring, event::Cstring,
                                        no_image_data::Cint, no_mask_data::Cint,
                                        C_NULL::Ptr{Cvoid})::Ptr{Image}
    if out == C_NULL
        throw(ArgumentError("Failed to load image"))
    end

    image = Image(out, nothing, [], [])

    finalizer(image) do x
        ccall((:image_free, libcrystfel), Cvoid, (Ptr{InternalImage},), x.internalptr)
    end

    return image
end


function Base.show(io::IO, mime::MIME"text/plain", image::Image)

    idata = unsafe_load(image.internalptr)
    @printf(io, "CrystFEL.Image(%p):\n\n", image.internalptr)

    println(io, "                    Serial number: ", idata.serial)
    write(io, "                         Filename: ")
    if idata.filename == C_NULL
        write(io, "<not set>")
    else
        write(io, unsafe_string(idata.filename))
    end
    write(io, "\n")

    write(io, "                         Frame ID: ")
    if idata.ev == C_NULL
        write(io, "<not set>")
    else
        write(io, unsafe_string(idata.ev))
    end
    write(io, "\n")

    write(io, "\n")
    println(io, "                       Wavelength: ", idata.lambda*1e10, " Å")
    println(io, "                        Bandwidth: ", idata.bw*100, " %")
    println(io, "                       Divergence: ", idata.div*1e3, " mrad")

    write(io, "\n")
    if idata.peaklist != C_NULL
        let npk = @ccall libcrystfel.image_feature_count(idata.peaklist::Ptr{InternalPeakList})::Cint
            println(io, "                  Number of peaks: ", npk)
        end
    else
        println(io, "                  Number of peaks: 0 (no peak list)")
    end

    println(io, "             Estimated resolution: ", 1e10/idata.peak_resolution, " Å")
    write(io, "                         Hit flag: ")
    if idata.hit != 0
        write(io, "set")
    else
        write(io, "not set")
    end
    write(io, "\n")

    write(io, "\n")
    println(io, "               Number of crystals: ", idata.n_crystals)
    println(io, " Number of indexing attempts made: ", idata.n_crystals)
    println(io, "             Indexed by algorithm: ", idata.indexed_by)
end


function Base.show(io::IO, image::Image)
    @printf(io, "CrystFEL.Image(%p)", image.internalptr)
end


end  # of module
