module Streams

import ..CrystFEL: libcrystfel
import ..CrystFEL.DataTemplates: DataTemplate, InternalDataTemplate
import ..CrystFEL.Images: Image, InternalImage
export Stream, chunkwrite, chunkread, allcrystals

# Represents the real C-side (opaque) structure.
mutable struct InternalStream end

# The Julia-side structure, needed to house the pointer to the C structure
mutable struct Stream
    internalptr::Ptr{InternalStream}
end


"""
    Stream(filename, "w", dtempl)

Opens a CrystFEL stream for writing.  Note that you must provide a `DataTemplate`,
which is needed to translate "panel coordinates" to "file coordinates".

Corresponds to CrystFEL C API routine `stream_open_for_write`.
"""
function Stream(filename, mode::AbstractString, dtempl::DataTemplate)

    if mode == "w"
        out = @ccall libcrystfel.stream_open_for_write(filename::Cstring,
                               dtempl.internalptr::Ptr{InternalDataTemplate})::Ptr{InternalStream}
        if out == C_NULL
            throw(ErrorException("Failed to open stream for reading"))
        end

        @ccall libcrystfel.stream_write_data_template(out::Ptr{InternalStream},
                               dtempl.internalptr::Ptr{InternalDataTemplate})::Cvoid

        finalizer(close, Stream(out))

    elseif mode =="r"
        throw(ArgumentError("To open a stream for reading, don't provide the DataTemplate"))

    else
        throw(ArgumentError("Unrecognised CrystFEL stream mode"*mode))
    end
end


"""
    Stream(filename, "r")

Opens a CrystFEL stream for reading.

Close the stream with `close` when you've finished (this will happen
automatically when the `Stream` object is finalized).

Corresponds to CrystFEL C API routine `stream_open_for_read`.
"""
function Stream(filename, mode::AbstractString)

    if mode == "r"
        out = @ccall libcrystfel.stream_open_for_read(filename::Cstring)::Ptr{InternalStream}
        if out == C_NULL
            throw(ErrorException("Failed to open stream for reading"))
        end
        finalizer(close, Stream(out))

    elseif mode == "w"
        throw(ArgumentError("To open a stream for writing, you must provide "
                            *"a DataTemplate: use Stream(filename, \"w\", dtempl)"))

    else
        throw(ArgumentError("Unrecognised CrystFEL stream mode"*mode))
    end
end


function Base.close(st::Stream)
    if st.internalptr != C_NULL
        @ccall libcrystfel.stream_close(st.internalptr::Ptr{InternalStream})::Cvoid
        st.internalptr = C_NULL
    end
end


function streamflags(peaks, reflections, imagedata)
    flags = 0
    if reflections
        flags |= 2
    end
    if peaks
        flags |= 4
    end
    if imagedata
        flags |= 8
    end
    return flags
end


function chunkwrite(st::Stream, image::Image; peaks=true, reflections=true)
    st.internalptr == C_NULL && throw(ErrorException("Stream is closed"))
    flags = streamflags(peaks, reflections, false)
    @ccall libcrystfel.stream_write_chunk(st.internalptr::Ptr{InternalStream},
                                          image.internalptr::Ptr{InternalImage},
                                          flags::Cint)::Cvoid
end


function chunkread(st::Stream; peaks=true, reflections=true, datageom=true)

    st.internalptr == C_NULL && throw(ErrorException("Stream is closed"))

    flags = streamflags(peaks, reflections, datageom)
    out = @ccall libcrystfel.stream_read_chunk(st.internalptr::Ptr{InternalStream},
                                               flags::Cint)::Ptr{InternalImage}
    out == C_NULL && return nothing

    finalizer(Image(out, nothing, [], [])) do x
        ccall((:image_free, libcrystfel), Cvoid, (Ptr{InternalImage},), x.internalptr)
    end

end


function allcrystals(st)
    Channel() do ch
        while true
            image = chunkread(st, peaks=false, reflections=true, datageom=false)
            image === nothing && break
            for cr in image.crystals
                put!(ch, (cr.crystal, cr.reflections))
            end
        end
    end
end


end  # of module
