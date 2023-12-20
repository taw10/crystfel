module Streams

import ..CrystFEL: libcrystfel
import ..CrystFEL.DataTemplates: DataTemplate, InternalDataTemplate
export Stream

# Represents the real C-side (opaque) structure.
mutable struct InternalStream end

# The Julia-side structure, needed to house the pointer to the C structure
mutable struct Stream
    internalptr::Ptr{InternalStream}
end


function Stream(filename, mode::AbstractString, dtempl::DataTemplate)
    if mode == "w"
        out = @ccall libcrystfel.stream_open_for_write(filename::Cstring,
                                                       dtempl.internalptr::Ptr{InternalDataTemplate})::Ptr{InternalStream}
        if out == C_NULL
            throw(ErrorException("Failed to open stream for reading"))
        end
        return Stream(out)
    elseif mode =="r"
        throw(ArgumentError("To open a stream for reading, don't provide the DataTemplate"))
    else
        throw(ArgumentError("Unrecognised CrystFEL stream mode"*mode))
    end
end


function Stream(filename, mode::AbstractString)
    if mode == "r"
        out = @ccall libcrystfel.stream_open_for_read(filename::Cstring)::Ptr{InternalStream}
        if out == C_NULL
            throw(ErrorException("Failed to open stream for reading"))
        end
        return Stream(out)
    elseif mode == "w"
        throw(ArgumentError("To open a stream for writing, you must provide "
                            *"a DataTemplate: use Stream(filename, \"w\", dtempl)"))
    else
        throw(ArgumentError("Unrecognised CrystFEL stream mode"*mode))
    end
end


end  # of module
