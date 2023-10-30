module DataTemplates

import ..CrystFEL: libcrystfel
export DataTemplate, InternalDataTemplate, loaddatatemplate

# Represents the real C-side (opaque) structure.
mutable struct InternalDataTemplate end

# The Julia-side structure, needed to house the pointer to the C structure
# Without this, we would only ever have a Ptr{DataTemplate}, not a DataTemplate.
mutable struct DataTemplate
    internalptr::Ptr{InternalDataTemplate}
end

"""
    loaddatatemplate(filename)

Creates a CrystFEL DataTemplate by loading a geometry file.

Corresponds to CrystFEL C API function `data_template_new_from_file()`.
"""
function loaddatatemplate(filename::AbstractString)

    out = ccall((:data_template_new_from_file, libcrystfel),
                Ptr{InternalDataTemplate}, (Cstring,), filename)
    if out == C_NULL
        throw(ArgumentError("Failed to load geometry file"))
    end

    dt = DataTemplate(out)

    finalizer(dt) do x
        ccall((:data_template_free, libcrystfel),
              Cvoid, (Ptr{InternalDataTemplate},), x.internalptr)
    end

    return dt
end

end  # of module
