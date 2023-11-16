module DataTemplates

import ..CrystFEL: libcrystfel
export DataTemplate, InternalDataTemplate, loaddatatemplate
export wavelength, cameralength

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


"""
    wavelength(datatemplate)

Retrieves the radiation wavelength from a `DataTemplate`, if possible.

It's not always possible to do this.  Some geometry files declare the
wavelength as a parameter to be retrieved fom each image's metadata.  A
program using this routine should take this possibility into account.

Corresponds to CrystFEL C API function `data_template_get_wavelength_if_possible`.
"""
function wavelength(dtempl::DataTemplate)
    wl = ccall((:data_template_get_wavelength_if_possible, libcrystfel),
               Cdouble, (Ptr{InternalDataTemplate},), dtempl.internalptr)
    if isnan(wl)
        return nothing
    else
        return wl
    end
end


"""
    cameralength(datatemplate)

Retrieves the camera length from a `DataTemplate`, if possible.

It's not always possible to do this.  Some geometry files declare the
detector position(s) as a parameter to be retrieved fom each image's metadata.  A
program using this routine should take this possibility into account.

Corresponds to CrystFEL C API function `data_template_get_clen_if_possible`.
"""
function cameralength(dtempl::DataTemplate)
    clen = ccall((:data_template_get_clen_if_possible, libcrystfel),
                 Cdouble, (Ptr{InternalDataTemplate},), dtempl.internalptr)
    if isnan(clen)
        return nothing
    else
        return clen
    end
end


end  # of module
