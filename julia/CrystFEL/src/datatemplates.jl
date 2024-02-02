module DataTemplates

import ..CrystFEL: libcrystfel
export DataTemplate, InternalDataTemplate, loaddatatemplate
export wavelength, cameralength
export translategroup!, rotategroup!

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


"""
    translategroup!(datatemplate, groupname, xshift, yshift, zshift)

Modifies `DataTemplate` by moving the specified panel group by the specified
amount (in metres).

Corresponds to CrystFEL C API function `data_template_translate_group`.
"""
function translategroup!(dtempl::DataTemplate, groupname, xshift, yshift, zshift)
    r = @ccall libcrystfel.data_template_translate_group_m(dtempl.internalptr::Ptr{InternalDataTemplate},
                                                             groupname::Cstring,
                                                             xshift::Cdouble,
                                                             yshift::Cdouble,
                                                             zshift::Cdouble)::Cint
    if r != 0
        throw(ErrorException("Failed to shift DataTemplate"))
    end

end


"""
    rotategroup!(datatemplate, groupname, angle, axis)

Modifies `DataTemplate` by rotating the specified panel group by the specified
amount (in degrees) about the specified xaxis (:x, :y or :z).

Corresponds to CrystFEL C API function `data_template_rotate_group`.
"""
function rotategroup!(dtempl::DataTemplate, groupname, angle, axis)
    r = @ccall libcrystfel.data_template_rotate_group(dtempl.internalptr::Ptr{InternalDataTemplate},
                                                      groupname::Cstring,
                                                      deg2rad(angle)::Cdouble,
                                                      String(axis)[1]::Cchar)::Cint
    if r != 0
        throw(ErrorException("Failed to rotate DataTemplate"))
    end

end

end  # of module
