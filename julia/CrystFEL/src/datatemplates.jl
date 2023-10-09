module DataTemplates

export DataTemplate, loaddatatemplate

mutable struct InternalDataTemplate end

mutable struct DataTemplate
    internalptr::Ptr{InternalDataTemplate}
end

function loaddatatemplate(filename::AbstractString)

    out = ccall((:data_template_new_from_file, :libcrystfel),
                Ptr{InternalDataTemplate}, (Cstring,), filename)
    if out == C_NULL
        throw(OutOfMemoryError())
    end

    dt = DataTemplate(out)

    finalizer(dt) do x
        ccall((:data_template_free, :libcrystfel),
              Cvoid, (Ptr{InternalDataTemplate},), x.internalptr)
    end

    return dt
end

end  # of module
