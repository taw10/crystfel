module PeakLists

import ..CrystFEL: libcrystfel
export PeakList, InternalPeakList

mutable struct InternalPeak
    fs::Cdouble
    ss::Cdouble
    panelnumber::Cint
    intensity::Cdouble
    name::Cstring
end

mutable struct InternalPeakList end

struct PeakList
    internalptr::Ptr{InternalPeakList}
end


function PeakList()
    out = @ccall libcrystfel.image_feature_list_new()::Ptr{InternalPeakList}
    if out == C_NULL
        throw(ArgumentError("Failed to create peak list"))
    end
    PeakList(out)
end


function Base.length(peaklist::PeakList)
    @ccall libcrystfel.image_feature_count(peaklist.internalptr::Ptr{InternalPeakList})::Cint
end

Base.firstindex(peaklist::PeakList) = 1
Base.lastindex(peaklist::PeakList) = length(peaklist)


function Base.push!(peaklist::PeakList, fs, ss, panelnumber, intensity, name=nothing)

    rname = isnothing(name) ? C_NULL : @ccall strdup(name::Cstring)::Cstring

    @ccall libcrystfel.image_add_feature(peaklist.internalptr::Ptr{InternalPeakList},
                                         fs::Cdouble, ss::Cdouble, panelnumber::Cint,
                                         intensity::Cdouble, rname::Cstring)::Cvoid
end


function Base.getindex(peaklist::PeakList, n)
    out = @ccall(libcrystfel.image_get_feature(peaklist.internalptr::Ptr{InternalPeakList},
                                               (n-1)::Cint)::Ptr{InternalPeak})
    if out == C_NULL
        throw(BoundsError(peaklist, n))
    end
    pdata = unsafe_load(out)
    if pdata.name == C_NULL
        nname = nothing
    else
        nname = unsafe_string(pdata.name)
    end
    return (fs=pdata.fs, ss=pdata.ss, panelnumber=pdata.panelnumber,
            intensity=pdata.intensity, name=nname)
end


function Base.iterate(peaklist::PeakList)
    return peaklist[1],(1,length(peaklist))
end


function Base.iterate(peaklist::PeakList, state)
    let nxt = state[1]+1
        len = state[2]
        if nxt == len+1
            return nothing
        else
            return peaklist[nxt],(nxt,state[2])
        end
    end
end


end   # of module
