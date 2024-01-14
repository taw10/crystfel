module PeakLists

using Printf
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

mutable struct PeakList
    internalptr::Ptr{InternalPeakList}
end


function PeakList()
    out = @ccall libcrystfel.image_feature_list_new()::Ptr{InternalPeakList}
    if out == C_NULL
        throw(ArgumentError("Failed to create peak list"))
    end
    finalizer(PeakList(out)) do pl
        @ccall libcrystfel.image_feature_list_free(pl.internalptr::Ptr{InternalPeakList})::Cvoid
    end
end


function Base.deepcopy(peaklist::PeakList)
    out = @ccall libcrystfel.image_feature_list_copy(peaklist.internalptr::Ptr{InternalPeakList})::Ptr{InternalPeakList}
    if out == C_NULL
        throw(ErrorException("Failed to copy peak list"))
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
    if length(peaklist) > 0
        return peaklist[1],(1,length(peaklist))
    else
        return nothing
    end
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


function Base.show(io::IO, ::MIME"text/plain", peaklist::PeakList)
    println(io, "Peak list with ", length(peaklist), " peaks")
    print(io, "     fs      ss  panel  intensity  name")
    let n = 0
        for pk in Iterators.take(peaklist, 11)
            if n == 10
                # We have printed 10 already, and are here again.  Truncate...
                print(io, "\n      ⋮       ⋮      ⋮          ⋮  ⋮")
                break
            end
            write(io, "\n")
            @printf(io, "%7.2f %7.2f %6i %10.2f  %s",
                    pk.fs, pk.ss, pk.panelnumber, pk.intensity, pk.name)
            n += 1
        end
    end
end

end   # of module
