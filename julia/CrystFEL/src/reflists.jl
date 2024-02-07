module RefLists

using Printf
import ..CrystFEL: libcrystfel
import ..CrystFEL.Symmetry: SymOpList, InternalSymOpList, symmetry_name
export RefList, loadreflist, savereflist!
export Reflection, UnmergedReflection, MergedReflection
export InternalRefList


# The internal libcrystfel structures, not exposed directly
# We only ever have e.g. a Ptr{InternalRefList}, never a real InternalRefList
mutable struct InternalRefList end
mutable struct InternalReflection end
mutable struct InternalRefListIterator end

# The Julian exposed types
abstract type Reflection end

mutable struct RefList{T<:Reflection}
    internalptr::Ptr{InternalRefList}
    symmetry::SymOpList
end

mutable struct RefListIterator
    lastrefl::Ptr{InternalReflection}
    internalptr::Ptr{InternalRefListIterator}
end

mutable struct MergedReflection <: Reflection
    internalptr::Ptr{InternalReflection}
end

mutable struct UnmergedReflection <: Reflection
    internalptr::Ptr{InternalReflection}
end


function RefList{MergedReflection}(sym::SymOpList)
    out = @ccall libcrystfel.reflist_new()::Ptr{InternalRefList}
    if out == C_NULL
        throw(ErrorException("Failed to create RefList"))
    end
    finalizer(RefList{MergedReflection}(out, sym)) do x
        @ccall libcrystfel.reflist_free(x.internalptr::Ptr{InternalRefList})::Cvoid
    end
end


function Base.iterate(reflist::RefList{T}) where T

    rli = Ref{Ptr{InternalRefListIterator}}(C_NULL)
    refl = ccall((:first_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalRefList},Ref{Ptr{InternalRefListIterator}}),
                 reflist.internalptr, rli)

    if refl == C_NULL
        throw(ArgumentError("Failed to find first reflection in list"))
    end

    iter = RefListIterator(refl,rli[])
    finalizer(iter) do x
        ccall((:free_reflistiterator, libcrystfel),
              Cvoid, (Ptr{InternalRefListIterator},), x.internalptr)
    end

    return T(refl),iter

end


function Base.iterate(::RefList{T}, iter) where T

    refl = ccall((:next_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalReflection},Ptr{InternalRefListIterator}),
                 iter.lastrefl, iter.internalptr)

    if refl == C_NULL
        iter.internalptr = C_NULL   # libcrystfel already freed it
        return nothing
    end

    iter.lastrefl = refl

    return T(refl),iter

end


Base.IteratorEltype(::RefList{T}) where T = T
Base.isdone(iter::RefListIterator) = ((iter.internalptr == C_NULL) && (iter.lastrefl != C_NULL))
Base.length(reflist::RefList) = ccall((:num_reflections, libcrystfel),
                                      Cint, (Ptr{InternalRefList},), reflist.internalptr)


function Base.getindex(reflist::RefList{T}, indices) where T

    refl = @ccall libcrystfel.find_refl(reflist.internalptr::Ptr{InternalRefList},
                                        indices[1]::Cint,
                                        indices[2]::Cint,
                                        indices[3]::Cint)::Ptr{InternalReflection}

    if refl == C_NULL
        return nothing
    else
        return T(refl)
    end
end


function loadreflist(filename::AbstractString)

    psym = Ref{Cstring}()
    out = ccall((:read_reflections_2, libcrystfel),
                Ptr{InternalRefList}, (Cstring,Ref{Cstring}),
                filename, psym)
    if out == C_NULL
        throw(ArgumentError("Failed to load reflection list"))
    end

    symmetryname = unsafe_string(psym[])
    return RefList{MergedReflection}(out, SymOpList(symmetryname))

end


function savereflist!(reflist::RefList{MergedReflection}, filename::AbstractString)

    r = @ccall libcrystfel.write_reflist_2(filename::Cstring,
                                           reflist.internalptr::Ptr{InternalRefList},
                                           reflist.symmetry.internalptr::Ptr{InternalSymOpList})::Cint

    if r != 0
        throw(ErrorException("Failed to save reflection list"))
    end

end


function Base.show(io::IO, ::MIME"text/plain", reflist::RefList{MergedReflection})
    println(io, "Merged reflection list in point group ", symmetry_name(reflist.symmetry))
    print(io, "   h    k    l  intensity  σ(intens) nmeas")
    let n = 0
        for refl in Iterators.take(reflist, 11)
            if n == 10
                # We have printed 10 already, and are here again.  Truncate...
                print(io, "\n   ⋮    ⋮    ⋮          ⋮          ⋮     ⋮")
                break
            end
            let ind = refl.indices
                write(io, "\n")
                @printf(io, "%4i %4i %4i %10.2f %10.2f %5i", ind[1], ind[2], ind[3],
                        refl.intensity, refl.sigintensity, refl.nmeasurements)
                n += 1
            end
        end
    end
end


function Base.show(io::IO, ::MIME"text/plain", reflist::RefList{UnmergedReflection})
    println(io, "Unmerged reflection list in point group ", symmetry_name(reflist.symmetry))
    print(io, "   h    k    l  intensity  σ(intens)      fs      ss  panel")
    let n = 0
        for refl in Iterators.take(reflist, 11)
            if n == 10
                # We have printed 10 already, and are here again.  Truncate...
                print(io, "\n   ⋮    ⋮    ⋮          ⋮          ⋮       ⋮       ⋮      ⋮")
                break
            end
            let ind = refl.indices,
                pos = refl.detectorposition
                write(io, "\n")
                @printf(io, "%4i %4i %4i %10.2f %10.2f %7.2f %7.2f %6i",
                        ind[1], ind[2], ind[3], refl.intensity, refl.sigintensity,
                        pos.fs, pos.ss, pos.panelnumber)
                n += 1
            end
        end
    end
end


function detectorpos(refl::Reflection)
    pfs = Ref{Cdouble}(0)
    pss = Ref{Cdouble}(0)
    pn = ccall((:get_panel_number, libcrystfel),
               Cint, (Ptr{InternalReflection},), refl.internalptr)
    ccall((:get_detector_pos, libcrystfel),
          Cint, (Ptr{InternalReflection},Ref{Cdouble},Ref{Cdouble}),
          refl.internalptr, pfs, pss)
    (fs=pfs[], ss=pss[], panelnumber=pn)
end


function indices(refl::Reflection)
    h = Ref{Cint}(0)
    k = Ref{Cint}(0)
    l = Ref{Cint}(0)
    ccall((:get_indices, libcrystfel),
          Cint, (Ptr{InternalReflection},Ref{Cint},Ref{Cint},Ref{Cint}),
          refl.internalptr, h, k, l)
    (h[], k[], l[])
end


function symmetricindices(refl::Reflection)
    h = Ref{Cint}(0)
    k = Ref{Cint}(0)
    l = Ref{Cint}(0)
    ccall((:get_symmetric_indices, libcrystfel),
          Cint, (Ptr{InternalReflection},Ref{Cint},Ref{Cint},Ref{Cint}),
          refl.internalptr, h, k, l)
    (h[], k[], l[])
end


function Base.getproperty(refl::Reflection, name::Symbol)
    if name === :intensity
        ccall((:get_intensity, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :sigintensity
        ccall((:get_esd_intensity, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :partiality
        ccall((:get_partiality, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :khalf
        ccall((:get_khalf, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :kpred
        ccall((:get_kpred, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :excitationerror
        ccall((:get_exerr, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :lorentzfactor
        ccall((:get_lorentz, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :phase
        ccall((:get_phase, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :peak
        ccall((:get_peak, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :meanbackground
        ccall((:get_mean_bg, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :temp1
        ccall((:get_temp1, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :temp2
        ccall((:get_temp2, libcrystfel),
              Cdouble, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :nmeasurements
        ccall((:get_redundancy, libcrystfel),
              Cint, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :flag
        ccall((:get_flag, libcrystfel),
              Cint, (Ptr{InternalReflection},), refl.internalptr)
    elseif name === :detectorposition
        detectorpos(refl)
    elseif name === :indices
        indices(refl)
    elseif name === :symmetricindices
        symmetricindices(refl)
    else
        getfield(refl, name)
    end
end


function Base.setproperty!(refl::Reflection, name::Symbol, val)
    if name === :intensity
        @ccall libcrystfel.set_intensity(refl.internalptr::Ptr{InternalReflection},
                                         val::Cdouble)::Cvoid
    elseif name === :sigintensity
        @ccall libcrystfel.set_esd_intensity(refl.internalptr::Ptr{InternalReflection},
                                             val::Cdouble)::Cvoid
    elseif name === :partiality
        @ccall libcrystfel.set_partiality(refl.internalptr::Ptr{InternalReflection},
                                          val::Cdouble)::Cvoid
    elseif name === :khalf
        @ccall libcrystfel.set_khalf(refl.internalptr::Ptr{InternalReflection},
                                     val::Cdouble)::Cvoid
    elseif name === :kpred
        @ccall libcrystfel.set_kpred(refl.internalptr::Ptr{InternalReflection},
                                     val::Cdouble)::Cvoid
    elseif name === :excitationerror
        @ccall libcrystfel.set_exerr(refl.internalptr::Ptr{InternalReflection},
                                     val::Cdouble)::Cvoid
    elseif name === :lorentzfactor
        @ccall libcrystfel.set_lorentz(refl.internalptr::Ptr{InternalReflection},
                                       val::Cdouble)::Cvoid
    elseif name === :phase
        @ccall libcrystfel.set_phase(refl.internalptr::Ptr{InternalReflection},
                                     val::Cdouble)::Cvoid
    elseif name === :peak
        @ccall libcrystfel.set_peak(refl.internalptr::Ptr{InternalReflection},
                                    val::Cdouble)::Cvoid
    elseif name === :meanbackground
        @ccall libcrystfel.set_mean_bg(refl.internalptr::Ptr{InternalReflection},
                                       val::Cdouble)::Cvoid
    elseif name === :temp1
        @ccall libcrystfel.set_temp1(refl.internalptr::Ptr{InternalReflection},
                                     val::Cdouble)::Cvoid
    elseif name === :temp2
        @ccall libcrystfel.set_temp2(refl.internalptr::Ptr{InternalReflection},
                                     val::Cdouble)::Cvoid
    elseif name === :nmeasurements
        @ccall libcrystfel.set_redundancy(refl.internalptr::Ptr{InternalReflection},
                                          val::Cint)::Cvoid
    elseif name === :flag
        @ccall libcrystfel.set_flag(refl.internalptr::Ptr{InternalReflection},
                                    val::Cint)::Cvoid
    elseif name === :detectorposition
        @ccall libcrystfel.set_detector_pos(refl.internalptr::Ptr{InternalReflection},
                                            val.fs::Cdouble, val.ss::Cdouble)::Cvoid
        @ccall libcrystfel.set_panel_number(refl.internalptr::Ptr{InternalReflection},
                                            val.panelnumber::Cint)::Cvoid
    elseif name === :indices
        throw(ErrorException("Cannot set the indices of a Reflection"))
    elseif name === :symmetricindices
        @ccall libcrystfel.set_symmetric_indices(refl.internalptr::Ptr{InternalReflection},
                                                 val[1]::Cint, val[2]::Cint, val[3]::Cint)::Cvoid
    else
        setfield!(refl, name, val)
    end
end


function Base.propertynames(::UnmergedReflection; private=false)
    names = (:intensity,:sigintensity,:partiality,:khalf,:kpred,:lorentzfactor,
             :excitationerror,:phase,:peak,:meanbackground,:temp1,:temp2,
             :nmeasurements,:flag,:detectorposition,:indices,:symmetricindices)
    if private
        tuple(names..., :internalptr)
    else
        names
    end
end


function Base.propertynames(::MergedReflection; private=false)
    names = (:intensity,:sigintensity,:phase,:temp1,:temp2,
             :nmeasurements,:flag,:indices,:symmetricindices)
    if private
        tuple(names..., :internalptr)
    else
        names
    end
end


function Base.show(io::IO, refl::MergedReflection)
    write(io, "MergedReflection(")
    show(io, refl.indices)
    write(io, ", intensity=")
    show(io, refl.intensity)
    write(io, ", σ(intensity)=")
    show(io, refl.sigintensity)
    write(io, ", nmeasurements=")
    show(io, refl.nmeasurements)
    write(io, ")")
end


function Base.show(io::IO, refl::UnmergedReflection)
    write(io, "UnmergedReflection(")
    show(io, refl.indices)
    write(io, " at ");
    show(io, refl.detectorposition)
    write(io, ", intensity=")
    show(io, refl.intensity)
    write(io, ", σ(intensity)=")
    show(io, refl.sigintensity)
    write(io, ")")
end


end  # of module
