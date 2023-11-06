module RefLists

using Printf
import ..CrystFEL: libcrystfel
import ..CrystFEL.Symmetry: SymOpList, InternalSymOpList, symmetry_name
export RefList, loadreflist, reflections


# The internal libcrystfel structures, not exposed directly
# We only ever have e.g. a Ptr{InternalRefList}, never a real InternalRefList
mutable struct InternalRefList end
mutable struct InternalReflection end
mutable struct InternalRefListIterator end

# The Julian exposed types
mutable struct RefList
    internalptr::Ptr{InternalRefList}
    symmetry::SymOpList
end

mutable struct Reflection
    internalptr::Ptr{InternalReflection}
end

mutable struct RefListIterator
    reflist::RefList
    lastrefl::Ptr{InternalReflection}
    internalptr::Ptr{InternalRefListIterator}
end


function reflections(reflist::RefList)
    iter = RefListIterator(reflist, C_NULL, C_NULL)
    finalizer(iter) do x
        if x.internalptr != C_NULL
            ccall((:free_reflistiterator, libcrystfel),
                  Cvoid, (Ptr{InternalRefListIterator},), x.internalptr)
        end
    end
    return iter
end


function Base.iterate(iter::RefListIterator)

    rli = Ref{Ptr{InternalRefListIterator}}(C_NULL)
    refl = ccall((:first_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalRefList},Ref{Ptr{InternalRefListIterator}}),
                 iter.reflist.internalptr, rli)

    if refl == C_NULL
        throw(ArgumentError("Failed to find first reflection in list"))
    end

    iter.lastrefl = refl
    iter.internalptr = rli[]

    return Reflection(refl),iter

end


function Base.iterate(iter::RefListIterator, _)

    refl = ccall((:next_refl, libcrystfel),
                 Ptr{InternalReflection}, (Ptr{InternalReflection},Ptr{InternalRefListIterator}),
                 iter.lastrefl, iter.internalptr)

    if refl == C_NULL
        iter.internalptr = C_NULL   # libcrystfel already freed it
        return nothing
    end

    iter.lastrefl = refl

    return Reflection(refl),iter

end


Base.IteratorEltype(::RefListIterator) = Reflection

Base.isdone(iter::RefListIterator) = ((iter.internalptr == C_NULL) && (iter.lastrefl != C_NULL))


function loadreflist(filename::AbstractString)

    psym = Ref{Cstring}()
    out = ccall((:read_reflections_2, libcrystfel),
                Ptr{InternalRefList}, (Cstring,Ref{Cstring}),
                filename, psym)
    if out == C_NULL
        throw(ArgumentError("Failed to load reflection list"))
    end

    symmetryname = unsafe_string(psym[])
    return RefList(out, SymOpList(symmetryname))

end


function Base.show(io::IO, ::MIME"text/plain", reflist::RefList)
    println(io, "Reflection list in point group ", symmetry_name(reflist.symmetry))
    print(io, "   h    k    l  intensity")
    let n = 0
        for refl in Iterators.take(reflections(reflist), 11)
            if n == 10
                # We have printed 10 already, and are here again.  Truncate...
                print(io, "\n   ⋮    ⋮    ⋮          ⋮")
                break
            end
            let ind = refl.indices
                write(io, "\n")
                @printf(io, "%4i %4i %4i %10.2f", ind[1], ind[2], ind[3], refl.intensity)
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


function Base.propertynames(::Reflection; private=false)
    names = (:intensity,:sigintensity,:partiality,:khalf,:kpred,:lorentzfactor,
             :excitationerror,:phase,:peak,:meanbackground,:temp1,:temp2,
             :nmeasurements,:flag,:detectorposition,:indices,:symmetricindices)
    if private
        tuple(names..., :internalptr)
    else
        names
    end
end


function Base.show(io::IO, refl::Reflection)
    write(io, "Reflection(")
    show(io, refl.indices)
    write(io, ", intensity=")
    show(io, refl.intensity)
    write(io, ")")
end

end  # of module
