module DetGeoms
export Panel, DetGeom

guardian = []

function protect(guardian, obj)
    push!(guardian, obj)
    obj
end

function unprotect(guardian, obj)
    let pos = findfirst(==(obj), guardian)
        if pos !== nothing
            deleteat!(guardian, pos)
        end
    end
end


mutable struct Panel
    name::Cstring
    cx::Cdouble
    cy::Cdouble
    cz::Cdouble
    pixel_pitch::Cdouble
    adu_per_photon::Cdouble
    max_adu::Cdouble
    fsx::Cdouble
    fsy::Cdouble
    fsz::Cdouble
    ssx::Cdouble
    ssy::Cdouble
    ssz::Cdouble
    w::Cint
    h::Cint
    group::Ptr{Cvoid}
end


"""
    Panel(name, width, height, (cnx, cny), clen, (fsx,fsy,fsz), (ssx,ssy,ssz), pixelsize, aduperphoton)

Create a panel for a CrystFEL `DetGeom`.

* `cnx` and `cny`: Corner position in pixel units
* `clen`: Corner z-position in meters
* `(fsx,fsy,fsz)`: Fast scan vector in pixel units
* `(ssx,ssy,ssz)`: Slow scan vector in pixel units
* `pixelsize`: Conversion factor from pixels to meters
* `aduperphoton`: Detector units per quantum, for error estimation

Additional keyword arguments:

* `max_adu`=Inf: Saturation value
* `group`: Panel group (for hierarchy)
"""
function Panel(name, width, height, corner::Tuple{Real, Real}, clen,
               fs::Tuple{Real,Real,Real}, ss::Tuple{Real,Real,Real},
               pixel_pitch, adu_per_photon,
               max_adu=Inf, group=C_NULL)

    myname = protect(guardian, deepcopy(name))

    p = Panel(pointer(myname),
              corner[1], corner[2], clen/pixel_pitch,
              pixel_pitch, adu_per_photon, max_adu,
              fs[1], fs[2], fs[3],
              ss[1], ss[2], ss[3],
              width, height, group)

    finalizer(p) do x
        unprotect(guardian, myname)
    end

end


function Base.show(io::IO, p::Panel)
    write(io, "Panel(")
    write(io, "name=\"")
    write(io, unsafe_string(p.name))
    write(io, "\", center=(")
    show(io, p.cx); write(io, ", "); show(io, p.cy); write(io, ", "); show(io, p.cz)
    write(io, "), fs=(")
    show(io, p.fsx); write(io, ", "); show(io, p.fsy); write(io, ", "); show(io, p.fsz)
    write(io, "), ss=(")
    show(io, p.ssx); write(io, ", "); show(io, p.ssy); write(io, ", "); show(io, p.ssz)
    write(io, "), size=(")
    show(io, p.w); write(io, ", "); show(io, p.h)
    write(io, "))")
end


"""
    DetGeom(panels; topgroup=g)

Create a CrystFEL `DetGeom` from a vector of `Panel`s.  Optionally set the
panel group which should be the top of the hierarchy.
"""
function DetGeom(panels; topgroup=nothing)
end

end  # of module
