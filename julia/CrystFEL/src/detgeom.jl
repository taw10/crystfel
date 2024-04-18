module DetGeoms
export Panel

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


mutable struct DetGeom
    panels::Ptr{Panel}
    n_panels::Cint
    top_group::Ptr{Cvoid}
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


end  # of module
