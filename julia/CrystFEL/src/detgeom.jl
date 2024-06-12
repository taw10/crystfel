module DetGeoms
export DetGeom, DetGeomPanel, InternalDetGeom, copydetgeom

mutable struct InternalPanel
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


mutable struct InternalDetGeom
    panels::Ptr{InternalPanel}
    n_panels::Cint
    top_group::Ptr{Cvoid}
end


mutable struct DetGeomPanel
    name
    cx
    cy
    cz
    pixel_pitch
    adu_per_photon
    max_adu
    fsx
    fsy
    fsz
    ssx
    ssy
    ssz
    w
    h
end


mutable struct DetGeom
    panels::Vector{DetGeomPanel}
end




function copy_dg_panel(p::InternalPanel)
    DetGeomPanel(unsafe_string(p.name),
                 p.cx, p.cy, p.cz,
                 p.pixel_pitch, p.adu_per_photon, p.max_adu,
                 p.fsx, p.fsy, p.fsz,
                 p.ssx, p.ssy, p.ssz,
                 p.w, p.h)
end


function copydetgeom(dg::Ptr{InternalDetGeom})
    dgdata = unsafe_load(dg, 1)
    detgeom = DetGeom(DetGeomPanel[])
    for i in 1:dgdata.n_panels
        paneldata = unsafe_load(dgdata.panels, i)
        push!(detgeom.panels, copy_dg_panel(paneldata))
    end
    detgeom
end


function Base.show(io::IO, p::DetGeomPanel)
    write(io, "Panel(")
    write(io, "name=\"")
    write(io, p.name)
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
