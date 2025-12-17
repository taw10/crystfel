module DetGeoms
export DetGeom, DetGeomPanel, DetGeomGroup, InternalDetGeom, copydetgeom
export findpanel, findgroup

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


mutable struct InternalDetGeomPanelGroup
    name::Cstring
    n_children::Cint
    parent::Ptr{InternalDetGeomPanelGroup}
    serial::Cint
    cx::Cdouble
    cy::Cdouble
    cz::Cdouble
    children::Ptr{Ptr{InternalDetGeomPanelGroup}}
    panel::Ptr{InternalPanel}
end


mutable struct InternalDetGeom
    panels::Ptr{InternalPanel}
    n_panels::Cint
    top_group::Ptr{InternalDetGeomPanelGroup}
end


mutable struct DetGeomPanel
    name
    corner
    fsvec
    ssvec
    pixel_pitch
    adu_per_photon
    max_adu
    w
    h
end

struct DetGeomGroup
    serial
    name
    children
    center    # For rotation
end


mutable struct DetGeom
    panels::Vector{DetGeomPanel}
    topgroup::DetGeomGroup
end


function copy_dg_panel(p::InternalPanel)
    DetGeomPanel(unsafe_string(p.name),
                 [p.cx, p.cy, p.cz],
                 [p.fsx, p.fsy, p.fsz],
                 [p.ssx, p.ssy, p.ssz],
                 p.pixel_pitch, p.adu_per_photon, p.max_adu,
                 p.w, p.h)
end


function copydggroup(panels, grp)
    grpdata = unsafe_load(grp, 1)
    name = unsafe_string(grpdata.name)
    children = DetGeomGroup[]
    if grpdata.n_children > 0
        for i in 1:grpdata.n_children
            cdata = unsafe_load(grpdata.children, i)
            push!(children, copydggroup(panels, cdata))
        end
        DetGeomGroup(grpdata.serial, name, children,
                     [grpdata.cx, grpdata.cy, grpdata.cz])
    else
        name = unsafe_string(grpdata.name)
        pn = findfirst(x->x.name==name, panels)
        DetGeomGroup(grpdata.serial, name, [panels[pn]],
                     [grpdata.cx, grpdata.cy, grpdata.cz])
    end
end


function copydetgeom(dg::Ptr{InternalDetGeom})
    dgdata = unsafe_load(dg, 1)
    panels = DetGeomPanel[]
    for i in 1:dgdata.n_panels
        paneldata = unsafe_load(dgdata.panels, i)
        push!(panels, copy_dg_panel(paneldata))
    end
    topgroup = copydggroup(panels, dgdata.top_group)
    DetGeom(panels, topgroup)
end

findpanel(dg::DetGeom, name) = dg.panels[findfirst(x->x.name==name, dg.panels)]

function seq100(ser, head)
    if ser == 0
        ()
    elseif ser < 100
        (head...,ser % 100)
    else
        seq100(ser ÷ 100, (head..., ser % 100))
    end
end

findgroup(dg::DetGeom, ser) = findgroup(dg.topgroup, ser, seq100(ser÷100, ()))

function findgroup(group::DetGeomGroup, ser, path)
    if length(path) == 0
        @assert(group.serial == ser)
        return group
    end
    findgroup(group.children[path[1]], ser, path[2:end])
end


function Base.show(io::IO, p::DetGeomPanel)
    write(io, "Panel(")
    write(io, "name=\"")
    write(io, p.name)
    write(io, "\", corner=")
    show(io, p.corner)
    write(io, "), fs=(")
    show(io, p.fsvec)
    write(io, "), ss=(")
    show(io, p.ssvec)
    write(io, "), size=(")
    show(io, p.w); write(io, ", "); show(io, p.h)
    write(io, "))")
end

function showgroup(io, p::DetGeomPanel, prefix)
    println(io, prefix, "Panel ", p.name, " (",p.w,"×",p.h,")")
end

function showgroup(io, grp::DetGeomGroup, prefix)
    println(io, prefix, "Group ", grp.name, " (serial number ", grp.serial, ")")
    for child in grp.children
        showgroup(io, child, prefix*" ")
    end
end

function Base.show(io::IO, ::MIME"text/plain", dg::DetGeom)
    println(io, "Detector geometry structure with ", length(dg.panels), " panels")
    showgroup(io, dg.topgroup, " ")
end

end  # of module
