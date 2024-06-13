using CrystFEL
using GLMakie

f = Figure(size=(1024,1024))
ax = Axis3(f[1,1], aspect=(1,1,1),
           xlabel="x (meters)", ylabel="y (meters)", zlabel="z (meters)")
sl = Slider(f[2,1], range=0:0.01:1)

dtempl = loaddatatemplate("cxidb-21.geom")
cfimage = CrystFEL.Image(dtempl)

motions = let fh = open("millepede.res", "r")
    readline(fh)   # discard header line
    motions = map(eachline(fh)) do line
        sgroupser,sdelta,_ = split(line)
        (parse(Int,sgroupser), parse(Float64,sdelta))
    end
    motions
end

function translate!(panel::DetGeomPanel, vec)
    panel.corner += vec/panel.pixel_pitch
end

function translate!(group::DetGeomGroup, vec)
    for child in group.children
        translate!(child, vec)
    end
end

function rotate2d!(p, c, ang)
    nx = c[1] + (p[1]-c[1])*cos(ang) - (p[2]-c[2])*sin(ang);
    ny = c[2] + (p[1]-c[1])*sin(ang) + (p[2]-c[2])*cos(ang);
    p[1] = nx
    p[2] = ny
end

function rotatex!(panel::DetGeomPanel, center, angle)
    rotate2d!(view(panel.corner, 2:3), center[2:3]/panel.pixel_pitch, angle);
    rotate2d!(view(panel.fsvec, 2:3), [0,0], angle);
    rotate2d!(view(panel.ssvec, 2:3), [0,0], angle);
end

function rotatey!(panel::DetGeomPanel, center, angle)
    rotate2d!(view(panel.corner, 3:-2:1), center[3:-2:1]/panel.pixel_pitch, angle);
    rotate2d!(view(panel.fsvec, 3:-2:1), [0,0], angle);
    rotate2d!(view(panel.ssvec, 3:-2:1), [0,0], angle);
end

function rotatez!(panel::DetGeomPanel, center, angle)
    rotate2d!(view(panel.corner, 1:2), center[1:2]/panel.pixel_pitch, angle);
    rotate2d!(view(panel.fsvec, 1:2), [0,0], angle);
    rotate2d!(view(panel.ssvec, 1:2), [0,0], angle);
end

function rotatex!(group::DetGeomGroup, center, angle)
    for child in group.children
        rotatex!(child, center, angle)
    end
end

function rotatey!(group::DetGeomGroup, center, angle)
    for child in group.children
        rotatey!(child, center, angle)
    end
end

function rotatez!(group::DetGeomGroup, center, angle)
    for child in group.children
        rotatez!(child, center, angle)
    end
end

triangles = lift(sl.value) do trscale

    dgmoved = cfimage.detgeom   # Makes a fresh copy of underlying C structure

    # Apply all `motions` in the same way as align_detector
    for m in motions
        param = m[1] % 100
        group = findgroup(dgmoved, m[1]-param)
        if param == 1
            translate!(group, trscale*[-m[2],0,0])
        elseif param == 2
            translate!(group, trscale*[0,-m[2],0])
        elseif param == 3
            translate!(group, trscale*[0,0,-m[2]])
        elseif param == 4
            rotatex!(group, group.center, trscale*-m[2])
        elseif param == 5
            rotatey!(group, group.center, trscale*-m[2])
        elseif param == 6
            rotatez!(group, group.center, trscale*-m[2])
        end
    end

    # Generate triangle mesh from panel descriptions
    points = Point3f[]
    for panel in dgmoved.panels

        corner = panel.pixel_pitch*panel.corner
        fs = panel.pixel_pitch*panel.w*panel.fsvec
        ss = panel.pixel_pitch*panel.h*panel.ssvec

        push!(points, corner)
        push!(points, corner+fs)
        push!(points, corner+fs+ss)
        push!(points, corner+fs+ss)
        push!(points, corner+ss)
        push!(points, corner)

    end

    return points

end

mesh!(triangles, color="#7f7f7f")
f
