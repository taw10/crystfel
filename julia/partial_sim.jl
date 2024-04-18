using CrystFEL

dtempl = loaddatatemplate("input.geom")
cell = UnitCell(OrthorhombicLattice, PrimitiveCell, 40.0, 50.0, 60.0)
full = loadreflist("input.hkl")
let st = Stream("partials.stream", "w", dtempl)
    for i in 1:10
        println("Writing image ", i)
        image = Image(dtempl)
        image.serial = i
        image.filename = "simulation_" * string(i)
        image.ev = "//"
        cr = Crystal(rotatecell(cell))
        reflections = predictreflections(cr, image)
        calculatepartialities!(reflections, cr, image, model=XSphereModel)
        for refl in reflections
            f = full[asymmetricindices(full.symmetry, refl.indices...)...]
            if f !== nothing
                refl.intensity = f.intensity * refl.partiality * refl.lorentzfactor
            else
                # No matching "full" intensity - can't do anything
                # Reflections with zero measurements won't be written to file
                refl.nmeasurements = 0
                println("Not found: ", refl.indices)
            end
        end
        push!(image, cr, reflections)
        chunkwrite(st, image)
    end
end
