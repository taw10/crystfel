using CrystFEL

dtempl = loaddatatemplate("julia/alignment-test.geom")
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
let st = Stream("partials.stream", "w", dtempl)
    for i in 1:10
        println("Writing image ", i)
        image = Image(dtempl)
        image.serial = i
        image.filename = "simulation_" * string(i)
        image.ev = "//"
        cr = Crystal(rotatecell(cell))
        cr.reflections = predictreflections(cr, image)
        push!(image, cr)
        chunkwrite(st, image)
    end
end
