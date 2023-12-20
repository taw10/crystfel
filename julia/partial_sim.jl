using CrystFEL

dtempl = loaddatatemplate("julia/alignment-test.geom")
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
let st = Stream("partials.stream", "w", dtempl)
    for i in 1:10
        println("Writing image ", i)
        image = Image(dtempl)
        cr = Crystal(rotatecell(cell))
        reflist = predictreflections(cr, image)
        chunkwrite(st, image)
    end
end
