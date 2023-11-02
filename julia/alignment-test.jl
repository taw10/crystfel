using CrystFEL

dtempl = loaddatatemplate("julia/alignment-test.geom")
image = Image(dtempl)

cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)

truth = predictreflections(image, cell)
