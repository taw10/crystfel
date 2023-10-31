using CrystFEL

dtempl = loaddatatemplate("julia/alignment-test.geom")
image = Image(dtempl)

cell = UnitCell(123, 45, 80, 90, 97, 90,
                centering='P', latticetype=MonoclinicLattice, uniqueaxis='b')

truth = predictreflections(image, cell)
