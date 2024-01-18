using CrystFEL
using Random

# "Simulate" a diffraction pattern from the reflections
function sketch_pattern(image, cr)
    reflist = predictreflections(cr, image)
    peaklist = PeakList()
    for refl in reflist
        if randn() > 0
            let dpos = refl.detectorposition
                push!(peaklist, dpos.fs, dpos.ss, dpos.panelnumber, 100.0)
            end
        end
    end
    return peaklist
end


dtempl = loaddatatemplate("julia/alignment-test.geom")
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
indexer = Indexer("asdf", dtempl, cell, retry=false, multilattice=false, refine=true)

image = Image(dtempl)
cr = Crystal(rotatecell(cell))
image.peaklist = sketch_pattern(image, cr)
index(image, indexer)
