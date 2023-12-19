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


function simulate_and_index(cell, image_true, dtempl_moved, mille, n)

    indexer = Indexer("asdf", dtempl_moved, cell, retry=false, multilattice=false, refine=true)

    for _ in 1:n

        cr = Crystal(rotatecell(cell))
        peaklist = sketch_pattern(image_true, cr)
        image_moved = Image(dtempl_moved)

        image_moved.peaklist = peaklist
        index(image_moved, indexer, mille=mille)

    end
end


dtempl_true = loaddatatemplate("julia/alignment-test.geom")
image_true = Image(dtempl_true)
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
dtempl_moved = loaddatatemplate("julia/alignment-test-moved.geom")
mille = Mille("mille.dat")
simulate_and_index(cell, image_true, dtempl_moved, mille, 100)
