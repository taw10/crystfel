using CrystFEL
using Random

# Create empty image for simulation purposes
dtempl = loaddatatemplate("julia/alignment-test.geom")
image = Image(dtempl)

# Create a crystal and calculate its reflections
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
cr = Crystal(cell)  # FIXME: Random rotation
truth = predictreflections(cr, image)

# "Simulate" a diffraction pattern from the reflections
function sketch_pattern(reflist)
    peaks = PeakList()
    for refl in reflist
        if randn() > 0
            let dpos = refl.detectorposition
                push!(peaks, dpos.fs, dpos.ss, dpos.panelnumber, 100.0)
            end
        end
    end
    return peaks
end

image.peaklist = sketch_pattern(truth)

# Index the pattern
indexer = Indexer("asdf", dtempl, cell, retry=false, multilattice=false, refine=true)
index(image, indexer)


# Utility routine for visualising peaks
function plotpanel(image, pn)
    x = []
    y = []
    for pk in image.peaklist
        if pk.panelnumber == pn
            push!(x, pk.fs)
            push!(y, pk.ss)
        end
    end
    scatter(x, y)
end
