using CrystFEL
using Random

# Create empty image for simulation purposes
dtempl = loaddatatemplate("julia/alignment-test.geom")
image = Image(dtempl)

# Create a crystal and calculate its reflections
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
cr = Crystal(cell)  # FIXME: Random rotation
truth = predictreflections(cr, image)

# Sketch a diffraction pattern
image.peaklist = PeakList()
for refl in truth
    if randn() > 3
        let dpos = refl.detectorposition
            push!(image.peaklist, dpos.fs, dpos.ss, dpos.panelnumber, 100.0)
        end
    end
end

# Index the pattern
indexer = Indexer("asdf-cell", dtempl, cell, retry=false, multilattice=false, refine=true)
index(image, indexer)
