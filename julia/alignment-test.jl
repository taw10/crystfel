using CrystFEL
using Random
using Plots

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

    for i in 1:n

        cr = Crystal(rotatecell(cell))
        peaklist = sketch_pattern(image_true, cr)
        image_moved = Image(dtempl_moved)

        image_moved.peaklist = peaklist
        index(image_moved, indexer, mille=mille)

        if i % 100 == 0
            print("*")
        else
            print(".")
        end

    end
    println("")

end


dtempl_true = loaddatatemplate("julia/alignment-test.geom")
image_true = Image(dtempl_true)
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
dtempl_moved = loaddatatemplate("julia/alignment-test-moved.geom")


function plotresiduals(filename)

    l = @layout([q0 q1; q2 q3])
    a = collect(Iterators.flatten(t))

    function ploth(n, offs, label)
        histogram(map(x->x.residual,filter(x->in(n, keys(x.globalgradients)), a[offs:3:end])), label=label)
    end

    q0 = ploth(101, 1, "q0")
    q1 = ploth(201, 1, "q1")
    q2 = ploth(301, 1, "q2")
    q3 = ploth(401, 1, "q3")
    plot(q0, q1, q2, q3, layout=l, plot_title="Fast scan residual")

    #q0 = ploth(102, 2, "q0")
    #q1 = ploth(202, 2, "q1")
    #q2 = ploth(302, 2, "q2")
    #q3 = ploth(402, 2, "q3")
    #plot(q0, q1, q2, q3, layout=l, plot_title="Slow scan residual", reuse=false)

end
