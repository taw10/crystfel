using CrystFEL
using Random
using Plots
using MillepedeII

# "Simulate" a diffraction pattern
function sketch_pattern(image, cr)
    reflist = predictreflections(cr, image)
    peaklist = PeakList()
    for refl in reflist
        if randn() > 0
            let dpos = refl.detectorposition
                push!(peaklist, round(dpos.fs), round(dpos.ss), dpos.panelnumber, 100.0)
            end
        end
    end
    return peaklist
end


function simulate_and_index(cell, image_true, dtempl_moved, mille, n)

    indexer = Indexer("asdf", dtempl_moved, cell, retry=false, multilattice=false, refine=true)

    nidx = 0
    for i in 1:n

        # Create a diffraction pattern for a random orientation
        cr = Crystal(rotatecell(cell))
        peaklist = sketch_pattern(image_true, cr)

        # Make an image with the correct spot positions,
        # but with an incorrect geometry
        image_moved = Image(dtempl_moved)
        image_moved.peaklist = peaklist

        # Index the pattern (and store Mille data),
        # based on the incorrect geometry
        index(image_moved, indexer, mille=mille)

        if image_moved.n_crystals > 0
            nidx += 1
        end

        if i % 100 == 0
            print("*")
        else
            print(".")
        end

    end
    println("")

    println("Indexed ", nidx, " out of ", n, " frames")

end


dtempl_true = loaddatatemplate("scripts/alignment-test.geom")
image_true = Image(dtempl_true)
cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
dtempl_moved = loaddatatemplate("scripts/alignment-test.geom")
translategroup!(dtempl_moved, "q1", 200e-6, 0, 0)
let mille = Mille("mille.dat")
    simulate_and_index(cell, image_true, dtempl_moved, mille, 100)
    close(mille)
end


function ploth(a, n, offs, label)
    histogram(map(x->x.residual,
                  filter(x->in(n, keys(x.globalgradients)),
                         a[offs:3:end])),
              label=label)
end


function plotfs(filename)

    t = loadmille(filename)
    l = @layout([q0 q1; q2 q3])
    a = collect(Iterators.flatten(t))

    q0 = ploth(a, 101, 1, "q0")
    q1 = ploth(a, 201, 1, "q1")
    q2 = ploth(a, 301, 1, "q2")
    q3 = ploth(a, 401, 1, "q3")
    plot(q0, q1, q2, q3, layout=l, plot_title="Fast scan residual / px")

end


function plotss(filename)

    t = loadmille(filename)
    l = @layout([q0 q1; q2 q3])
    a = collect(Iterators.flatten(t))

    q0 = ploth(a, 102, 2, "q0")
    q1 = ploth(a, 202, 2, "q1")
    q2 = ploth(a, 302, 2, "q2")
    q3 = ploth(a, 402, 2, "q3")
    plot(q0, q1, q2, q3, layout=l, plot_title="Slow scan residual / px")

end


function plotr(filename)
    t = loadmille(filename)
    l = @layout([q0 q1; q2 q3])
    a = collect(Iterators.flatten(t))
    meas = histogram(map(x->x.residual, a[3:3:end]))
    plot(meas, plot_title="Excitation error residual / UNITS")
end
