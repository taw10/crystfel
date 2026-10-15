using CrystFEL
using Random

function addpeaks(image, cr, peaklist, ch)
    reflist = predictreflections(cr, image)
    for refl in reflist
        if randn() > ch
            let dpos = refl.detectorposition
                push!(peaklist, round(dpos.fs), round(dpos.ss), dpos.panelnumber, 100.0)
            end
        end
    end
end


function writedata(fh, name, peaklist)

    println(fh, "int "*name*"_n_peaks = ", length(peaklist), ";")

    println(fh, "float "*name*"_peak_fs[] = {")
    for p in peaklist
        println(fh, p.fs, ",")
    end
    println(fh, "};")

    println(fh, "float "*name*"_peak_ss[] = {")
    for p in peaklist
        println(fh, p.ss, ",")
    end
    println(fh, "};")

    println(fh, "int "*name*"_peak_pn[] = {")
    for p in peaklist
        println(fh, p.panelnumber, ",")
    end
    println(fh, "};")

end


dtempl = loaddatatemplate("tests/indextest.geom")
image = Image(dtempl)

let fh = open("tests/indextest.h", "w")

    # One crystal
    peaklist = PeakList()
    cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
    addpeaks(image, Crystal(rotatecell(cell)), peaklist, 0.0)
    writedata(fh, "t1", peaklist)

    # Two crystals
    peaklist = PeakList()
    cell = UnitCell(MonoclinicLattice, PrimitiveCell, 123, 45, 80, 90, 97, 90)
    addpeaks(image, Crystal(rotatecell(cell)), peaklist, 0.0)
    addpeaks(image, Crystal(rotatecell(cell)), peaklist, 0.0)
    writedata(fh, "t2", peaklist)

    # One crystal, small cell
    peaklist = PeakList()
    cell = UnitCell(TriclinicLattice, PrimitiveCell, 12.3, 4.5, 8.0, 120, 97, 90)
    addpeaks(image, Crystal(rotatecell(cell), profileradius=5e7), peaklist, -1000)
    writedata(fh, "t3", peaklist)

    # Two crystals, small cell
    peaklist = PeakList()
    cell = UnitCell(TriclinicLattice, PrimitiveCell, 12.3, 4.5, 8.0, 120, 97, 90)
    addpeaks(image, Crystal(rotatecell(cell), profileradius=5e7), peaklist, -1000)
    addpeaks(image, Crystal(rotatecell(cell), profileradius=5e7), peaklist, -1000)
    writedata(fh, "t4", peaklist)
end
