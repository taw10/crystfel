using CrystFEL
using LinearAlgebra

function anglebetween(v1, v2)
    let v1n = norm(v1), v2n = norm(v2)
        return 2*atan(norm(v1*v2n - v2*v1n),
                      norm(v1*v2n + v2*v1n))
    end
end

let st = Stream("/home/twhite/experiments/cxidb-193/short.stream", "r"),
    merged = mergereflections(allcrystals(st), SymOpList("mmm")) do refl,crystal

        polfrac = 1.0
        polangle = 0.0

        lp = transpose(crystal.cell.reciprocalcartesian) * refl.symmetricindices
        tt = anglebetween([0,0,1], lp+[0,0,refl.kpred])
        phi = atan(lp[2], lp[1]) - polangle
		pol =         polfrac*(1.0 - cos(phi)*cos(phi)*sin(tt)*sin(tt)) +
		        (1.0-polfrac)*(1.0 - sin(phi)*sin(phi)*sin(tt)*sin(tt))

        return refl.intensity / pol

    end
    savereflist!(merged, "merged.hkl")
end
