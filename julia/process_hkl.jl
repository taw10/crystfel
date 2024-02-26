using CrystFEL
using LinearAlgebra


macro addmeasurement(measurement, weight,
                     mean, sumweight, wksp)
    return quote
        delta = $(esc(measurement)) - $(esc(mean))
        newsumweight = $(esc(sumweight)) + $(esc(weight))
        R = delta * $(esc(weight)) / newsumweight
        $(esc(mean)) += R
        $(esc(wksp)) += $(esc(sumweight)) * delta * R
        $(esc(sumweight)) = newsumweight
    end
end

cstddev(nmeas, work1, work2) = sqrt(work2/work1)/sqrt(nmeas)


function mergereflections(correction, crystalrefls, sym)

    merged = RefList{MergedReflection}(sym)

    for (cr,reflections) in crystalrefls

        for refl in reflections

            indices = asymmetricindices(sym, refl.indices)
            model_version = get!(merged, indices)
            @addmeasurement(correction(refl, cr), 1.0,
                            model_version.intensity,
                            model_version.temp1,
                            model_version.temp2)
            model_version.nmeasurements += 1

        end
    end

    for refl in merged
        if refl.nmeasurements > 1
            refl.sigintensity = cstddev(refl.nmeasurements, refl.temp1, refl.temp2)
        else
            # Cannot delete a reflection from a list (especially not
            # while iterating), but setting nmeasurements to zero
            # prevents it from being written to the output file.
            refl.nmeasurements = 0
        end
    end

    return merged

end


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
