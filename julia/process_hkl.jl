using CrystFEL


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
            @addmeasurement(correction(refl.intensity, cr), 1.0,
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


let st = Stream("input.stream", "r")
    merged = mergereflections((i,cr)->i, allcrystals(st), SymOpList("2/m_uab"))
    savereflist!(merged, "merged.hkl")
end
