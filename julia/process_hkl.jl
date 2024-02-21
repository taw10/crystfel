using CrystFEL


function mergereflections(correction, crystalrefls, sym)

    merged = RefList{MergedReflection}(sym)

    for (cr,reflections) in crystalrefls

        for refl in reflections

            indices = asymmetricindices(sym, refl.indices)
            model_version = get!(merged, indices)

            w = 1.0
            mean = model_version.intensity
            sumweight = model_version.temp1
            M2 = model_version.temp2
            temp = w + sumweight
            delta = correction(refl.intensity, cr) - mean
            R = delta * w / temp
            model_version.intensity = mean + R
            model_version.temp1 = temp
            model_version.temp2 = M2 + sumweight * delta * R
            model_version.nmeasurements += 1

        end
    end

    for refl in merged
        if refl.nmeasurements > 1
            refl.sigintensity = sqrt(refl.temp2/refl.temp1)/sqrt(refl.nmeasurements)
        else
            refl.nmeasurements = 0
        end
    end

    return merged

end


let st = Stream("input.stream", "r")
    merged = mergereflections((i,cr)->i, allcrystals(st), SymOpList("2/m_uab"))
    savereflist!(merged, "merged.hkl")
end
