using CrystFEL

st = Stream("input.stream", "r")
sym = SymOpList("2/m_uab")
merged = RefList{MergedReflection}(sym)

for (cr,reflections) in allcrystals(st)

    for refl in reflections

        indices = asymmetricindices(sym, refl.indices)
        model_version = merged[indices]
        if model_version === nothing
            model_version = push!(merged, indices)
        end

        w = 1.0
        mean = model_version.intensity
        sumweight = model_version.temp1
        M2 = model_version.temp2
        temp = w + sumweight
        delta = refl.intensity - mean
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

savereflist!(merged, "merged.hkl")
