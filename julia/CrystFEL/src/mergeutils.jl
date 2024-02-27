module MergeUtils

import ..CrystFEL: libcrystfel
using ..CrystFEL.RefLists
using ..CrystFEL.Symmetry
import ..CrystFEL.UnitCells: InternalUnitCell
export @addmeasurement, cstddev, mergereflections


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


struct Polarisation
    fraction::Cdouble
    angle::Cdouble
    disable::Cint
end

function polarisation_correction!(reflist, cell, polfrac, polangle)
    pol = Polarisation(polfrac, rad2deg(polangle), 0)
    @ccall libcrystfel.polarisation_correction(reflist.internalptr::Ptr{InternalRefList},
                                               cell.internalptr::Ptr{InternalUnitCell},
                                               pol::Ref{Polarisation})::Cvoid
end


function mergereflections(correction, crystalrefls, sym)

    merged = RefList{MergedReflection}(sym)

    for (cr,reflections) in crystalrefls

        polarisation_correction!(reflections, cr.cell, 1.0, 0.0)

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

mergereflections(crystalrefls, sym) = mergereflections((refl,cr)->refl.intensity, crystalrefls, sym)


end  # of module
