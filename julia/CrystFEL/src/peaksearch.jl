module PeakSearch

import ..CrystFEL: libcrystfel
import ..CrystFEL.Images: InternalImage
import ..CrystFEL.PeakLists: PeakList, InternalPeakList

export zaefpeaks, peakfinder8, peakfinder9


function tf10(val)
    if val
        return 1
    else
        return 0
    end
end


function zaefpeaks(image; threshold=100, mingrad=100000, minsnr=5,
        radiusinn=4, radiusmid=5, radiusout=7, usesaturated=true)
    out = @ccall libcrystfel.search_peaks(image.internalptr::Ptr{InternalImage},
                                          threshold::Cfloat,
                                          mingrad::Cfloat,
                                          minsnr::Cfloat,
                                          radiusinn::Cdouble,
                                          radiusmid::Cdouble,
                                          radiusout::Cdouble,
                                          tf10(usesaturated)::Cint)::Ptr{InternalPeakList}
    if out == C_NULL
        throw(ErrorException("Peak search failed"))
    end
    PeakList(out)
end



function peakfinder8(image; threshold=100, minsnr=5, minpix=2, maxpix=200,
        localbg=3, minres=0, maxres=5000, usesaturated=true, maxpeaks=2000)
    out = @ccall libcrystfel.peakfinder8(image.internalptr::Ptr{InternalImage},
                                         maxpeaks::Cint,
                                         threshold::Cfloat,
                                         minsnr::Cfloat,
                                         minpix::Cint,
                                         maxpix::Cint,
                                         localbg::Cint,
                                         minres::Cint,
                                         maxres::Cint,
                                         tf10(usesaturated)::Cint,
                                         0::Cint,
                                         C_NULL::Ptr{Cvoid})::Ptr{InternalPeakList}
    if out == C_NULL
        throw(ErrorException("Peak search failed"))
    end
    PeakList(out)
end


function peakfinder9(image; minsnrbig=7, minsnrpeak=6, minsnrwhole=5, minbgsig=11,
        brightpxcut=-Inf, window=5)
    out = @ccall libcrystfel.search_peaks_peakfinder9(image.internalptr::Ptr{InternalImage},
                                                      minsnrbig::Cfloat,
                                                      minsnrpeak::Cfloat,
                                                      minsnrwhole::Cfloat,
                                                      minbgsig::Cfloat,
                                                      brightpxcut::Cfloat,
                                                      window::Cint)::Ptr{InternalPeakList}
    if out == C_NULL
        throw(ErrorException("Peak search failed"))
    end
    PeakList(out)
end


end   # of module
