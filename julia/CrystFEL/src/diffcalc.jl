module DiffractionCalculations

import ..CrystFEL: libcrystfel
import ..CrystFEL.Images: InternalImage
import ..CrystFEL.Crystals: InternalCrystal
import ..CrystFEL.RefLists: RefList, UnmergedReflection, InternalRefList
import ..CrystFEL.Symmetry: SymOpList
export predictreflections


function predictreflections(cr, image; maxres=1e10)
    refls = ccall((:predict_to_res, libcrystfel),
                  Ptr{InternalRefList},
                  (Ptr{InternalCrystal}, Ptr{InternalImage}, Cdouble),
                  cr.internalptr, image.internalptr, maxres)
    sym = SymOpList("1")
    return RefList{UnmergedReflection}(refls, sym)
end


end   # of module
