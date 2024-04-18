module DiffractionCalculations

import ..CrystFEL: libcrystfel
import ..CrystFEL.Images: InternalImage, Image
import ..CrystFEL.Crystals: InternalCrystal, Crystal
import ..CrystFEL.RefLists: RefList, UnmergedReflection, InternalRefList
import ..CrystFEL.Symmetry: SymOpList
export predictreflections, calculatepartialities!
export PartialityModel, UnityModel, XSphereModel, OffsetModel, RandomModel, GeneralGaussianModel


"""
Enumeration of the available partiality models.
"""
@enum PartialityModel begin
    UnityModel
    XSphereModel
    OffsetModel
    RandomModel
    GeneralGaussianModel
end


function predictreflections(cr::Crystal, image::Image; maxres=1e10)

    refls = @ccall libcrystfel.predict_to_res(cr.internalptr::Ptr{InternalCrystal},
                                              image.internalptr::Ptr{InternalImage},
                                              maxres::Cdouble)::Ptr{InternalRefList}
    sym = SymOpList("1")
    return RefList{UnmergedReflection}(refls, sym)
end


function calculatepartialities!(reflist::RefList{UnmergedReflection},
        cr::Crystal, image::Image; model=XSphereModel, maxres=1e10)

    @ccall libcrystfel.calculate_partialities(reflist.internalptr::Ptr{InternalRefList},
                                              cr.internalptr::Ptr{InternalCrystal},
                                              image.internalptr::Ptr{InternalImage},
                                              model::Cint)::Cvoid
end


end   # of module
