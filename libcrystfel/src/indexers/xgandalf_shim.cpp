/*
 * xgandalf_shim.cpp
 *
 * Copyright © 2019-2026 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
 *   2026 Thomas White <taw@physics.org>
 *
 * Adapted from XGANDALF.
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <Eigen/Dense>
#include <xgandalf/IndexerPlain.h>
#include <xgandalf/Lattice.h>

#include "xgandalf_shim.h"

namespace xgandalf {

    ExperimentSettings *ExperimentSettings_new_nolatt(float beamEenergy_eV,
                                                      float detectorDistance_m,
                                                      float detectorRadius_m,
                                                      float divergenceAngle_deg,
                                                      float nonMonochromaticity,
                                                      float minRealLatticeVectorLength_A,
                                                      float maxRealLatticeVectorLength_A,
                                                      float reflectionRadius_1_per_A)
    {
        return new ExperimentSettings(beamEenergy_eV,
                                      detectorDistance_m,
                                      detectorRadius_m,
                                      divergenceAngle_deg,
                                      nonMonochromaticity,
                                      minRealLatticeVectorLength_A,
                                      maxRealLatticeVectorLength_A,
                                      reflectionRadius_1_per_A);
    }

    ExperimentSettings *ExperimentSettings_new(float beamEenergy_eV,
                                               float detectorDistance_m,
                                               float detectorRadius_m,
                                               float divergenceAngle_deg,
                                               float nonMonochromaticity,
                                               const Lattice_t sampleReciprocalLattice_1A,
                                               float tolerance,
                                               float reflectionRadius_1_per_A)
    {
        const Lattice_t& l = sampleReciprocalLattice_1A;

        Eigen::Matrix3f lattice;
        lattice << l.ax, l.bx, l.cx, l.ay, l.by, l.cy, l.az, l.bz, l.cz;

        return new ExperimentSettings(beamEenergy_eV,
                                      detectorDistance_m,
                                      detectorRadius_m,
                                      divergenceAngle_deg,
                                      nonMonochromaticity,
                                      lattice,
                                      tolerance,
                                      reflectionRadius_1_per_A);
    }

    void ExperimentSettings_delete(ExperimentSettings *experimentSettings)
    {
        delete experimentSettings;
    }

    void allocReciprocalPeaks(reciprocalPeaks_1_per_A_t *reciprocalPeaks_1_per_A)
    {
        reciprocalPeaks_1_per_A->coordinates_x = new float[MAX_PEAK_COUNT_FOR_INDEXER];
        reciprocalPeaks_1_per_A->coordinates_y = new float[MAX_PEAK_COUNT_FOR_INDEXER];
        reciprocalPeaks_1_per_A->coordinates_z = new float[MAX_PEAK_COUNT_FOR_INDEXER];
    }

    void freeReciprocalPeaks(reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A)
    {
        delete[] reciprocalPeaks_1_per_A.coordinates_x;
        delete[] reciprocalPeaks_1_per_A.coordinates_y;
        delete[] reciprocalPeaks_1_per_A.coordinates_z;
    }

    IndexerPlain *IndexerPlain_new(ExperimentSettings *experimentSettings)
    {
        return new IndexerPlain(*experimentSettings);
    }

    void IndexerPlain_delete(IndexerPlain *indexerPlain)
    {
        delete indexerPlain;
    }

    void IndexerPlain_setSamplingPitch(IndexerPlain *indexerPlain, samplingPitch_t samplingPitch)
    {
        IndexerPlain::SamplingPitch pitch;
        switch (samplingPitch)
        {
            case SAMPLING_PITCH_extremelyLoose:
                pitch = IndexerPlain::SamplingPitch::extremelyLoose;
                break;
            case SAMPLING_PITCH_loose:
                pitch = IndexerPlain::SamplingPitch::loose;
                break;
            case SAMPLING_PITCH_standard:
                pitch = IndexerPlain::SamplingPitch::standard;
                break;
            case SAMPLING_PITCH_dense:
                pitch = IndexerPlain::SamplingPitch::dense;
                break;
            case SAMPLING_PITCH_extremelyDense:
                pitch = IndexerPlain::SamplingPitch::extremelyDense;
                break;

            case SAMPLING_PITCH_standardWithSeondaryMillerIndices:
                pitch = IndexerPlain::SamplingPitch::standardWithSeondaryMillerIndices;
                break;
            case SAMPLING_PITCH_denseWithSeondaryMillerIndices:
                pitch = IndexerPlain::SamplingPitch::denseWithSeondaryMillerIndices;
                break;
            case SAMPLING_PITCH_extremelyDenseWithSeondaryMillerIndices:
                pitch = IndexerPlain::SamplingPitch::extremelyDenseWithSeondaryMillerIndices;
                break;
            default:
                pitch = IndexerPlain::SamplingPitch::standard;
                break;
        }

        indexerPlain->setSamplingPitch(pitch);
    }

    void IndexerPlain_setGradientDescentIterationsCount(IndexerPlain *indexerPlain,
            gradientDescentIterationsCount_t gradientDescentIterationsCount)
    {
        IndexerPlain::GradientDescentIterationsCount iterationsCount;

        switch (gradientDescentIterationsCount)
        {
            case GRADIENT_DESCENT_ITERATION_COUNT_verryFew:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::exremelyFew;
                break;
            case GRADIENT_DESCENT_ITERATION_COUNT_few:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::few;
                break;
            case GRADIENT_DESCENT_ITERATION_COUNT_standard:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::standard;
                break;
            case GRADIENT_DESCENT_ITERATION_COUNT_many:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::many;
                break;
            case GRADIENT_DESCENT_ITERATION_COUNT_manyMany:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::manyMany;
                break;
            case GRADIENT_DESCENT_ITERATION_COUNT_extremelyMany:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::extremelyMany;
                break;

            default:
                iterationsCount = IndexerPlain::GradientDescentIterationsCount::standard;
                break;
        }

        indexerPlain->setGradientDescentIterationsCount(iterationsCount);
    }

    void IndexerPlain_setRefineWithExactLattice(IndexerPlain *indexerPlain, int flag)
    {
        indexerPlain->setRefineWithExactLattice((bool)flag);
    }

    void IndexerPlain_setMaxPeaksToUseForIndexing(IndexerPlain *indexerPlain, int maxPeaksToUseForIndexing)
    {
        indexerPlain->setMaxPeaksToUseForIndexing(maxPeaksToUseForIndexing);
    }

    void IndexerPlain_index(IndexerPlain *indexerPlain,
                            Lattice_t *assembledLattices,
                            int *assembledLatticesCount,
                            int maxAssambledLatticesCount,
                            reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A,
                            int *peakCountOnLattices)
    {
        Eigen::Matrix3Xf reciprocalPeaks_1_per_A_matrix(3, reciprocalPeaks_1_per_A.peakCount);
        for (int i = 0; i < reciprocalPeaks_1_per_A.peakCount; i++)
        {
            reciprocalPeaks_1_per_A_matrix.col(i) << reciprocalPeaks_1_per_A.coordinates_x[i],
                                                     reciprocalPeaks_1_per_A.coordinates_y[i],
                                                     reciprocalPeaks_1_per_A.coordinates_z[i];
        }

        std::vector<Lattice> assembledLatticesVector;
        std::vector<int> peakCountOnLatticesVector;
        indexerPlain->index(assembledLatticesVector, reciprocalPeaks_1_per_A_matrix, peakCountOnLatticesVector);

        for (*assembledLatticesCount = 0;
             (size_t)*assembledLatticesCount < assembledLatticesVector.size()
                    && *assembledLatticesCount < maxAssambledLatticesCount;
             (*assembledLatticesCount)++)
        {
            Eigen::Matrix3f basis = assembledLatticesVector[*assembledLatticesCount].getBasis();
            assembledLattices[*assembledLatticesCount].ax = basis(0, 0);
            assembledLattices[*assembledLatticesCount].ay = basis(1, 0);
            assembledLattices[*assembledLatticesCount].az = basis(2, 0);
            assembledLattices[*assembledLatticesCount].bx = basis(0, 1);
            assembledLattices[*assembledLatticesCount].by = basis(1, 1);
            assembledLattices[*assembledLatticesCount].bz = basis(2, 1);
            assembledLattices[*assembledLatticesCount].cx = basis(0, 2);
            assembledLattices[*assembledLatticesCount].cy = basis(1, 2);
            assembledLattices[*assembledLatticesCount].cz = basis(2, 2);

            if (peakCountOnLattices != NULL)
            {
                peakCountOnLattices[*assembledLatticesCount] = peakCountOnLatticesVector[*assembledLatticesCount];
            }
        }
    }

    void reorderLattice(const Lattice_t *prototype, Lattice_t *lattice)
    {
        Eigen::Matrix3f prototypeBasis;
        prototypeBasis << prototype->ax, prototype->bx, prototype->cx,
                          prototype->ay, prototype->by, prototype->cy,
                          prototype->az, prototype->bz, prototype->cz;
        Lattice prototypeLattice(prototypeBasis);

        Eigen::Matrix3f basis;
        basis << lattice->ax, lattice->bx, lattice->cx,
                 lattice->ay, lattice->by, lattice->cy,
                 lattice->az, lattice->bz, lattice->cz;
        Lattice latticeWrapper(basis);

        latticeWrapper.reorder(prototypeLattice.getBasisVectorNorms(),
                               prototypeLattice.getBasisVectorAngles_deg());
        basis = latticeWrapper.getBasis();

        lattice->ax = basis(0, 0);
        lattice->ay = basis(1, 0);
        lattice->az = basis(2, 0);
        lattice->bx = basis(0, 1);
        lattice->by = basis(1, 1);
        lattice->bz = basis(2, 1);
        lattice->cx = basis(0, 2);
        lattice->cy = basis(1, 2);
        lattice->cz = basis(2, 2);
    }

    void reduceLattice(Lattice_t *lattice, LatticeTransform_t *appliedReductionTransform)
    {
        Eigen::Matrix3f basis, reducedBasis;
        basis << lattice->ax, lattice->bx, lattice->cx,
                 lattice->ay, lattice->by, lattice->cy,
                 lattice->az, lattice->bz, lattice->cz;
        Lattice latticeWrapper(basis);

        reducedBasis = latticeWrapper.minimize().getBasis();

        lattice->ax = reducedBasis(0, 0);
        lattice->ay = reducedBasis(1, 0);
        lattice->az = reducedBasis(2, 0);
        lattice->bx = reducedBasis(0, 1);
        lattice->by = reducedBasis(1, 1);
        lattice->bz = reducedBasis(2, 1);
        lattice->cx = reducedBasis(0, 2);
        lattice->cy = reducedBasis(1, 2);
        lattice->cz = reducedBasis(2, 2);

        Eigen::Matrix3f transformMatrix = basis.inverse() * reducedBasis;

        appliedReductionTransform->matrixElement_0_0 = transformMatrix(0, 0);
        appliedReductionTransform->matrixElement_1_0 = transformMatrix(1, 0);
        appliedReductionTransform->matrixElement_2_0 = transformMatrix(2, 0);
        appliedReductionTransform->matrixElement_0_1 = transformMatrix(0, 1);
        appliedReductionTransform->matrixElement_1_1 = transformMatrix(1, 1);
        appliedReductionTransform->matrixElement_2_1 = transformMatrix(2, 1);
        appliedReductionTransform->matrixElement_0_2 = transformMatrix(0, 2);
        appliedReductionTransform->matrixElement_1_2 = transformMatrix(1, 2);
        appliedReductionTransform->matrixElement_2_2 = transformMatrix(2, 2);
    }

    void restoreLattice(Lattice_t *lattice, LatticeTransform_t *appliedReductionTransform)
    {
        Eigen::Matrix3f transformMatrix;

        transformMatrix(0, 0) = appliedReductionTransform->matrixElement_0_0;
        transformMatrix(1, 0) = appliedReductionTransform->matrixElement_1_0;
        transformMatrix(2, 0) = appliedReductionTransform->matrixElement_2_0;
        transformMatrix(0, 1) = appliedReductionTransform->matrixElement_0_1;
        transformMatrix(1, 1) = appliedReductionTransform->matrixElement_1_1;
        transformMatrix(2, 1) = appliedReductionTransform->matrixElement_2_1;
        transformMatrix(0, 2) = appliedReductionTransform->matrixElement_0_2;
        transformMatrix(1, 2) = appliedReductionTransform->matrixElement_1_2;
        transformMatrix(2, 2) = appliedReductionTransform->matrixElement_2_2;

        Eigen::Matrix3f basis, reducedBasis;
        reducedBasis << lattice->ax, lattice->bx, lattice->cx,
                        lattice->ay, lattice->by, lattice->cy,
                        lattice->az, lattice->bz, lattice->cz;

        basis = reducedBasis * transformMatrix.inverse();

        lattice->ax = basis(0, 0);
        lattice->ay = basis(1, 0);
        lattice->az = basis(2, 0);
        lattice->bx = basis(0, 1);
        lattice->by = basis(1, 1);
        lattice->bz = basis(2, 1);
        lattice->cx = basis(0, 2);
        lattice->cy = basis(1, 2);
        lattice->cz = basis(2, 2);
    }

} // namespace xgandalf
