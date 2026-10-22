/*
 * xgandalf_shim.h
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
 */


#ifndef XGANDALF_SHIM_H
#define XGANDALF_SHIM_H

#define MAX_PEAK_COUNT_FOR_INDEXER 20000
#define MAX_PEAK_COUNT_FOR_PROJECTION 80000

#ifdef __cplusplus
    extern "C" {
        namespace xgandalf {
#endif

typedef struct
{
    float *coordinates_x;
    float *coordinates_y;
    float *coordinates_z;
    int peakCount;
} xgreciprocalPeaks_1_per_A_t;

typedef struct
{
    float ax;
    float ay;
    float az;

    float bx;
    float by;
    float bz;

    float cx;
    float cy;
    float cz;
} xgLattice_t;

typedef struct
{
    float matrixElement_0_0;
    float matrixElement_1_0;
    float matrixElement_2_0;
    float matrixElement_0_1;
    float matrixElement_1_1;
    float matrixElement_2_1;
    float matrixElement_0_2;
    float matrixElement_1_2;
    float matrixElement_2_2;
} xgLatticeTransform_t;

void xgallocReciprocalPeaks(xgreciprocalPeaks_1_per_A_t *reciprocalPeaks_1_per_A);
void xgfreeReciprocalPeaks(xgreciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A);

typedef struct ExperimentSettings xgExperimentSettings;

xgExperimentSettings *xgExperimentSettings_new_nolatt(float beamEenergy_eV,
                                                      float detectorDistance_m,
                                                      float detectorRadius_m,
                                                      float divergenceAngle_deg,
                                                      float nonMonochromaticity,
                                                      float minRealLatticeVectorLength_A,
                                                      float maxRealLatticeVectorLength_A,
                                                      float reflectionRadius_1_per_A);

xgExperimentSettings *xgExperimentSettings_new(float beamEenergy_eV,
                                               float detectorDistance_m,
                                               float detectorRadius_m,
                                               float divergenceAngle_deg,
                                               float nonMonochromaticity,
                                               const xgLattice_t sampleReciprocalLattice_1A,
                                               float tolerance,
                                               float reflectionRadius_1_per_A);

void xgExperimentSettings_delete(xgExperimentSettings *experimentSettings);

typedef enum {
    SAMPLING_PITCH_extremelyLoose = 0,
    SAMPLING_PITCH_loose = 1,
    SAMPLING_PITCH_standard = 2,
    SAMPLING_PITCH_dense = 3,
    SAMPLING_PITCH_extremelyDense = 4,

    SAMPLING_PITCH_standardWithSeondaryMillerIndices = 5,
    SAMPLING_PITCH_denseWithSeondaryMillerIndices = 6,
    SAMPLING_PITCH_extremelyDenseWithSeondaryMillerIndices = 7,

    SAMPLING_PITCH_lastEnum
} xgsamplingPitch_t;

typedef enum {
    GRADIENT_DESCENT_ITERATION_COUNT_verryFew = 0,
    GRADIENT_DESCENT_ITERATION_COUNT_few = 1,
    GRADIENT_DESCENT_ITERATION_COUNT_standard = 2,
    GRADIENT_DESCENT_ITERATION_COUNT_many = 3,
    GRADIENT_DESCENT_ITERATION_COUNT_manyMany = 4,
    GRADIENT_DESCENT_ITERATION_COUNT_extremelyMany = 5,

    GRADIENT_DESCENT_ITERATION_COUNT_custom = 6,

    GRADIENT_DESCENT_ITERATION_COUNT_lastEnum
} xggradientDescentIterationsCount_t;

typedef struct IndexerPlain xgIndexerPlain;

xgIndexerPlain *xgIndexerPlain_new(xgExperimentSettings* experimentSettings);
void xgIndexerPlain_delete(xgIndexerPlain* indexerPlain);

void xgIndexerPlain_setSamplingPitch(xgIndexerPlain *indexerPlain,
                                   xgsamplingPitch_t samplingPitch);
void xgIndexerPlain_setGradientDescentIterationsCount(xgIndexerPlain *indexerPlain,
                                                    xggradientDescentIterationsCount_t gradientDescentIterationsCount);
void xgIndexerPlain_setRefineWithExactLattice(xgIndexerPlain *indexerPlain, int flag);
void xgIndexerPlain_setMaxPeaksToUseForIndexing(xgIndexerPlain *indexerPlain,
                                                int maxPeaksToUseForIndexing);

void xgIndexerPlain_index(xgIndexerPlain *indexerPlain,
                          xgLattice_t *assembledLattices,
                          int *assembledLatticesCount,
                          int maxAssambledLatticesCount,
                          xgreciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A,
                          int *peakCountOnLattices);

void xgreorderLattice(const xgLattice_t* prototype, xgLattice_t* lattice);
void xgreduceLattice(xgLattice_t* lattice, xgLatticeTransform_t* appliedReductionTransform);
void xgrestoreLattice(xgLattice_t* lattice, xgLatticeTransform_t* appliedReductionTransform);

#ifdef __cplusplus
}
} /* extern "C" */
#endif

#endif /* XGANDALF_SHIM_H */
