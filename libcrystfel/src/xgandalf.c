/*
 * xgandalf.c
 *
 * Interface to XGANDALF indexer
 *
 * Copyright © 2017-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017-2018 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "xgandalf.h"

#ifdef HAVE_XGANDALF
#include <stdlib.h>

#include "utils.h"
#include "cell-utils.h"
#include "peaks.h"

#include "xgandalf/adaptions/crystfel/Lattice.h"
#include "xgandalf/adaptions/crystfel/ExperimentSettings.h"
#include "xgandalf/adaptions/crystfel/IndexerPlain.h"

struct xgandalf_private_data {
	IndexerPlain *indexer;
	reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A;

	IndexingMethod indm;
	UnitCell *cellTemplate;
	Lattice_t sampleRealLatticeReduced_A;   //same as cellTemplate
};

#define FAKE_DETECTOR_DISTANCE (0.1)
#define FAKE_DETECTOR_RADIUS (0.1)
#define FAKE_BEAM_ENERGY (1)
#define FAKE_DIVERGENCE_ANGLE_DEG (0.05)
#define FAKE_NON_MONOCHROMATICITY (0.005)
#define FAKE_REFLECTION_RADIUS (0.0001)

#define MAX_ASSEMBLED_LATTICES_COUNT (10)

static void reduceCell(UnitCell* cell);
static void makeRightHanded(UnitCell* cell);

int run_xgandalf(struct image *image, void *ipriv)
{
	struct xgandalf_private_data *xgandalf_private_data =
	                                      (struct xgandalf_private_data*) ipriv;
	reciprocalPeaks_1_per_A_t *reciprocalPeaks_1_per_A =
	                          &(xgandalf_private_data->reciprocalPeaks_1_per_A);

	int peakCountMax = image_feature_count(image->features);
	reciprocalPeaks_1_per_A->peakCount = 0;
	for (int i = 0; i < peakCountMax && i < MAX_PEAK_COUNT_FOR_INDEXER; i++) {
		struct imagefeature *f;
		f = image_get_feature(image->features, i);
		if (f == NULL) {
			continue;
		}

		reciprocalPeaks_1_per_A->coordinates_x[reciprocalPeaks_1_per_A->peakCount]
		                                                        = f->rx * 1e-10;
		reciprocalPeaks_1_per_A->coordinates_y[reciprocalPeaks_1_per_A->peakCount]
		                                                        = f->ry * 1e-10;
		reciprocalPeaks_1_per_A->coordinates_z[reciprocalPeaks_1_per_A->peakCount]
		                                                        = f->rz * 1e-10;
		reciprocalPeaks_1_per_A->peakCount++;
	}

	Lattice_t assembledLattices[MAX_ASSEMBLED_LATTICES_COUNT];
	int assembledLatticesCount;
	IndexerPlain_index(xgandalf_private_data->indexer,
	                   assembledLattices,
	                   &assembledLatticesCount,
	                   MAX_ASSEMBLED_LATTICES_COUNT,
	                   *reciprocalPeaks_1_per_A,
	                   NULL);

	if (assembledLatticesCount > 0) { //no multi-lattice at the moment
		assembledLatticesCount = 1;
	}

	int goodLatticesCount = assembledLatticesCount;
	for (int i = 0; i < assembledLatticesCount && i < 1; i++) {
		reorderLattice(&(xgandalf_private_data->sampleRealLatticeReduced_A),
		               &assembledLattices[i]);

		UnitCell *uc;
		uc = cell_new();

		Lattice_t *l = &assembledLattices[i];

		cell_set_cartesian(uc, l->ax * 1e-10, l->ay * 1e-10, l->az * 1e-10,
		                       l->bx * 1e-10, l->by * 1e-10, l->bz * 1e-10,
		                       l->cx * 1e-10, l->cy * 1e-10, l->cz * 1e-10);
		makeRightHanded(uc);

		if (validate_cell(uc)) {
			STATUS("Problem with returned cell!\n");
		}

		Crystal *cr = crystal_new();
		if (cr == NULL) {
			ERROR("Failed to allocate crystal.\n");
			return 0;
		}
		crystal_set_cell(cr, uc);
		image_add_crystal(image, cr);

	}

	return goodLatticesCount;
}

void *xgandalf_prepare(IndexingMethod *indm,
	                   UnitCell *cell,
	                   struct xgandalf_options *xgandalf_opts)
{
	struct xgandalf_private_data *xgandalf_private_data =
	                              malloc(sizeof(struct xgandalf_private_data));
	allocReciprocalPeaks(&(xgandalf_private_data->reciprocalPeaks_1_per_A));
	xgandalf_private_data->indm = *indm;
	xgandalf_private_data->cellTemplate = NULL;

	float tolerance = xgandalf_opts->tolerance;
	samplingPitch_t samplingPitch = xgandalf_opts->sampling_pitch;
	gradientDescentIterationsCount_t gradientDescentIterationsCount =
	                                xgandalf_opts->grad_desc_iterations;

	if (*indm & INDEXING_USE_CELL_PARAMETERS) {
		xgandalf_private_data->cellTemplate = cell;

		UnitCell* primitiveCell = uncenter_cell(cell, NULL);
		reduceCell(primitiveCell);

		double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
		int ret = cell_get_reciprocal(primitiveCell, &asx, &asy, &asz,
		                                             &bsx, &bsy, &bsz,
		                                             &csx, &csy, &csz);
		if (ret != 0) {
			ERROR("cell_get_reciprocal did not finish properly!");
		}

		Lattice_t sampleReciprocalLattice_1_per_A = {
		        .ax = asx * 1e-10, .ay = asy * 1e-10, .az = asz * 1e-10,
		        .bx = bsx * 1e-10, .by = bsy * 1e-10, .bz = bsz * 1e-10,
		        .cx = csx * 1e-10, .cy = csy * 1e-10, .cz = csz * 1e-10 };

		double ax, ay, az, bx, by, bz, cx, cy, cz;
		ret = cell_get_cartesian(primitiveCell, &ax, &ay, &az,
		                                        &bx, &by, &bz,
		                                        &cx, &cy, &cz);
		if (ret != 0) {
			ERROR("cell_get_cartesian did not finish properly!");
		}
		Lattice_t sampleRealLatticeReduced_A = {
		        .ax = ax * 1e10, .ay = ay * 1e10, .az = az * 1e10,
		        .bx = bx * 1e10, .by = by * 1e10, .bz = bz * 1e10,
		        .cx = cx * 1e10, .cy = cy * 1e10, .cz = cz * 1e10 };
		xgandalf_private_data->sampleRealLatticeReduced_A =
		                                             sampleRealLatticeReduced_A;

		ExperimentSettings *experimentSettings =
				ExperimentSettings_new(FAKE_BEAM_ENERGY,
		                               FAKE_DETECTOR_DISTANCE,
		                               FAKE_DETECTOR_RADIUS,
		                               FAKE_DIVERGENCE_ANGLE_DEG,
		                               FAKE_NON_MONOCHROMATICITY,
		                               sampleReciprocalLattice_1_per_A,
		                               tolerance,
		                               FAKE_REFLECTION_RADIUS);

		xgandalf_private_data->indexer = IndexerPlain_new(experimentSettings);
		IndexerPlain_setSamplingPitch(xgandalf_private_data->indexer,
		                              samplingPitch);
		IndexerPlain_setGradientDescentIterationsCount(
		                                        xgandalf_private_data->indexer,
		                                        gradientDescentIterationsCount);

		if (xgandalf_opts->no_deviation_from_provided_cell) {
			IndexerPlain_setRefineWithExactLattice(
			                                     xgandalf_private_data->indexer,
			                                     1);
		}

		ExperimentSettings_delete(experimentSettings);
		cell_free(primitiveCell);
	}
	else {
		Lattice_t sampleRealLatticeReduced_A = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		xgandalf_private_data->sampleRealLatticeReduced_A =
		                                             sampleRealLatticeReduced_A;

		ExperimentSettings *experimentSettings =
		   ExperimentSettings_new_nolatt(FAKE_BEAM_ENERGY,
		                                 FAKE_DETECTOR_DISTANCE,
		                                 FAKE_DETECTOR_RADIUS,
		                                 FAKE_DIVERGENCE_ANGLE_DEG,
		                                 FAKE_NON_MONOCHROMATICITY,
		                                 xgandalf_opts->minLatticeVectorLength_A,
		                                 xgandalf_opts->maxLatticeVectorLength_A,
		                                 FAKE_REFLECTION_RADIUS);

		xgandalf_private_data->indexer = IndexerPlain_new(experimentSettings);
		IndexerPlain_setSamplingPitch(xgandalf_private_data->indexer,
		                              samplingPitch);
		IndexerPlain_setGradientDescentIterationsCount(
		                                       xgandalf_private_data->indexer,
		                                       gradientDescentIterationsCount);

		ExperimentSettings_delete(experimentSettings);
	}

	/* Flags that XGANDALF knows about */
	*indm &= INDEXING_METHOD_MASK
	        | INDEXING_USE_CELL_PARAMETERS;

	return xgandalf_private_data;
}

void xgandalf_cleanup(void *pp)
{
	struct xgandalf_private_data *xgandalf_private_data =
			                                 (struct xgandalf_private_data*) pp;

	freeReciprocalPeaks(xgandalf_private_data->reciprocalPeaks_1_per_A);
	IndexerPlain_delete(xgandalf_private_data->indexer);
}

static void reduceCell(UnitCell *cell)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };

	reduceLattice(&l);

	cell_set_cartesian(cell, l.ax, l.ay, l.az,
	                         l.bx, l.by, l.bz,
	                         l.cx, l.cy, l.cz);

	makeRightHanded(cell);
}

static void makeRightHanded(UnitCell *cell)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	if (!right_handed(cell)) {
		cell_set_cartesian(cell, -ax, -ay, -az, -bx, -by, -bz, -cx, -cy, -cz);
	}
}

const char *xgandalf_probe(UnitCell *cell)
{
    return "xgandalf";
}

#else

int run_xgandalf(struct image *image, void *ipriv)
{
	ERROR("This copy of CrystFEL was compiled without XGANDALF support.\n");
	return 0;
}

void *xgandalf_prepare(IndexingMethod *indm, UnitCell *cell,
		struct xgandalf_options *xgandalf_opts)
{
	ERROR("This copy of CrystFEL was compiled without XGANDALF support.\n");
	ERROR("To use XGANDALF indexing, recompile with XGANDALF.\n");
	return NULL;
}

void xgandalf_cleanup(void *pp)
{
}

const char *xgandalf_probe(UnitCell *cell)
{
    return NULL;
}

#endif // HAVE_XGANDALF
