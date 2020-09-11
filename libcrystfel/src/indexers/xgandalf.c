/*
 * xgandalf.c
 *
 * Interface to XGANDALF indexer
 *
 * Copyright © 2017-2020 Deutsches Elektronen-Synchrotron DESY,
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xgandalf.h"

#include <stdlib.h>

#include "utils.h"
#include "cell-utils.h"
#include "peaks.h"

#ifdef HAVE_XGANDALF
#include "xgandalf/adaptions/crystfel/Lattice.h"
#include "xgandalf/adaptions/crystfel/ExperimentSettings.h"
#include "xgandalf/adaptions/crystfel/IndexerPlain.h"
#endif

/** \file xgandalf.h */

struct xgandalf_options {
	unsigned int sampling_pitch;
	unsigned int grad_desc_iterations;
	float tolerance;
	unsigned int no_deviation_from_provided_cell;
	float minLatticeVectorLength_A;
	float maxLatticeVectorLength_A;
	int maxPeaksForIndexing;
};

#ifdef HAVE_XGANDALF

struct xgandalf_private_data {
	IndexerPlain *indexer;
	reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A;

	IndexingMethod indm;
	UnitCell *cellTemplate;
	Lattice_t sampleRealLattice_A;   //same as cellTemplate
	IntegerMatrix *centeringTransformation;
	LatticeTransform_t latticeReductionTransform;
};

#define FAKE_DETECTOR_DISTANCE (0.1)
#define FAKE_DETECTOR_RADIUS (0.1)
#define FAKE_BEAM_ENERGY (1)
#define FAKE_DIVERGENCE_ANGLE_DEG (0.05)
#define FAKE_NON_MONOCHROMATICITY (0.005)
#define FAKE_REFLECTION_RADIUS (0.0001)

#define MAX_ASSEMBLED_LATTICES_COUNT (10)

static void reduceCell(UnitCell* cell, LatticeTransform_t* appliedReductionTransform);
static void restoreCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform);
static void makeRightHanded(UnitCell* cell);

int run_xgandalf(struct image *image, void *ipriv)
{
	int i;
	struct xgandalf_private_data *xgandalf_private_data = (struct xgandalf_private_data*) ipriv;
	reciprocalPeaks_1_per_A_t *reciprocalPeaks_1_per_A = &(xgandalf_private_data->reciprocalPeaks_1_per_A);

	int peakCountMax = image_feature_count(image->features);
	reciprocalPeaks_1_per_A->peakCount = 0;
	for ( i = 0; i < peakCountMax && i < MAX_PEAK_COUNT_FOR_INDEXER; i++) {
		struct imagefeature *f;
		f = image_get_feature(image->features, i);
		if (f == NULL) {
			continue;
		}

		reciprocalPeaks_1_per_A->coordinates_x[reciprocalPeaks_1_per_A->peakCount] = f->rx * 1e-10;
		reciprocalPeaks_1_per_A->coordinates_y[reciprocalPeaks_1_per_A->peakCount] = f->ry * 1e-10;
		reciprocalPeaks_1_per_A->coordinates_z[reciprocalPeaks_1_per_A->peakCount] = f->rz * 1e-10;
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
	for ( i = 0; i < assembledLatticesCount && i < 1; i++) {
		reorderLattice(&(xgandalf_private_data->sampleRealLattice_A),
		                 &assembledLattices[i]);

		UnitCell *uc;
		uc = cell_new();

		Lattice_t *l = &assembledLattices[i];

		cell_set_cartesian(uc, l->ax * 1e-10, l->ay * 1e-10, l->az * 1e-10,
		                       l->bx * 1e-10, l->by * 1e-10, l->bz * 1e-10,
		                       l->cx * 1e-10, l->cy * 1e-10, l->cz * 1e-10);
		makeRightHanded(uc);

		if(xgandalf_private_data->cellTemplate != NULL){
			restoreCell(uc, &xgandalf_private_data->latticeReductionTransform);

			UnitCell *new_cell_trans = cell_transform_intmat(uc, xgandalf_private_data->centeringTransformation);
			cell_free(uc);
			uc = new_cell_trans;

			cell_set_lattice_type(new_cell_trans, cell_get_lattice_type(xgandalf_private_data->cellTemplate));
			cell_set_centering(new_cell_trans, cell_get_centering(xgandalf_private_data->cellTemplate));
			cell_set_unique_axis(new_cell_trans, cell_get_unique_axis(xgandalf_private_data->cellTemplate));
		}

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

void *xgandalf_prepare(IndexingMethod *indm, UnitCell *cell,
                       struct xgandalf_options *xgandalf_opts)
{
	struct xgandalf_private_data *xgandalf_private_data = malloc(sizeof(struct xgandalf_private_data));
	allocReciprocalPeaks(&(xgandalf_private_data->reciprocalPeaks_1_per_A));
	xgandalf_private_data->indm = *indm;
	xgandalf_private_data->cellTemplate = NULL;
	xgandalf_private_data->centeringTransformation = NULL;

	float tolerance = xgandalf_opts->tolerance;
	samplingPitch_t samplingPitch = xgandalf_opts->sampling_pitch;
	gradientDescentIterationsCount_t gradientDescentIterationsCount = xgandalf_opts->grad_desc_iterations;

	if (*indm & INDEXING_USE_CELL_PARAMETERS) {

		xgandalf_private_data->cellTemplate = cell;

		UnitCell* primitiveCell = uncenter_cell(cell, &xgandalf_private_data->centeringTransformation, NULL);

		reduceCell(primitiveCell, &xgandalf_private_data->latticeReductionTransform);

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
		Lattice_t sampleRealLattice_A = {
		        .ax = ax * 1e10, .ay = ay * 1e10, .az = az * 1e10,
		        .bx = bx * 1e10, .by = by * 1e10, .bz = bz * 1e10,
		        .cx = cx * 1e10, .cy = cy * 1e10, .cz = cz * 1e10 };
		xgandalf_private_data->sampleRealLattice_A = sampleRealLattice_A;

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

		if (xgandalf_opts->no_deviation_from_provided_cell) {
			IndexerPlain_setRefineWithExactLattice(xgandalf_private_data->indexer, 1);
		}

		ExperimentSettings_delete(experimentSettings);
		cell_free(primitiveCell);

	} else {

		Lattice_t sampleRealLattice_A = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		xgandalf_private_data->sampleRealLattice_A = sampleRealLattice_A;

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

		ExperimentSettings_delete(experimentSettings);
	}

	IndexerPlain_setSamplingPitch(xgandalf_private_data->indexer,
	        samplingPitch);
	IndexerPlain_setGradientDescentIterationsCount(xgandalf_private_data->indexer,
	        gradientDescentIterationsCount);
	IndexerPlain_setMaxPeaksToUseForIndexing(xgandalf_private_data->indexer,
			xgandalf_opts->maxPeaksForIndexing);

	/* Flags that XGANDALF knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_USE_CELL_PARAMETERS;

	return xgandalf_private_data;
}


void xgandalf_cleanup(void *pp)
{
	struct xgandalf_private_data *xgandalf_private_data = pp;

	freeReciprocalPeaks(xgandalf_private_data->reciprocalPeaks_1_per_A);
	IndexerPlain_delete(xgandalf_private_data->indexer);
	if(xgandalf_private_data->centeringTransformation != NULL){
		intmat_free(xgandalf_private_data->centeringTransformation);
	}
	free(xgandalf_private_data);
}

static void reduceCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };

	reduceLattice(&l, appliedReductionTransform);

	cell_set_cartesian(cell, l.ax, l.ay, l.az,
	                         l.bx, l.by, l.bz,
	                         l.cx, l.cy, l.cz);

	makeRightHanded(cell);
}

static void restoreCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform)
{

	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };

	restoreLattice(&l, appliedReductionTransform);

	cell_set_cartesian(cell, l.ax, l.ay, l.az,
	        l.bx, l.by, l.bz,
	        l.cx, l.cy, l.cz);

	makeRightHanded(cell);
}

static void makeRightHanded(UnitCell *cell)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	if ( !right_handed(cell) ) {
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

static void xgandalf_show_help()
{
	printf("Parameters for the TakeTwo indexing algorithm:\n"
"     --xgandalf-sampling-pitch\n"
"                           Sampling pitch: 0 (loosest) to 4 (most dense)\n"
"                            or with secondary Miller indices: 5 (loosest) to\n"
"                            7 (most dense).  Default: 6\n"
"     --xgandalf-grad-desc-iterations\n"
"                           Gradient descent iterations: 0 (few) to 5 (many)\n"
"                            Default: 4\n"
"     --xgandalf-fast-execution       Shortcut to set\n"
"                                     --xgandalf-sampling-pitch=2\n"
"                                     --xgandalf-grad-desc-iterations=3\n"
"     --xgandalf-tolerance  Relative tolerance of the lattice vectors.\n"
"                            Default is 0.02\n"
"     --xgandalf-no-deviation-from-provided-cell\n"
"                           Force the fitted cell to have the same lattice\n"
"                            parameters as the provided one\n"
"     --xgandalf-min-lattice-vector-length\n"
"                           Minimum possible lattice vector length in A.\n"
"                            Default: 30 A\n"
"     --xgandalf-max-lattice-vector-length\n"
"                           Maximum possible lattice vector length in A.\n"
"                            Default: 250 A\n"
"     --xgandalf-max-peaks\n"
"                           Maximum number of peaks used for indexing.\n"
"                           All peaks are used for refinement.\n"
"                            Default: 250\n"
);
}


int xgandalf_default_options(XGandalfOptions **opts_ptr)
{
	XGandalfOptions *opts;

	opts = malloc(sizeof(struct xgandalf_options));
	if ( opts == NULL ) return ENOMEM;

	opts->sampling_pitch = 6;
	opts->grad_desc_iterations = 4;
	opts->tolerance = 0.02;
	opts->no_deviation_from_provided_cell = 0;
	opts->minLatticeVectorLength_A = 30;
	opts->maxLatticeVectorLength_A = 250;
	opts->maxPeaksForIndexing = 250;

	*opts_ptr = opts;
	return 0;
}


static error_t xgandalf_parse_arg(int key, char *arg,
                                  struct argp_state *state)
{
	struct xgandalf_options **opts_ptr = state->input;
	int r;

	switch ( key ) {

		case ARGP_KEY_INIT :
		r = xgandalf_default_options(opts_ptr);
		if ( r ) return r;
		break;

		case 1 :
		xgandalf_show_help();
		return EINVAL;

		case 2 :
		if (sscanf(arg, "%u", &(*opts_ptr)->sampling_pitch) != 1) {
			ERROR("Invalid value for --xgandalf-sampling-pitch\n");
			return EINVAL;
		}
		break;

		case 3 :
		if (sscanf(arg, "%u", &(*opts_ptr)->grad_desc_iterations) != 1) {
			ERROR("Invalid value for --xgandalf-grad-desc-iterations\n");
			return EINVAL;
		}
		break;

		case 4 :
		if (sscanf(arg, "%f", &(*opts_ptr)->tolerance) != 1) {
			ERROR("Invalid value for --xgandalf-tolerance\n");
			return EINVAL;
		}
		break;

		case 5 :
		(*opts_ptr)->no_deviation_from_provided_cell = 1;
		break;

		case 6 :
		if (sscanf(arg, "%f", &(*opts_ptr)->minLatticeVectorLength_A) != 1) {
			ERROR("Invalid value for --xgandalf-min-lattice-vector-length\n");
			return EINVAL;
		}
		break;

		case 7 :
		if (sscanf(arg, "%f", &(*opts_ptr)->maxLatticeVectorLength_A) != 1) {
			ERROR("Invalid value for --xgandalf-max-lattice-vector-length\n");
			return EINVAL;
		}
		break;

		case 8 :
		(*opts_ptr)->sampling_pitch = 2;
		(*opts_ptr)->grad_desc_iterations = 3;
		break;

		case 9 :
		if (sscanf(arg, "%i", &(*opts_ptr)->maxPeaksForIndexing) != 1) {
			ERROR("Invalid value for --xgandalf-max-peaks\n");
			return EINVAL;
		}
		break;

	}

	return 0;
}


static struct argp_option xgandalf_options[] = {

	{"help-xgandalf", 1, NULL, OPTION_NO_USAGE,
	 "Show options for XGANDALF indexing algorithm", 99},

	{"xgandalf-sampling-pitch", 2, "pitch", OPTION_HIDDEN, NULL},
	{"xgandalf-sps", 2, "pitch", OPTION_HIDDEN, NULL},

	{"xgandalf-grad-desc-iterations", 3, "n", OPTION_HIDDEN, NULL},
	{"xgandalf-gdis", 3, "n", OPTION_HIDDEN, NULL},

	{"xgandalf-tolerance", 4, "t", OPTION_HIDDEN, NULL},
	{"xgandalf-tol", 4, "t", OPTION_HIDDEN, NULL},

	{"xgandalf-no-deviation-from-provided-cell", 5, NULL, OPTION_HIDDEN, NULL},
	{"xgandalf-ndfpc", 5, NULL, OPTION_HIDDEN, NULL},

	{"xgandalf-min-lattice-vector-length", 6, "len", OPTION_HIDDEN, NULL},
	{"xgandalf-min-lvl", 6, "len", OPTION_HIDDEN, NULL},

	{"xgandalf-max-lattice-vector-length", 7, "len", OPTION_HIDDEN, NULL},
	{"xgandalf-max-lvl", 7, "len", OPTION_HIDDEN, NULL},

	{"xgandalf-fast-execution", 8, NULL, OPTION_HIDDEN, NULL},

	{"xgandalf-max-peaks", 9, "n", OPTION_HIDDEN, NULL},

	{0}
};


struct argp xgandalf_argp = { xgandalf_options, xgandalf_parse_arg,
                              NULL, NULL, NULL, NULL, NULL };
