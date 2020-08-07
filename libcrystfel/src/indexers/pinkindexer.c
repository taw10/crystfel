/*
 * pinkindexer.c
 *
 * Interface to PinkIndexer
 *
 * Copyright © 2017-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017-2019 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
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

#include "pinkindexer.h"


#include <stdlib.h>
#include <sys/errno.h>
#include <argp.h>

#include "utils.h"
#include "cell-utils.h"
#include "peaks.h"

struct pinkIndexer_options {
	unsigned int considered_peaks_count;
	unsigned int angle_resolution;
	unsigned int refinement_type;
	float maxResolutionForIndexing_1_per_A;
	float tolerance;
	int multi;
	int thread_count;
	int min_peaks;
	int no_check_indexed;
	float reflectionRadius; /* In m^-1 */
	float customPhotonEnergy;
	float customBandwidth;
	float maxRefinementDisbalance;
};

#ifdef HAVE_PINKINDEXER

#include <pinkIndexer/adaptions/crystfel/Lattice.h>
#include <pinkIndexer/adaptions/crystfel/ExperimentSettings.h>
#include <pinkIndexer/adaptions/crystfel/PinkIndexer.h>

#define MAX_MULTI_LATTICE_COUNT 8

struct pinkIndexer_private_data {
	PinkIndexer *pinkIndexer;
	reciprocalPeaks_1_per_A_t reciprocalPeaks_1_per_A;
	float *intensities;

	IndexingMethod indm;
	UnitCell *cellTemplate;
	int threadCount;
	int multi;
	int min_peaks;

	int no_check_indexed;

	float maxRefinementDisbalance;

	IntegerMatrix *centeringTransformation;
	LatticeTransform_t latticeReductionTransform;
};

//static void reduceCell(UnitCell* cell, LatticeTransform_t* appliedReductionTransform);
//static void restoreCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform);
static void reduceReciprocalCell(UnitCell* cell, LatticeTransform_t* appliedReductionTransform);
static void restoreReciprocalCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform);
static void makeRightHanded(UnitCell* cell);

int run_pinkIndexer(struct image *image, void *ipriv)
{
	struct pinkIndexer_private_data* pinkIndexer_private_data = (struct pinkIndexer_private_data*) ipriv;
	reciprocalPeaks_1_per_A_t* reciprocalPeaks_1_per_A = &(pinkIndexer_private_data->reciprocalPeaks_1_per_A);
	float *intensities = pinkIndexer_private_data->intensities;

	int peakCountMax = image_feature_count(image->features);
	if (peakCountMax < 5) {
		int goodLatticesCount = 0;
		return goodLatticesCount;
	}
	reciprocalPeaks_1_per_A->peakCount = 0;
	for (int i = 0; i < peakCountMax && i < MAX_PEAK_COUNT_FOR_INDEXER; i++) {
		struct imagefeature *f;
		f = image_get_feature(image->features, i);
		if (f == NULL) {
			continue;
		}

		reciprocalPeaks_1_per_A->coordinates_x[reciprocalPeaks_1_per_A->peakCount] = f->rz * 1e-10;
		reciprocalPeaks_1_per_A->coordinates_y[reciprocalPeaks_1_per_A->peakCount] = f->rx * 1e-10;
		reciprocalPeaks_1_per_A->coordinates_z[reciprocalPeaks_1_per_A->peakCount] = f->ry * 1e-10;
		intensities[reciprocalPeaks_1_per_A->peakCount] = (float) (f->intensity);
		reciprocalPeaks_1_per_A->peakCount++;
	}
	int indexed = 0;
	Lattice_t indexedLattice[MAX_MULTI_LATTICE_COUNT];
	float center_shift[MAX_MULTI_LATTICE_COUNT][2];



	do {
		int peakCount = reciprocalPeaks_1_per_A->peakCount;
		int matchedPeaksCount = PinkIndexer_indexPattern(pinkIndexer_private_data->pinkIndexer,
		        &(indexedLattice[indexed]), center_shift[indexed], reciprocalPeaks_1_per_A, intensities,
		        pinkIndexer_private_data->maxRefinementDisbalance,
		        pinkIndexer_private_data->threadCount);

		if(matchedPeaksCount == -1){
			STATUS("WARNING: Indexing solution was rejected due to too large disbalance of the refinement."
			        "If you see this message often, check the documentation for the parameter "
			        "--pinkIndexer-max-refinement-disbalance\n");

			matchedPeaksCount = 0;
		}

		printf("matchedPeaksCount %d from %d\n",matchedPeaksCount,peakCount);
		if ((matchedPeaksCount >= 25 && matchedPeaksCount >= peakCount * 0.30)
		        || matchedPeaksCount >= peakCount * 0.4
		        || matchedPeaksCount >= 70
		        || pinkIndexer_private_data->no_check_indexed == 1)
		                {
			UnitCell *uc;
			uc = cell_new();

			Lattice_t *l = &(indexedLattice[indexed]);

			cell_set_reciprocal(uc, l->ay * 1e10, l->az * 1e10, l->ax * 1e10,
			        l->by * 1e10, l->bz * 1e10, l->bx * 1e10,
			        l->cy * 1e10, l->cz * 1e10, l->cx * 1e10);

			restoreReciprocalCell(uc, &pinkIndexer_private_data->latticeReductionTransform);

			UnitCell *new_cell_trans = cell_transform_intmat(uc, pinkIndexer_private_data->centeringTransformation);
			cell_free(uc);
			uc = new_cell_trans;

			cell_set_lattice_type(new_cell_trans, cell_get_lattice_type(pinkIndexer_private_data->cellTemplate));
			cell_set_centering(new_cell_trans, cell_get_centering(pinkIndexer_private_data->cellTemplate));
			cell_set_unique_axis(new_cell_trans, cell_get_unique_axis(pinkIndexer_private_data->cellTemplate));

			if (validate_cell(uc)) {
				ERROR("pinkIndexer: problem with returned cell!\n");
			}

			Crystal * cr = crystal_new();
			if (cr == NULL) {
				ERROR("Failed to allocate crystal.\n");
				return 0;
			}
			crystal_set_cell(cr, uc);
			crystal_set_det_shift(cr, center_shift[indexed][0], center_shift[indexed][1]);
			image_add_crystal(image, cr);
			indexed++;

		} else {
			break;
		}
	} while (pinkIndexer_private_data->multi
	        && indexed <= MAX_MULTI_LATTICE_COUNT
	        && reciprocalPeaks_1_per_A->peakCount >= pinkIndexer_private_data->min_peaks);

	return indexed;
}

void *pinkIndexer_prepare(IndexingMethod *indm, UnitCell *cell,
                          struct pinkIndexer_options *pinkIndexer_opts,
                          const DataTemplate *dtempl)
{
	if ( beam->photon_energy_from != NULL && pinkIndexer_opts->customPhotonEnergy <= 0) {
		ERROR("For pinkIndexer, the photon_energy must be defined as a "
		      "constant in the geometry file or as a parameter (see --pinkIndexer-override-photon-energy)\n");
		return NULL;
	}
	if ( (det->panels[0].clen_from != NULL) && pinkIndexer_opts->refinement_type ==
	        REFINEMENT_TYPE_firstFixedThenVariableLatticeParametersCenterAdjustmentMultiSeed) {
		ERROR("Using center refinement makes it necessary to have the detector distance fixed in the geometry file!");
		return NULL;
	}
	if(cell == NULL){
		ERROR("Pink indexer needs a unit cell file to be specified!");
		return NULL;
	}

	struct pinkIndexer_private_data* pinkIndexer_private_data = malloc(sizeof(struct pinkIndexer_private_data));
	allocReciprocalPeaks(&(pinkIndexer_private_data->reciprocalPeaks_1_per_A));
	pinkIndexer_private_data->intensities = malloc(MAX_PEAK_COUNT_FOR_INDEXER * sizeof(float));
	pinkIndexer_private_data->indm = *indm;
	pinkIndexer_private_data->cellTemplate = cell;
	pinkIndexer_private_data->threadCount = pinkIndexer_opts->thread_count;
	pinkIndexer_private_data->multi = pinkIndexer_opts->multi;
	pinkIndexer_private_data->min_peaks = pinkIndexer_opts->min_peaks;
	pinkIndexer_private_data->no_check_indexed = pinkIndexer_opts->no_check_indexed;
	pinkIndexer_private_data->maxRefinementDisbalance = pinkIndexer_opts->maxRefinementDisbalance;

	UnitCell* primitiveCell = uncenter_cell(cell, &pinkIndexer_private_data->centeringTransformation, NULL);

	//reduceCell(primitiveCell, &pinkIndexer_private_data->latticeReductionTransform);
	reduceReciprocalCell(primitiveCell, &pinkIndexer_private_data->latticeReductionTransform);

	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	int ret = cell_get_reciprocal(primitiveCell, &asx, &asy, &asz, &bsx, &bsy, &bsz, &csx, &csy, &csz);
	if (ret != 0) {
		ERROR("cell_get_reciprocal did not finish properly!");
	}

	Lattice_t lattice = { .ax = asz * 1e-10, .ay = asx * 1e-10, .az = asy * 1e-10,
	        .bx = bsz * 1e-10, .by = bsx * 1e-10, .bz = bsy * 1e-10,
	        .cx = csz * 1e-10, .cy = csx * 1e-10, .cz = csy * 1e-10 };

	float detectorDistance_m;
	if ( det->panels[0].clen_from != NULL ) {
		detectorDistance_m =  0.25;  /* fake value */
	} else {
		detectorDistance_m = det->panels[0].clen + det->panels[0].coffset;
	}

	float beamEenergy_eV = beam->photon_energy;
	float nonMonochromaticity = beam->bandwidth*5;
	if(pinkIndexer_opts->customPhotonEnergy > 0){
		beamEenergy_eV = pinkIndexer_opts->customPhotonEnergy;
	}
	if(pinkIndexer_opts->customBandwidth >= 0){
		nonMonochromaticity = pinkIndexer_opts->customBandwidth;
	}

	float reflectionRadius_1_per_A;
	if (pinkIndexer_opts->reflectionRadius < 0) {
		reflectionRadius_1_per_A = 0.02
		        * sqrt(lattice.ax * lattice.ax + lattice.ay * lattice.ay + lattice.az * lattice.az);
	}
	else {
		reflectionRadius_1_per_A = pinkIndexer_opts->reflectionRadius * 1e10;  /* m^-1 to A^-1*/
	}

	if(beamEenergy_eV > 75000 && nonMonochromaticity < 0.02 && reflectionRadius_1_per_A < 0.0005){
		STATUS("Trying to index electron diffraction? It might be helpful to set a higher reflection radius (see documentation for --pinkIndexer-reflection-radius)")
	}

	float divergenceAngle_deg = 0.01; //fake

	float tolerance = pinkIndexer_opts->tolerance;
	Lattice_t sampleReciprocalLattice_1_per_A = lattice;
	float detectorRadius_m = 0.03; //fake, only for prediction
	ExperimentSettings* experimentSettings = ExperimentSettings_new(beamEenergy_eV, detectorDistance_m,
	        detectorRadius_m, divergenceAngle_deg, nonMonochromaticity, sampleReciprocalLattice_1_per_A, tolerance,
	        reflectionRadius_1_per_A);

	consideredPeaksCount_t consideredPeaksCount = pinkIndexer_opts->considered_peaks_count;
	angleResolution_t angleResolution = pinkIndexer_opts->angle_resolution;
	refinementType_t refinementType = pinkIndexer_opts->refinement_type;
	float maxResolutionForIndexing_1_per_A = pinkIndexer_opts->maxResolutionForIndexing_1_per_A;
	pinkIndexer_private_data->pinkIndexer = PinkIndexer_new(experimentSettings, consideredPeaksCount, angleResolution,
	        refinementType,
	        maxResolutionForIndexing_1_per_A);

	ExperimentSettings_delete(experimentSettings);
	cell_free(primitiveCell);

	/* Flags that pinkIndexer knows about */
	*indm &= INDEXING_METHOD_MASK
	        | INDEXING_USE_CELL_PARAMETERS;

	return pinkIndexer_private_data;
}

//static void reduceCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform)
//{
//	double ax, ay, az, bx, by, bz, cx, cy, cz;
//	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
//
//	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };
//
//	reduceLattice(&l, appliedReductionTransform);
//
//	cell_set_cartesian(cell, l.ax, l.ay, l.az,
//	        l.bx, l.by, l.bz,
//	        l.cx, l.cy, l.cz);
//
//	makeRightHanded(cell);
//}
//
//static void restoreCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform)
//{
//
//	double ax, ay, az, bx, by, bz, cx, cy, cz;
//	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
//
//	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };
//
//	restoreLattice(&l, appliedReductionTransform);
//
//	cell_set_cartesian(cell, l.ax, l.ay, l.az,
//	        l.bx, l.by, l.bz,
//	        l.cx, l.cy, l.cz);
//
//	makeRightHanded(cell);
//}

static void reduceReciprocalCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };

	reduceLattice(&l, appliedReductionTransform);

	cell_set_reciprocal(cell, l.ax, l.ay, l.az,
	        l.bx, l.by, l.bz,
	        l.cx, l.cy, l.cz);

	makeRightHanded(cell);
}

static void restoreReciprocalCell(UnitCell *cell, LatticeTransform_t* appliedReductionTransform)
{

	double ax, ay, az, bx, by, bz, cx, cy, cz;
	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	Lattice_t l = { ax, ay, az, bx, by, bz, cx, cy, cz };

	restoreLattice(&l, appliedReductionTransform);

	cell_set_reciprocal(cell, l.ax, l.ay, l.az,
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

void pinkIndexer_cleanup(void *pp)
{
	struct pinkIndexer_private_data* pinkIndexer_private_data = (struct pinkIndexer_private_data*) pp;

	freeReciprocalPeaks(pinkIndexer_private_data->reciprocalPeaks_1_per_A);
	free(pinkIndexer_private_data->intensities);
	intmat_free(pinkIndexer_private_data->centeringTransformation);
	PinkIndexer_delete(pinkIndexer_private_data->pinkIndexer);
}

const char *pinkIndexer_probe(UnitCell *cell)
{
	return "pinkIndexer";
}

#else /* HAVE_PINKINDEXER */

int run_pinkIndexer(struct image *image, void *ipriv)
{
	ERROR("This copy of CrystFEL was compiled without PINKINDEXER support.\n");
	return 0;
}

extern void *pinkIndexer_prepare(IndexingMethod *indm, UnitCell *cell,
		struct pinkIndexer_options *pinkIndexer_opts,
		const DataTemplate *dtempl)
{
	ERROR("This copy of CrystFEL was compiled without PINKINDEXER support.\n");
	ERROR("To use PINKINDEXER indexing, recompile with PINKINDEXER.\n");
	return NULL;
}

void pinkIndexer_cleanup(void *pp)
{
}

const char *pinkIndexer_probe(UnitCell *cell)
{
	return NULL;
}

#endif /* HAVE_PINKINDEXER */

static void pinkIndexer_show_help()
{
	printf(
"Parameters for the PinkIndexer indexing algorithm:\n"
"     --pinkIndexer-considered-peaks-count=n\n"
"                           Considered peaks count, 0 (fewest) to 4 (most)\n"
"                            Default: 4\n"
"     --pinkIndexer-angle-resolution=n\n"
"                           Angle resolution, 0 (loosest) to 4 (most dense)\n"
"                            Default: 2\n"
"     --pinkIndexer-refinement-type=n\n"
"                           Refinement type, 0 (none) to 5 (most accurate)\n"
"                            Default: 1\n"
"     --pinkIndexer-tolerance=n\n"
"                           Relative tolerance of the lattice vectors.\n"
"                            Default 0.06\n"
"     --pinkIndexer-reflection-radius=n\n"
"                           Radius of the reflections in reciprocal space.\n"
"                            Specified in 1/A.  Default is 2%% of a*.\n"
"     --pinkIndexer-max-resolution-for-indexing=n\n"
"                           Measured in 1/A\n"
"     --pinkIndexer-multi   Use pinkIndexers own multi indexing.\n"
"     --pinkIndexer-thread-count=n\n"
"                           Thread count for internal parallelization \n"
"                            Default: 1\n"
"     --pinkIndexer-no-check-indexed\n"
"                           Disable internal check for correct indexing\n"
"                            solutions\n"
"     --pinkIndexer-max-refinement-disbalance=n\n"
"                           Maximum disbalance after refinement:\n"
"                            0 (no disbalance) to 2 (extreme disbalance), default 0.4\n"
"     --pinkIndexer-override-photon-energy=ev\n"
"                           Mean energy in eV to use for indexing.\n"
"     --pinkIndexer-override-bandwidth=n\n"
"                           Bandwidth in (delta energy)/(mean energy) to use for indexing.\n"
"     --pinkIndexer-override-visible-energy-range=min-max\n"
"                           Overrides photon energy and bandwidth according to a range of \n"
"                           energies that have high enough intensity to produce \"visible\" \n"
"                           Bragg spots on the detector.\n"
"                           Min and max range borders are separated by a minus sign (no whitespace).\n"
	);
}


static error_t pinkindexer_parse_arg(int key, char *arg,
                                     struct argp_state *state)
{
	float tmp, tmp2;
	struct pinkIndexer_options **opts_ptr = state->input;

	switch ( key ) {

		case ARGP_KEY_INIT :
		*opts_ptr = malloc(sizeof(struct pinkIndexer_options));
		if ( *opts_ptr == NULL ) return ENOMEM;
		(*opts_ptr)->considered_peaks_count = 4;
		(*opts_ptr)->angle_resolution = 2;
		(*opts_ptr)->refinement_type = 1;
		(*opts_ptr)->tolerance = 0.06;
		(*opts_ptr)->maxResolutionForIndexing_1_per_A = +INFINITY;
		(*opts_ptr)->thread_count = 1;
		(*opts_ptr)->multi = 0;
		(*opts_ptr)->no_check_indexed = 0;
		(*opts_ptr)->min_peaks = 2;
		(*opts_ptr)->reflectionRadius = -1;
		(*opts_ptr)->customPhotonEnergy = -1;
		(*opts_ptr)->customBandwidth = -1;
		(*opts_ptr)->maxRefinementDisbalance = 0.4;
		break;

		case 1 :
		pinkIndexer_show_help();
		return EINVAL;

		case 2 :
		if (sscanf(arg, "%u", &(*opts_ptr)->considered_peaks_count) != 1)
		{
			ERROR("Invalid value for "
			      "--pinkIndexer-considered-peaks-count\n");
			return EINVAL;
		}
		break;

		case 3 :
		if (sscanf(arg, "%u", &(*opts_ptr)->angle_resolution) != 1)
		{
			ERROR("Invalid value for "
			      "--pinkIndexer-angle_resolution\n");
			return EINVAL;
		}
		break;

		case 4 :
		if (sscanf(arg, "%u", &(*opts_ptr)->refinement_type) != 1)
		{
			ERROR("Invalid value for "
			      "--pinkIndexer-refinement-type\n");
			return EINVAL;
		}
		break;

		case 5 :
		if (sscanf(arg, "%d", &(*opts_ptr)->thread_count) != 1)
		{
			ERROR("Invalid value for --pinkIndexer-thread-count\n");
			return EINVAL;
		}
		break;

		case 6 :
		if (sscanf(arg, "%f", &(*opts_ptr)->maxResolutionForIndexing_1_per_A) != 1)
		{
			ERROR("Invalid value for "
			      "--pinkIndexer-max-resolution-for-indexing\n");
			return EINVAL;
		}
		break;

		case 7 :
		if (sscanf(arg, "%f", &(*opts_ptr)->tolerance) != 1)
		{
			ERROR("Invalid value for --pinkIndexer-tolerance\n");
			return EINVAL;
		}
		break;

		case 8 :
		(*opts_ptr)->multi = 1;
		break;

		case 9 :
		(*opts_ptr)->no_check_indexed = 1;
		break;

		case 10 :
		if (sscanf(arg, "%f", &tmp) != 1) {
			ERROR("Invalid value for --pinkIndexer-reflection-radius\n");
			return EINVAL;
		}
		(*opts_ptr)->reflectionRadius = tmp / 1e10; /* A^-1 to m^-1 */
		break;

		case 11 :
		if (sscanf(arg, "%f", &(*opts_ptr)->customPhotonEnergy) != 1)
		{
			ERROR("Invalid value for --pinkIndexer-override-photon-energy\n");
			return EINVAL;
		}
		break;

		case 12 :
		if (sscanf(arg, "%f", &(*opts_ptr)->customBandwidth) != 1)
		{
			ERROR("Invalid value for --pinkIndexer-override-bandwidth\n");
			return EINVAL;
		}
		break;
		case 13 :
		if (sscanf(arg, "%f-%f", &tmp, &tmp2) != 2)
		{
			ERROR("Invalid value for --pinkIndexer-override-visible-energy-range\n");
			return EINVAL;
		}
		(*opts_ptr)->customPhotonEnergy = (tmp + tmp2)/2;
		(*opts_ptr)->customBandwidth = (tmp2 - tmp)/(*opts_ptr)->customPhotonEnergy;
		if((*opts_ptr)->customBandwidth < 0){
			(*opts_ptr)->customBandwidth *= -1;
		}
		break;
		case 14 :
		if (sscanf(arg, "%f", &(*opts_ptr)->maxRefinementDisbalance) != 1)
		{
			ERROR("Invalid value for --pinkIndexer-max-refinement-disbalance\n");
			return EINVAL;
		}
	}

	return 0;
}


static struct argp_option pinkindexer_options[] = {

	{"help-pinkindexer", 1, NULL, OPTION_NO_USAGE,
	 "Show options for PinkIndexer indexing algorithm", 99},

	{"pinkIndexer-considered-peaks-count", 2, "n", OPTION_HIDDEN, NULL},
	{"pinkIndexer-cpc", 2, "n", OPTION_HIDDEN, NULL},

	{"pinkIndexer-angle-resolution", 3, "ang", OPTION_HIDDEN, NULL},
	{"pinkIndexer-ar", 3, "ang", OPTION_HIDDEN, NULL},

	{"pinkIndexer-refinement-type", 4, "t", OPTION_HIDDEN, NULL},
	{"pinkIndexer-rt", 4, "t", OPTION_HIDDEN, NULL},

	{"pinkIndexer-thread-count", 5, "n", OPTION_HIDDEN, NULL},
	{"pinkIndexer-tc", 5, "n", OPTION_HIDDEN, NULL},

	{"pinkIndexer-max-resolution-for-indexing", 6, "res", OPTION_HIDDEN, NULL},
	{"pinkIndexer-mrfi", 6, "res", OPTION_HIDDEN, NULL},

	{"pinkIndexer-tolerance", 7, "tol", OPTION_HIDDEN, NULL},
	{"pinkIndexer-tol", 7, "tol", OPTION_HIDDEN, NULL},

	{"pinkIndexer-multi", 8, NULL, OPTION_HIDDEN, NULL},

	{"pinkIndexer-no-check-indexed", 9, NULL, OPTION_HIDDEN, NULL},

	{"pinkIndexer-reflection-radius", 10, "r", OPTION_HIDDEN, NULL},

	{"pinkIndexer-override-photon-energy", 11, "ev", OPTION_HIDDEN, NULL},

	{"pinkIndexer-override-bandwidth", 12, "bw", OPTION_HIDDEN, NULL},

	{"pinkIndexer-override-visible-energy-range", 13, "overridenVisibleEnergyRange", OPTION_HIDDEN, NULL},

	{"pinkIndexer-max-refinement-disbalance", 14, "maxDisbalance", OPTION_HIDDEN, NULL},

	{0}
};


struct argp pinkIndexer_argp = { pinkindexer_options,
                                 pinkindexer_parse_arg,
                                 NULL, NULL, NULL, NULL, NULL };
