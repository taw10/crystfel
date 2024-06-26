/*
 * ffbidx.c
 *
 * Invoke the Fast Feedback Indexer library
 *
 * Copyright © 2023 Paul Scherrer Institute
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Filip Leonarski <filip.leonarski@psi.ch>
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

#include <libcrystfel-config.h>
#include <stdio.h>
#include <stdlib.h>

#include "cell-utils.h"
#include "ffbidx.h"

#ifdef HAVE_FFBIDX

#include <ffbidx/c_api.h>

struct ffbidx_private_data {
    UnitCell *cellTemplate;
    struct ffbidx_options opts;
};

static void makeRightHanded(UnitCell *cell)
{
    // From xgandalf.c
    double ax, ay, az, bx, by, bz, cx, cy, cz;
    cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

    if ( !right_handed(cell) ) {
        cell_set_cartesian(cell, -ax, -ay, -az, -bx, -by, -bz, -cx, -cy, -cz);
    }
}

int run_ffbidx(struct image *image, void *ipriv) {
    int npk;
    int i;

    struct ffbidx_private_data *prv_data = (struct ffbidx_private_data *) ipriv;

    npk = image_feature_count(image->features);
    if ( npk < prv_data->opts.min_peaks )
        return 0;

    // Setup configuration
    struct config_runtime cruntime;
    struct config_persistent cpers;
    struct config_ifssr cifssr;

    set_defaults(&cpers, &cruntime, &cifssr);
    cpers.max_spots = prv_data->opts.max_peaks;

    cifssr.min_spots = prv_data->opts.min_peaks;
    cpers.max_output_cells = prv_data->opts.output_cells;
    cpers.num_candidate_vectors = prv_data->opts.num_candidate_vectors;

    struct input ffbidx_input;
    struct output ffbidx_output;

    ffbidx_input.cell.x = (float *) calloc(3, sizeof(float));
    ffbidx_input.cell.y = (float *) calloc(3, sizeof(float));
    ffbidx_input.cell.z = (float *) calloc(3, sizeof(float));

    ffbidx_input.spot.x = (float *) calloc(npk, sizeof(float));
    ffbidx_input.spot.y = (float *) calloc(npk, sizeof(float));
    ffbidx_input.spot.z = (float *) calloc(npk, sizeof(float));

    ffbidx_input.n_cells = 1;
    ffbidx_input.n_spots = npk;
    ffbidx_input.new_cells = true;
    ffbidx_input.new_spots = true;

    ffbidx_output.x = (float *) calloc(3 * cpers.max_output_cells, sizeof(float));
    ffbidx_output.y = (float *) calloc(3 * cpers.max_output_cells, sizeof(float));
    ffbidx_output.z = (float *) calloc(3 * cpers.max_output_cells, sizeof(float));
    ffbidx_output.score = (float *) calloc(cpers.max_output_cells, sizeof(float));
    ffbidx_output.n_cells = cpers.max_output_cells;

    for ( i=0; i<npk; i++ ) {

        struct imagefeature *f;
        double r[3];

        f = image_get_feature(image->features, i);
        if ( f == NULL ) {
            ERROR("Empty feature ???");
            continue;
        }

        detgeom_transform_coords(&image->detgeom->panels[f->pn],
                                 f->fs, f->ss, image->lambda,
                                 0.0, 0.0, r);
        ffbidx_input.spot.x[i] = r[0] * 1e-10;
        ffbidx_input.spot.y[i] = r[1] * 1e-10;
        ffbidx_input.spot.z[i] = r[2] * 1e-10;
    }

    double cell_internal_double[9];

    cell_get_cartesian(prv_data->cellTemplate,
                       &cell_internal_double[0],&cell_internal_double[3],&cell_internal_double[6],
                       &cell_internal_double[1],&cell_internal_double[4],&cell_internal_double[7],
                       &cell_internal_double[2],&cell_internal_double[5],&cell_internal_double[8]);

    for (int i = 0; i < 3; i++) {
        ffbidx_input.cell.x[i] = cell_internal_double[0 + i] * 1e10;
        ffbidx_input.cell.y[i] = cell_internal_double[3 + i] * 1e10;
        ffbidx_input.cell.z[i] = cell_internal_double[6 + i] * 1e10;

//        ffbidx_input.cell.x[i] = cell_internal_double[3*i] * 1e10;
//        ffbidx_input.cell.y[i] = cell_internal_double[3*i+1] * 1e10;
//        ffbidx_input.cell.z[i] = cell_internal_double[3*i+2] * 1e10;
    }

    // settings.cvc_threshold = prv_data->opts.threshold_for_solution;
    // settings.crt_num_sample_points = prv_data->opts.sample_points;

    int handle = create_indexer(&cpers, NULL, NULL);
    if (handle < 0) {
        ERROR("Error creating indexer\n");
        return 0;
        // handle error
    }

    if (index_start(handle, &ffbidx_input, &ffbidx_output, &cruntime, NULL, NULL)) {
        ERROR("Error running index_start\n");
        drop_indexer(handle);
        return 0;
    }

    if (index_end(handle, &ffbidx_output)) {
        ERROR("Error running index_end\n");
        drop_indexer(handle);
        return 0;
    }

    if (refine(handle, &ffbidx_input, &ffbidx_output, &cifssr, 0, 1)) {
        ERROR("Error running refine\n");
        drop_indexer(handle);
        return 0;
    }

    int bcell;
    if ((bcell = best_cell(handle, &ffbidx_output)) < 0) {
        ERROR("Error running best_cell\n");
        drop_indexer(handle);
        return 0;
    }
    UnitCell *uc;

    uc = cell_new();

    cell_set_cartesian(uc,
                        ffbidx_output.x[bcell * 3 + 0] * 1e-10, ffbidx_output.y[bcell * 3 + 0] * 1e-10, ffbidx_output.z[bcell * 3 + 0] * 1e-10,
                       ffbidx_output.x[bcell * 3 + 1] * 1e-10, ffbidx_output.y[bcell * 3 + 1] * 1e-10, ffbidx_output.z[bcell * 3 + 1] * 1e-10,
                       ffbidx_output.x[bcell * 3 + 2] * 1e-10, ffbidx_output.y[bcell * 3 + 2] * 1e-10, ffbidx_output.z[bcell * 3 + 2] * 1e-10);

    makeRightHanded(uc);

    cell_set_lattice_type(uc, cell_get_lattice_type(prv_data->cellTemplate));
    cell_set_centering(uc, cell_get_centering(prv_data->cellTemplate));
    cell_set_unique_axis(uc, cell_get_unique_axis(prv_data->cellTemplate));

    free(ffbidx_input.spot.x);
    free(ffbidx_input.spot.y);
    free(ffbidx_input.spot.z);

    free(ffbidx_input.cell.x);
    free(ffbidx_input.cell.y);
    free(ffbidx_input.cell.z);

    free(ffbidx_output.x);
    free(ffbidx_output.y);
    free(ffbidx_output.z);
    free(ffbidx_output.score);

    if (drop_indexer(handle) != 0) {
        ERROR("Error dropping indexer\n");
    }

    if ( validate_cell(uc) ) {
        ERROR("ffbidx: problem with returned cell!\n");
        cell_free(uc);
        return 0;
    }

    Crystal *cr = crystal_new();
    if ( cr == NULL ) {
        ERROR("Failed to allocate crystal.\n");
        return 0;
    }
    crystal_set_cell(cr, uc);
    crystal_set_det_shift(cr, 0, 0);
    image_add_crystal(image, cr);
    return 1;
}

void *ffbidx_prepare(IndexingMethod *indm, UnitCell *cell, struct ffbidx_options *opts) {
    if ( cell == NULL ) {
        ERROR("Unit cell information is required for fast feedback indexer.\n");
        return NULL;
    }

    struct ffbidx_private_data *prv_data = (struct ffbidx_private_data *) malloc(sizeof(struct ffbidx_private_data));

    prv_data->cellTemplate = cell;
    prv_data->opts = *opts;

    *indm &= INDEXING_METHOD_MASK | INDEXING_USE_CELL_PARAMETERS;
    return prv_data;
}

void ffbidx_cleanup(void *pp) {
    free(pp);
}

const char *ffbidx_probe(UnitCell *cell) {
    return "ffbidx";
}

#else

int run_ffbidx(struct image *image, void *ipriv)
{
	ERROR("This copy of CrystFEL was compiled without FFBIDX support.\n");
	return 0;
}


void *ffbidx_prepare(IndexingMethod *indm, UnitCell *cell, struct ffbidx_options *opts)
{
	ERROR("This copy of CrystFEL was compiled without FFBIDX support.\n");
	ERROR("To use FFBIDX indexing, recompile with FFBIDX.\n");
	return NULL;
}


void ffbidx_cleanup(void *pp)
{
}


const char *ffbidx_probe(UnitCell *cell)
{
	return NULL;
}

#endif

static void ffbidx_show_help()
{
    printf("Parameters for the fast feedback indexing algorithm:\n"
           "     --ffbidx-max-peaks\n"
           "                            Maximum number of peaks used for indexing.\n"
           "                            Default: 100\n"
           "     --ffbidx-min-peaks\n"
           "                            Maximum number of indexed peaks to accept solution.\n"
           "                            Default: 9\n"
           "     --ffbidx-threshold\n"
           "                            Threshold to accept solution as indexed.\n"
           "                            Default: 0.02\n"
           "     --ffbidx-output-cells\n"
           "                            Number of output cells.\n"
           "                            Default: 1\n"
           "     --ffbidx-sample-points\n"
           "                            Number of sample points.\n"
           "                            Default: 32768\n"

    );
}


int ffbidx_default_options(struct ffbidx_options **opts_ptr)
{
    struct ffbidx_options *opts;

    opts = malloc(sizeof(struct ffbidx_options));
    if ( opts == NULL ) return ENOMEM;

    opts->max_peaks = 100;
    opts->min_peaks = 9;
    opts->threshold_for_solution = 0.02f;
    opts->output_cells = 1;
    opts->sample_points = 32*1024;
    opts->num_candidate_vectors = 32;
    *opts_ptr = opts;
    return 0;
}


static error_t ffbidx_parse_arg(int key, char *arg, struct argp_state *state)
{
    struct ffbidx_options **opts_ptr = state->input;
    int r;

    switch ( key ) {
        case ARGP_KEY_INIT :
            r = ffbidx_default_options(opts_ptr);
            if ( r ) return r;
            break;

        case 1 :
            ffbidx_show_help();
            return EINVAL;

        case 2 :
            if (sscanf(arg, "%u", &(*opts_ptr)->max_peaks) != 1) {
                ERROR("Invalid value for --ffbidx-max-peaks\n");
                return EINVAL;
            }
            break;

        case 3 :
            if (sscanf(arg, "%u", &(*opts_ptr)->min_peaks) != 1) {
                ERROR("Invalid value for --ffbidx-min-peaks\n");
                return EINVAL;
            }
            break;

        case 4 :
            if (sscanf(arg, "%f", &(*opts_ptr)->threshold_for_solution) != 1) {
                ERROR("Invalid value for --ffbidx-threshold\n");
                return EINVAL;
            }
            if (((*opts_ptr)->threshold_for_solution <= 0.0f) || ((*opts_ptr)->threshold_for_solution > 1.0f)) {
                ERROR("Invalid value for --ffbidx-threshold; must be in range 0.0-1.0\n");
                return EINVAL;
            }
            break;
        case 5 :
            if (sscanf(arg, "%u", &(*opts_ptr)->output_cells) != 1) {
                ERROR("Invalid value for --ffbidx-output-cells\n");
                return EINVAL;
            }
            if (((*opts_ptr)->output_cells == 0) || ((*opts_ptr)->output_cells > 128)) {
                ERROR("Invalid value for --ffbidx-output-cells; must be in range 1-128\n");
                return EINVAL;
            }
            break;
        case 6 :
            if (sscanf(arg, "%u", &(*opts_ptr)->sample_points) != 1) {
                ERROR("Invalid value for --ffbidx-sample-points\n");
                return EINVAL;
            }
            break;
    }

    return 0;
}


static struct argp_option ffbidx_options[] = {
        {"help-ffbidx", 1, NULL, OPTION_NO_USAGE, "Show options for fast feedback indexing algorithm", 99},
        {"ffbidx-max-peaks", 2, "ffbidx_maxn", OPTION_HIDDEN, NULL},
        {"ffbidx-min-peaks", 3, "ffbidx_minn", OPTION_HIDDEN, NULL},
        {"ffbidx-threshold", 4, "ffbidx_threshold", OPTION_HIDDEN, NULL},
        {"ffbidx-output-cells", 5, "ffbidx_out_cells", OPTION_HIDDEN, NULL},
        {"ffbidx-sample-points", 6, "ffbidx_sample_points", OPTION_HIDDEN, NULL},
        {0}
};


struct argp ffbidx_argp = { ffbidx_options, ffbidx_parse_arg,
                              NULL, NULL, NULL, NULL, NULL };
