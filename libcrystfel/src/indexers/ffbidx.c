/*
 * ffbidx.c
 *
 * Invoke the Fast Feedback Indexer library
 *
 * Copyright © 2023-2024 Paul Scherrer Institute
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023-2024 Filip Leonarski <filip.leonarski@psi.ch>
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

    ffbidx_output.x = (float *) calloc(3 * prv_data->opts.cpers.max_output_cells, sizeof(float));
    ffbidx_output.y = (float *) calloc(3 * prv_data->opts.cpers.max_output_cells, sizeof(float));
    ffbidx_output.z = (float *) calloc(3 * prv_data->opts.cpers.max_output_cells, sizeof(float));
    ffbidx_output.score = (float *) calloc(prv_data->opts.cpers.max_output_cells, sizeof(float));
    ffbidx_output.n_cells = prv_data->opts.cpers.max_output_cells;

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
    }

    int handle = create_indexer(&prv_data->opts.cpers, NULL, NULL);
    if (handle < 0) {
        ERROR("Error creating indexer\n");
        return 0;
        // handle error
    }

    if (index_start(handle, &ffbidx_input, &ffbidx_output, &prv_data->opts.cruntime, NULL, NULL)) {
        ERROR("Error running index_start\n");
        drop_indexer(handle);
        return 0;
    }

    if (index_end(handle, &ffbidx_output)) {
        ERROR("Error running index_end\n");
        drop_indexer(handle);
        return 0;
    }

    if (refine(handle, &ffbidx_input, &ffbidx_output, &prv_data->opts.cifssr, 0, 1)) {
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

static void ffbidx_show_help(void)
{
    struct config_runtime r;
    struct config_persistent p;
    struct config_ifssr i;
    set_defaults(&p, &r, &i);
    printf("Parameters for the fast feedback indexing algorithm:\n");;
    printf("     --ffbidx-max-spots=\n");
    printf("                            Maximum number of peaks used for indexing.\n");
    printf("                            Default: %d\n", p.max_spots);
    printf("     --ffbidx-output-cells=\n");
    printf("                            Number of output cells.\n");
    printf("                            Default: %d\n", p.max_output_cells);
    printf("     --ffbidx-num-candidate-vectors=\n");
    printf("                            Number of candidate vectors (per input cell vector).\n");
    printf("                            Default: %d\n", p.num_candidate_vectors);
    printf("     --ffbidx-redundant-calculations\n");
    printf("     --ffbidx-no-redundant-calculations\n");
    printf("                            Compute candidates for all three cell vectors instead of just one\n");
    printf("                            Default: %s\n", p.redundant_computations
                                                               ? "--ffbidx-redundant-calculations"
                                                               : "--ffbidx-no-redundant-calculations");
    printf("     --ffbidx-length-threshold=\n");
    printf("                            Threshold for determining equal vector length (|va| - threshold < |vb| < |va| + threshold)\n");
    printf("                            Default: %f\n", r.length_threshold);

    printf("     --ffbidx-triml=\n");
    printf("                            Lower trim value for distance to nearest integer objective value - 0 < triml < trimh\n");
    printf("                            Default: %f\n", r.triml);

    printf("     --ffbidx-trimh=\n");
    printf("                            Higher trim value for distance to nearest integer objective value - triml < trimh < 0.5\n");
    printf("                            Default: %f\n", r.trimh);
    printf("     --ffbidx-delta=\n");
    printf("                            Log2 curve position: score = log2(trim(dist(x)) + delta)\n");
    printf("                            Default: %f\n", r.delta);
    printf("     --ffbidx-dist1=\n");
    printf("                            Maximum distance to int for single coordinate\n");
    printf("                            Default: %f\n", r.dist1);
    printf("     --ffbidx-dist3=\n");
    printf("                            Maximum distance to int for triple coordinates\n");
    printf("                            Default: %f\n", r.dist3);
    printf("     --ffbidx-num-halfsphere-points=\n");
    printf("                            Number of sample points on half sphere for finding vector candidates\n");
    printf("                            Default: %d\n", r.num_halfsphere_points);
    printf("     --ffbidx-num-angle-points=\n");
    printf("                            Number of sample points in rotation space for finding cell candidates (0: auto)\n");
    printf("                            Default: %d\n", r.num_angle_points);
    printf("\n");
    printf("     --ffbidx-min-spots=\n");
    printf("                            Refinement: Maximum number of indexed peaks to accept solution.\n");
    printf("                            Default: %d\n", i.min_spots);
    printf("     --ffbidx-ifssr-max-dist=\n");
    printf("                            Refinement: max distance to reciprocal spots for inliers\n");
    printf("                            Default: %f\n", i.max_distance);

    printf("     --ffbidx-ifssr-max-iter=\n");
    printf("                            Refinement: max number of iterations\n");
    printf("                            Default: %d\n", i.max_iter);

    printf("     --ffbidx-ifssr-threshold-contraction=\n");
    printf("                            Refinement: contract error threshold by this value in every iteration\n");
    printf("                            Default: %f\n", i.threshold_contraction);
}


int ffbidx_default_options(struct ffbidx_options **opts_ptr)
{
    struct ffbidx_options *opts;

    opts = malloc(sizeof(struct ffbidx_options));
    if ( opts == NULL ) return ENOMEM;

    set_defaults(&opts->cpers, &opts->cruntime, &opts->cifssr);
    opts->cpers.max_input_cells = 1;

    *opts_ptr = opts;
    return 0;
}


int sscanf_uint(const char *arg, unsigned *val) {
    int tmp = 0;
    int ret = sscanf(arg, "%d", &tmp);
    if (ret != 1)
        return -1;
    if  (tmp < 0)
        return -1;
    *val = tmp;
    return 1;
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
            if (sscanf_uint(arg, &(*opts_ptr)->cpers.max_spots) != 1) {
                ERROR("Invalid value for --ffbidx-max-spots\n");
                return EINVAL;
            }
            break;
        case 3 :
            if (sscanf_uint(arg,  &(*opts_ptr)->cifssr.min_spots) != 1) {
                ERROR("Invalid value for --ffbidx-min-spots\n");
                return EINVAL;
            }
            break;
        case 4 :
            if (sscanf_uint(arg,  &(*opts_ptr)->cpers.max_output_cells) != 1) {
                ERROR("Invalid value for --ffbidx-output-cells\n");
                return EINVAL;
            }
            break;
        case 5:
            if (sscanf_uint(arg,  &(*opts_ptr)->cpers.num_candidate_vectors) != 1) {
                ERROR("Invalid value for --ffbidx-num-candidate-vectors\n");
                return EINVAL;
            }
            break;
        case 6:
            (*opts_ptr)->cpers.redundant_computations = true;
            break;
        case 7:
            (*opts_ptr)->cpers.redundant_computations = false;
            break;
        case 8:
            if (sscanf(arg, "%f", &(*opts_ptr)->cifssr.max_distance) != 1) {
                ERROR("Invalid value for --ffbidx-ifssr-max-dist\n");
                return EINVAL;
            }
            break;
        case 9:
            if (sscanf_uint(arg,  &(*opts_ptr)->cifssr.max_iter) != 1) {
                ERROR("Invalid value for --ffbidx-ifssr-max-iter\n");
                return EINVAL;
            }
            break;
        case 10:
            if (sscanf(arg, "%f", &(*opts_ptr)->cifssr.threshold_contraction) != 1) {
                ERROR("Invalid value for --ffbidx-ifssr-threshold-contraction\n");
                return EINVAL;
            }
            break;
        case 11:
            if (sscanf(arg, "%f", &(*opts_ptr)->cruntime.length_threshold) != 1) {
                ERROR("Invalid value for --ffbidx-length-threshold\n");
                return EINVAL;
            }
            break;
        case 12:
            if (sscanf(arg, "%f", &(*opts_ptr)->cruntime.triml) != 1) {
                ERROR("Invalid value for --ffbidx-triml\n");
                return EINVAL;
            }
            break;
        case 13:
            if (sscanf(arg, "%f", &(*opts_ptr)->cruntime.trimh) != 1) {
                ERROR("Invalid value for --ffbidx-trimh\n");
                return EINVAL;
            }
            break;
        case 14:
            if (sscanf(arg, "%f", &(*opts_ptr)->cruntime.delta) != 1) {
                ERROR("Invalid value for --ffbidx-delta\n");
                return EINVAL;
            }
            break;
        case 15:
            if (sscanf(arg, "%f", &(*opts_ptr)->cruntime.dist1) != 1) {
                ERROR("Invalid value for --ffbidx-dist1\n");
                return EINVAL;
            }
            break;
        case 16:
            if (sscanf(arg, "%f", &(*opts_ptr)->cruntime.dist3) != 1) {
                ERROR("Invalid value for --ffbidx-dist3\n");
                return EINVAL;
            }
            break;

        case 17:
            if (sscanf_uint(arg,  &(*opts_ptr)->cruntime.num_halfsphere_points) != 1) {
                ERROR("Invalid value for --ffbidx-num-halfsphere-points\n");
                return EINVAL;
            }
            break;
        case 18:
            if (sscanf_uint(arg,  &(*opts_ptr)->cruntime.num_angle_points) != 1) {
                ERROR("Invalid value for --ffbidx-num-angle-points\n");
                return EINVAL;
            }
            break;
        default:
            break;
    }

    char msg[256];
    memset(msg, 0, 256);
    struct error ffbidx_err;
    ffbidx_err.msg_len = 255;
    ffbidx_err.message = msg;
    if (check_config(&(*opts_ptr)->cpers,
                     &(*opts_ptr)->cruntime, &
                     (*opts_ptr)->cifssr,
                     &ffbidx_err) != 0) {
        ERROR(msg);
        return EINVAL;
    }

    return 0;
}


static struct argp_option ffbidx_options[] = {
        {"help-ffbidx", 1, NULL, OPTION_NO_USAGE, "Show options for fast feedback indexing algorithm", 99},
        {"ffbidx-max-spots", 2, "ffbidx_maxn", OPTION_HIDDEN, NULL},
        {"ffbidx-min-spots", 3, "ffbidx_minn", OPTION_HIDDEN, NULL},
        {"ffbidx-output-cells", 4, "ffbidx_out_cells", OPTION_HIDDEN, NULL},
        {"ffbidx-num-candidate-vectors", 5, "ffbidx_num_candidate_vectors", OPTION_HIDDEN, NULL},
        {"ffbidx-redundant-computations", 6, NULL, OPTION_HIDDEN, NULL},
        {"ffbidx-no-redundant-computations", 7, NULL, OPTION_HIDDEN, NULL},
        {"ffbidx-ifssr-max-dist", 8, "ffbidx_max_dist", OPTION_HIDDEN, NULL},
        {"ffbidx-ifssr-max-iter", 9, "ffbidx_max_iter", OPTION_HIDDEN, NULL},
        {"ffbidx-ifssr-threshold-contraction", 10, "ffbidx_thrc", OPTION_HIDDEN, NULL},
        {"ffbidx-length-threshold", 11, "ffbidx_lt", OPTION_HIDDEN, NULL},
        {"ffbidx-triml", 12, "ffbidx_lt", OPTION_HIDDEN, NULL},
        {"ffbidx-trimh", 13, "ffbidx_lh", OPTION_HIDDEN, NULL},
        {"ffbidx-delta", 14, "ffbidx_delta", OPTION_HIDDEN, NULL},
        {"ffbidx-dist1", 15, "ffbidx_dist1", OPTION_HIDDEN, NULL},
        {"ffbidx-dist3", 16, "ffbidx_dist3", OPTION_HIDDEN, NULL},
        {"ffbidx-num-halfsphere-points", 17, "ffbidx_nhp", OPTION_HIDDEN, NULL},
        {"ffbidx-num-angle-points", 18, "ffbidx_nap", OPTION_HIDDEN, NULL},
        {0}
};


struct argp ffbidx_argp = { ffbidx_options, ffbidx_parse_arg,
                              NULL, NULL, NULL, NULL, NULL };
