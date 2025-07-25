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
 * Info:
 * - This code uses the c_api.h defined at
 *   https://github.com/paulscherrerinstitute/fast-feedback-indexer/blob/main/indexer/src/ffbidx/c_api.h
 * - struct ffbidx_options is defined in crystfels index.h
 */

#include <libcrystfel-config.h>
#include <stdio.h>
#include <stdlib.h>

#include "cell-utils.h"
#include "ffbidx.h"

#ifdef HAVE_FFBIDX

#include <ffbidx/c_api.h>

#define ERR_MSG_LEN 128u

/* ffbidx_private_data state */
#define STATE_UNDEF 0
#define STATE_CREATED 1
#define STATE_USED 2

struct ffbidx_private_data {
    UnitCell *cellTemplate;     // cell template from crystfel, set in ffbidx_prepare()
    struct ffbidx_options opts; // indexer persistent/runtime/refinement configuration
    struct input in;            // indexer input
    struct output out;          // indexer output
    int handle;                 // indexer state handle
    char msg[ERR_MSG_LEN];      // space for error message
    struct error err_msg;       // error message space description
    int state;
};

static void ffbidx_error_msg(struct error *err, const char* activity)
{
    const size_t msg_len = 128;
    char msg[msg_len];
   if ((err != NULL) && (err->message != NULL)) {
        snprintf(msg, msg_len, "ffbidx error - %s: %s\n", activity, err->message);
        ERROR(msg);
    } else {
        ERROR("ffbidx error - %s\n", activity);
    }
}

static void ffbidx_error(int handle, const char* activity)
{
    if (handle>=0) {
        ffbidx_error_msg(error_message(handle), activity);
    } else {
        ERROR("ffbidx: Invalid handle in ffbidx_error()!\n");
    }
}

static void check_ffbidx_state(const struct ffbidx_private_data *priv, int state_expected)
{
    if (priv == NULL) {
        fprintf(stderr, "ffbidx - FATAL: private data pointer is NULL!\n");
        exit(-1);
    }
    if (priv->state != state_expected) {
        fprintf(stderr, "ffbidx - FATAL: state is %d, but should be %d!\n", priv->state, state_expected);
        exit(-1);
    }
}

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
    unsigned npk, i;
    int bcell;

    struct ffbidx_private_data *priv = (struct ffbidx_private_data *) ipriv;
    check_ffbidx_state(priv, STATE_CREATED);
    priv->state = STATE_USED;

    // how many spots to consider
    npk = (unsigned)image_feature_count(image->features);
    if (npk < priv->opts.cifssr.min_spots) {
        priv->state = STATE_CREATED;
        return 0;   // too few spots
    }

    if (npk > priv->opts.cpers.max_spots)
        npk = priv->opts.cpers.max_spots;

    priv->in.n_spots = npk;

    // convert spots to indexer spot input format
    for ( i=0; i<npk; i++ ) {

        struct imagefeature *f;
        double r[3];

        f = image_get_feature(image->features, i);
        if ( f == NULL ) {
            ERROR("ffbidx: Empty feature ???");
            continue;
        }

        detgeom_transform_coords(&image->detgeom->panels[f->pn],
                                 f->fs, f->ss, image->lambda,
                                 0.0, 0.0, r);
        priv->in.spot.x[i] = r[0] * 1e-10;
        priv->in.spot.y[i] = r[1] * 1e-10;
        priv->in.spot.z[i] = r[2] * 1e-10;
    }

    // start indexing
    if (index_start(priv->handle, &priv->in, &priv->out, &priv->opts.cruntime, NULL, NULL)) {
        ffbidx_error(priv->handle, "index_start");
        priv->state = STATE_CREATED;
        return 0;
    }

    // wait for indexing output
    if (index_end(priv->handle, &priv->out)) {
        ffbidx_error(priv->handle, "index_end");
        priv->state = STATE_CREATED;
        return 0;
    }

    // refine output cells
    if (refine(priv->handle, &priv->in, &priv->out, &priv->opts.cifssr, 0, 1)) {
        ffbidx_error(priv->handle, "refine");
        priv->state = STATE_CREATED;
        return 0;
    }

    // cell with best score
    if ((bcell = best_cell(priv->handle, &priv->out)) < 0) {
        ffbidx_error(priv->handle, "best_cell");
        priv->state = STATE_CREATED;
        return 0;
    }

    if (priv->opts.cifssr.max_distance < priv->out.score[bcell]) {
        priv->state = STATE_CREATED;
        return 0;   // cell not viable
    }

    // convert best cell to crystfel output format
    UnitCell *uc;
    uc = cell_new();

    cell_set_cartesian(uc,
                       priv->out.x[bcell * 3 + 0] * 1e-10, priv->out.y[bcell * 3 + 0] * 1e-10, priv->out.z[bcell * 3 + 0] * 1e-10,
                       priv->out.x[bcell * 3 + 1] * 1e-10, priv->out.y[bcell * 3 + 1] * 1e-10, priv->out.z[bcell * 3 + 1] * 1e-10,
                       priv->out.x[bcell * 3 + 2] * 1e-10, priv->out.y[bcell * 3 + 2] * 1e-10, priv->out.z[bcell * 3 + 2] * 1e-10);

    makeRightHanded(uc);

    cell_set_lattice_type(uc, cell_get_lattice_type(priv->cellTemplate));
    cell_set_centering(uc, cell_get_centering(priv->cellTemplate));
    cell_set_unique_axis(uc, cell_get_unique_axis(priv->cellTemplate));

    priv->state = STATE_CREATED;

    // validate and produce result in format expected by crystfel
    if ( validate_cell(uc) ) {
        ERROR("ffbidx: problem with returned cell!\n");
        cell_free(uc);
        return 0;
    }

    Crystal *cr = crystal_new();
    if ( cr == NULL ) {
        ERROR("ffbidx: Failed to allocate crystal.\n");
        return 0;
    }
    crystal_set_cell(cr, uc);
    crystal_set_det_shift(cr, 0, 0);
    image_add_crystal(image, cr);
    return 1;
}

void *ffbidx_prepare(IndexingMethod indm, UnitCell *cell, struct ffbidx_options *opts) {
    if ( cell == NULL ) {
        ERROR("ffbidx: Unit cell information is required but missing.\n");
        return NULL;
    }

    struct ffbidx_private_data *priv = (struct ffbidx_private_data *) calloc(1, sizeof(struct ffbidx_private_data));

    priv->state = STATE_UNDEF;
    priv->cellTemplate = cell;
    priv->opts = *opts;
    priv->err_msg.message = priv->msg;
    priv->err_msg.msg_len = ERR_MSG_LEN;

    struct config_persistent *cp = &priv->opts.cpers; // for convenientce

    // allocate cell input space for one cell
    priv->in.cell.x = (float *) calloc(3, sizeof(float));
    priv->in.cell.y = (float *) calloc(3, sizeof(float));
    priv->in.cell.z = (float *) calloc(3, sizeof(float));
    priv->in.n_cells = 1;
    cp->max_input_cells = 1;    // make sure max is 1

    // allocate spot input space for max_spots
    priv->in.spot.x = (float *) calloc(cp->max_spots, sizeof(float));
    priv->in.spot.y = (float *) calloc(cp->max_spots, sizeof(float));
    priv->in.spot.z = (float *) calloc(cp->max_spots, sizeof(float));

    priv->in.new_cells = true; // TODO: only set to true on first indexer run
    priv->in.new_spots = true;

    // allocate cell and score output space for max_output_cells
    priv->out.x = (float *) calloc(3 * cp->max_output_cells, sizeof(float));
    priv->out.y = (float *) calloc(3 * cp->max_output_cells, sizeof(float));
    priv->out.z = (float *) calloc(3 * cp->max_output_cells, sizeof(float));
    priv->out.score = (float *) calloc(cp->max_output_cells, sizeof(float));
    priv->out.n_cells = cp->max_output_cells;

    // convert given cell to indexer cell input format
    double cell_internal_double[9];

    cell_get_cartesian(priv->cellTemplate,
                       &cell_internal_double[0],&cell_internal_double[3],&cell_internal_double[6],
                       &cell_internal_double[1],&cell_internal_double[4],&cell_internal_double[7],
                       &cell_internal_double[2],&cell_internal_double[5],&cell_internal_double[8]);

    for (int i = 0; i < 3; i++) {
        priv->in.cell.x[i] = cell_internal_double[0 + i] * 1e10;
        priv->in.cell.y[i] = cell_internal_double[3 + i] * 1e10;
        priv->in.cell.z[i] = cell_internal_double[6 + i] * 1e10;
    }

    // create indexer state
    priv->handle = create_indexer(&priv->opts.cpers, MEMORY_PIN_STATIC, &priv->err_msg, NULL);
    if (priv->handle < 0) {
        ffbidx_error_msg(&priv->err_msg, "create_indexer");
        return NULL;
    }

    priv->state = STATE_CREATED;
    return priv;
}

void ffbidx_cleanup(void *ipriv) {
    struct ffbidx_private_data *priv = (struct ffbidx_private_data *) ipriv;
    check_ffbidx_state(priv, STATE_CREATED);
    priv->state = STATE_UNDEF;

    // drop indexer state
    if (drop_indexer(priv->handle) != 0) {
        ERROR("ffbidx: Error dropping indexer\n");
    }

    // deallocate input spot space
    free(priv->in.spot.x);
    free(priv->in.spot.y);
    free(priv->in.spot.z);

    // deallocate input cell space
    free(priv->in.cell.x);
    free(priv->in.cell.y);
    free(priv->in.cell.z);

    // deallocate output cell and score space
    free(priv->out.x);
    free(priv->out.y);
    free(priv->out.z);
    free(priv->out.score);

    // deallocate private data
    free(priv);
}

const char *ffbidx_probe(UnitCell *cell) {
    return "ffbidx";
}

// fill in default option values
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

// parse unsigned integer
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

#else // HAVE_FFBIDX is not set

int run_ffbidx(struct image *image, void *ipriv)
{
	ERROR("This copy of CrystFEL was compiled without FFBIDX support.\n");
	return 0;
}


void *ffbidx_prepare(IndexingMethod indm, UnitCell *cell, struct ffbidx_options *opts)
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

int ffbidx_default_options(struct ffbidx_options **opts_ptr)
{
    return 0;
}

#endif

static void ffbidx_show_help(void)
{
#ifdef HAVE_FFBIDX  // print ffbidx option descriptions

    struct config_runtime r;
    struct config_persistent p;
    struct config_ifssr i;
    set_defaults(&p, &r, &i);
    printf("Parameters for the fast feedback indexing algorithm:\n");

    printf("     --ffbidx-max-spots=\n");
    printf("                            Maximum number of peaks used for indexing.\n");
    printf("                            Default: %d\n", p.max_spots);

    printf("     --ffbidx-output-cells=\n");
    printf("                            Number of output cells.\n");
    printf("                            Default: %d\n", p.max_output_cells);

    printf("     --ffbidx-num-candidate-vectors=\n");
    printf("                            Number of candidate vectors (per input cell vector).\n");
    printf("                            Default: %d\n", p.num_candidate_vectors);

    printf("     --ffbidx-redundant-computations\n");
    printf("     --ffbidx-no-redundant-computations\n");
    printf("                            Compute candidates for all three cell vectors instead of just one\n");
    printf("                            Default: %s\n", p.redundant_computations
                                                               ? "--ffbidx-redundant-computations"
                                                               : "--ffbidx-no-redundant-computations");

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

    printf("     --ffbidx-min-spots=\n");
    printf("                            Vector refinement: minimum number of peaks to fit against\n");
    printf("                            Default: %d\n", r.min_spots);

    printf("     --ffbidx-num-halfsphere-points=\n");
    printf("                            Number of sample points on half sphere for finding vector candidates\n");
    printf("                            Default: %d\n", r.num_halfsphere_points);

    printf("     --ffbidx-num-angle-points=\n");
    printf("                            Number of sample points in rotation space for finding cell candidates (0: auto)\n");
    printf("                            Default: %d\n", r.num_angle_points);

    printf("\n");
    printf("     --ffbidx-ifssr-min-spots=\n");
    printf("                            Cell refinement: minimum number of peaks to fit against\n");
    printf("                            Default: %d\n", i.min_spots);

    printf("     --ffbidx-ifssr-max-dist=\n");
    printf("                            Cell refinement: max distance to reciprocal spots for inliers\n");
    printf("                            Default: %f\n", i.max_distance);

    printf("     --ffbidx-ifssr-max-iter=\n");
    printf("                            Cell refinement: max number of iterations\n");
    printf("                            Default: %d\n", i.max_iter);

    printf("     --ffbidx-ifssr-threshold-contraction=\n");
    printf("                            Cell refinement: contract error threshold by this value in every iteration\n");
    printf("                            Default: %f\n", i.threshold_contraction);

#else  // HAVE_FFBIDX is not set

    printf("This copy of CrystFEL was compiled without FFBIDX support.\n");

#endif
}

static error_t ffbidx_parse_arg(int key, char *arg, struct argp_state *state)
{
#ifdef HAVE_FFBIDX  // parse all ffbidx options

    struct ffbidx_options **opts_ptr = state->input;
    int r;

    switch ( key ) {
        case ARGP_KEY_INIT :
            r = ffbidx_default_options(opts_ptr);
            if ( r ) return r;
            break;

        case 1:
            ffbidx_show_help();
            return EINVAL;

        case 2:
            if (sscanf_uint(arg, &(*opts_ptr)->cpers.max_spots) != 1) {
                ERROR("Invalid value for --ffbidx-max-spots\n");
                return EINVAL;
            }
            break;
        case 3:
            if (sscanf_uint(arg,  &(*opts_ptr)->cifssr.min_spots) != 1) {
                ERROR("Invalid value for --ffbidx-ifssr-min-spots\n");
                return EINVAL;
            }
            break;
        case 4:
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
        case 19:
            if (sscanf_uint(arg,  &(*opts_ptr)->cruntime.min_spots) != 1) {
                ERROR("Invalid value for --ffbidx-min-spots\n");
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
                     &(*opts_ptr)->cruntime,
                     &(*opts_ptr)->cifssr,
                     &ffbidx_err) != 0) {
        ERROR(msg);
        return EINVAL;
    }

#else // HAVE_FFBIDX is not set, parse only --help-ffbidx

    switch ( key ) {
        case 1:
            ffbidx_show_help();
            return EINVAL;

        default:
            break;
    }

#endif

    return 0;
}


static struct argp_option ffbidx_options[] = {
        {"help-ffbidx", 1, NULL, OPTION_NO_USAGE, "Show options for fast feedback indexing algorithm", 99},

#ifdef HAVE_FFBIDX  // these options are only relevant if ffbidx is present
        {"ffbidx-max-spots", 2, "ffbidx_maxn", OPTION_HIDDEN, NULL},
        {"ffbidx-ifssr-min-spots", 3, "ffbidx_cminn", OPTION_HIDDEN, NULL},
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
        {"ffbidx-min-spots", 19, "ffbidx_vminn", OPTION_HIDDEN, NULL},
#endif

        {0}
};


struct argp ffbidx_argp = { ffbidx_options, ffbidx_parse_arg,
                            NULL, NULL, NULL, NULL, NULL };
