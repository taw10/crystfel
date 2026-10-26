/*
 * indextest.c
 *
 * Unit test for indexing
 *
 * Copyright © 2026 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2026 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <image.h>
#include <index.h>
#include <cell-utils.h>

#include "indextest.h"

int main(int argc, char *argv[])
{
    DataTemplate *dtempl;
    struct image *image;
    int i;
    IndexingPrivate *ipriv;
    UnitCell *cell;
    double wl, clen;
    int flags;
    struct taketwo_options *taketwoopts;
    struct felix_options *felix_opts;
    struct xgandalf_options *xgandalf_opts;
    struct pinkindexer_options *pinkIndexer_opts;
    struct fromfile_options *fromfile_opts;
    struct smallcell_options *smallcell_opts;
    struct asdf_options *asdf_opts;
    struct ffbidx_options *ffbidx_opts;
    char *algo;
    int n_peaks;
    float *peak_fs;
    float *peak_ss;
    int *peak_pn;
    int refine = 0;
    int multi = 0;
    int n_exp = 0;
    float tols[] = {0.1, 0.1, 0.1, deg2rad(5.0), deg2rad(5.0), deg2rad(5.0)};

    if ( argc < 4 ) {
        fprintf(stderr, "Need at least four arguments\n");
        return 1;
    }

    dtempl = data_template_new_from_file(argv[1]);
    if ( dtempl == NULL ) {
        fprintf(stderr, "Failed to load geometry file\n");
        return 1;
    }

    algo = argv[2];
    if ( strcmp(argv[3], "multi") == 0 ) {
        multi = 1;
    }
    if ( strcmp(argv[4], "refine") == 0 ) {
        refine = 1;
    }
    printf("Testing %s (%s %s), refine=%i, multi=%i\n", algo, argv[3], argv[4], refine, multi);

    image = image_create_for_simulation(dtempl);
    if ( image == NULL ) {
        fprintf(stderr, "Failed to create image\n");
        return 1;
    }

    /* Select input data depending on whether we are testing
     * multi-lattice indexing or not */
    if ( strcmp(algo, "smallcell") == 0 ) {
        cell = cell_new_from_parameters(12.3e-10, 4.5e-10, 8.0e-10,
                                        deg2rad(120), deg2rad(97), deg2rad(90));
        cell_set_lattice_type(cell, L_TRICLINIC);
        cell_set_centering(cell, 'P');
        if ( multi ) {
            peak_fs = t4_peak_fs;
            peak_ss = t4_peak_ss;
            peak_pn = t4_peak_pn;
            n_peaks = t4_n_peaks;
            n_exp = 2;
        } else {
            peak_fs = t3_peak_fs;
            peak_ss = t3_peak_ss;
            peak_pn = t3_peak_pn;
            n_peaks = t3_n_peaks;
            n_exp = 1;
        }
    } else {
        cell = cell_new_from_parameters(123e-10, 45e-10, 80e-10,
                                        deg2rad(90), deg2rad(97), deg2rad(90));
        cell_set_lattice_type(cell, L_MONOCLINIC);
        cell_set_unique_axis(cell, 'b');
        cell_set_centering(cell, 'P');
        if ( multi ) {
            peak_fs = t2_peak_fs;
            peak_ss = t2_peak_ss;
            peak_pn = t2_peak_pn;
            n_peaks = t2_n_peaks;
            n_exp = 2;
        } else {
            peak_fs = t1_peak_fs;
            peak_ss = t1_peak_ss;
            peak_pn = t1_peak_pn;
            n_peaks = t1_n_peaks;
            n_exp = 1;
        }
    }

    image->features = image_feature_list_new();
    for ( i=0; i<n_peaks; i++ ) {
        image_add_feature(image->features, peak_fs[i], peak_ss[i], peak_pn[i], 1.0, NULL);
    }

    wl = data_template_get_wavelength_if_possible(dtempl);
    clen = data_template_get_clen_if_possible(dtempl);

    flags = INDEXING_RETRY | INDEXING_CHECK_PEAKS | INDEXING_CHECK_CELL;
    if ( refine ) flags |= INDEXING_REFINE;
    if ( multi ) flags |= INDEXING_MULTI;
    default_method_options(&taketwoopts,
                           &xgandalf_opts,
                           &ffbidx_opts,
                           &pinkIndexer_opts,
                           &felix_opts,
                           &fromfile_opts,
                           &smallcell_opts,
                           &asdf_opts);

    /* Set XGandalf fast mode */
    xgandalf_opts->sampling_pitch = 2;
    xgandalf_opts->grad_desc_iterations = 3;

    ipriv = setup_indexing(algo, cell, tols, flags, wl, clen, 1,
                           taketwoopts, xgandalf_opts,
                           ffbidx_opts, pinkIndexer_opts,
                           felix_opts, NULL, smallcell_opts, asdf_opts);

    if ( ipriv == NULL ) {
        fprintf(stderr, "Failed to set up indexing\n");
        return 1;
    }

    index_pattern(image, ipriv);

    printf("Found %i crystals\n", image->n_crystals);
    for ( i=0; i<image->n_crystals; i++ ) {
        cell_print_full(crystal_get_cell(image->crystals[i].cr));
    }

    if ( image->n_crystals != n_exp ) {
        fprintf(stderr, "Got %i crystals (should be %i)\n", image->n_crystals, n_exp);
        return 1;
    }

    /* Junk left behind by Felix */
    unlink("gmon.out");

    return 0;
}
