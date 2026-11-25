/*
 * deltacchalf.c
 *
 * Check deltaCChalf calculation
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
#include <stdarg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <utils.h>
#include <symmetry.h>
#include <cell.h>
#include <cell-utils.h>
#include <fom.h>
#include <reflist.h>
#include <reflist-utils.h>

#include "../src/rejection.h"
#include "../src/merge.h"

static RefList *random_reflections(UnitCell *cell, double max_r, gsl_rng *rng,
                                   int *phmax, int *pkmax, int *plmax)
{
    double ax, ay, az;
    double bx, by, bz;
    double cx, cy, cz;
    double asx, asy, asz;
    double bsx, bsy, bsz;
    double csx, csy, csz;
    RefList *reflections;
    int hmax, kmax, lmax;
    signed int h, k, l;

    cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

    hmax = max_r * modulus(ax, ay, az);
    kmax = max_r * modulus(bx, by, bz);
    lmax = max_r * modulus(cx, cy, cz);
    *phmax = hmax; *pkmax = kmax; *plmax = lmax;

    if ( (hmax >= 512) || (kmax >= 512) || (lmax >= 512) ) {
        printf("Resolution is too high\n");
        return NULL;
    }

    reflections = reflist_new();

    cell_get_reciprocal(cell, &asx, &asy, &asz,
                              &bsx, &bsy, &bsz,
                              &csx, &csy, &csz);

    for ( h=-hmax; h<=hmax; h++ ) {
    for ( k=-kmax; k<=kmax; k++ ) {
    for ( l=-lmax; l<=lmax; l++ ) {

        Reflection *refl;

        if ( forbidden_reflection(cell, h, k, l) ) continue;
        if ( 2.0*resolution(cell, h, k, l) > max_r ) continue;

        refl = add_refl(reflections, h, k, l);
        set_intensity(refl, gsl_ran_gaussian_tail(rng, 0.01, 1000.0));
        set_esd_intensity(refl, 0.0);
        set_redundancy(refl, 1);
        set_partiality(refl, 1.0);
        set_lorentz(refl, 1.0);

    }
    }
    }

    return reflections;
}


static void add_crystal(struct crystal_refls *cr, int n_refl, RefList *all,
                        int hmax, int kmax, int lmax,
                        int junk, gsl_rng *rng, UnitCell *cell)
{
    int j;

    cr->cr = crystal_new();
    cr->refls = reflist_new();
    cr->image_owns_crystal = 1;
    cr->image_owns_refls = 1;

    crystal_set_cell(cr->cr, cell_new_from_cell(cell));

    for ( j=0; j<n_refl; j++ ) {

        Reflection *refl;
        Reflection *full;
        signed int h, k, l;
        double Ifull, sigI;

        do {
            h = gsl_rng_uniform_int(rng, hmax*2)-hmax;
            k = gsl_rng_uniform_int(rng, kmax*2)-kmax;
            l = gsl_rng_uniform_int(rng, lmax*2)-lmax;
            full = find_refl(all, h, k, l);
        } while ( full == NULL );

        refl = add_refl(cr->refls, h, k, l);

        Ifull = get_intensity(full);
        sigI = sqrt(Ifull);

        if ( !junk ) {
            set_intensity(refl, Ifull+gsl_ran_gaussian(rng, sigI));
        } else {
            set_intensity(refl, gsl_ran_gaussian(rng, sigI));
        }
        set_esd_intensity(refl, sigI);
        set_redundancy(refl, 1);
        set_partiality(refl, 1.0);
        set_lorentz(refl, 1.0);

    }
}


static struct crystal_refls *make_data(int n_cryst, UnitCell *cell, SymOpList *sym)
{
    int i;
    struct crystal_refls *cr;
    RefList *all;
    gsl_rng *rng;
    int hmax, kmax, lmax;

    rng = gsl_rng_alloc(gsl_rng_mt19937);

    FILE *fh;
    unsigned long int seed;

    fh = fopen("/dev/urandom", "r");
    if ( fh == NULL ) {
        ERROR("Failed to open /dev/urandom.  Results will not be random\n");
    }

    if ( fread(&seed, sizeof(seed), 1, fh) == 1 ) {
        gsl_rng_set(rng, seed);
    } else {
        ERROR("Failed to seed RNG\n");
    }
    fclose(fh);

    all = random_reflections(cell, 1/(2e-10), rng, &hmax, &kmax, &lmax);

    printf("Resolution range is %i %i %i\n", hmax, kmax, lmax);
    printf("%i reflections in complete data set\n", num_reflections(all));

    cr = malloc(n_cryst*sizeof(struct crystal_refls));
    if ( cr == NULL ) return NULL;

    for ( i=0; i<n_cryst-1; i++ ) {
        add_crystal(&cr[i], 100, all, hmax, kmax, lmax, 0, rng, cell);
    }

    add_crystal(&cr[n_cryst-1], 100, all, hmax, kmax, lmax, 1, rng, cell);

    return cr;
}


int main(int argc, char *argv[])
{
    int fail = 0;
    RefList *merged;
    int n_used = 0;
    ErrorModel *emodel;
    struct crystal_refls *cr;
    const int n_cryst = 1000;
    UnitCell *cell;
    SymOpList *sym;
    int i;
    int n_wrong = 0;

    cell = cell_new_from_parameters(50.0e-10, 60.0e-10, 110.0e-10,
                                    deg2rad(90.0), deg2rad(90.0), deg2rad(90.9));
    sym = get_pointgroup("mmm");
    cr = make_data(n_cryst, cell, sym);
    if ( cr == NULL ) return 1;

    emodel = error_model_new(EMODEL_EV06);
    merged = merge_intensities(cr, n_cryst, 1, 2, INFINITY, 1, 0, emodel, 1, &n_used);
    check_rejection(cr, n_cryst, merged, INFINITY, 0, 1);

    for ( i=0; i<n_cryst-1; i++ ) {
        if ( crystal_get_user_flag(cr[i].cr) ) {
            n_wrong++;
        }
    }

    if ( n_wrong ) {
        printf("%i crystals were wrongly marked as bad.\n", n_wrong);
        fail = 1;
    }

    if ( !crystal_get_user_flag(cr[n_cryst-1].cr) ) {
        printf("Crystal %i was wrongly marked as good.\n", n_cryst-1);
        fail = 1;
    }

    return fail;
}
