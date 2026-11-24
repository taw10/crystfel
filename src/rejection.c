/*
 * rejection.c
 *
 * Crystal rejection for scaling
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2019 Thomas White <taw@physics.org>
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


#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_statistics.h>

#include "crystal.h"
#include "reflist.h"
#include "rejection.h"
#include "cell-utils.h"
#include "post-refinement.h"
#include "merge.h"


static double mean_intensity(RefList *list)
{
    Reflection *refl;
    RefListIterator *iter;
    double total = 0.0;
    int n = 0;

    for ( refl = first_refl(list, &iter);
          refl != NULL;
          refl = next_refl(refl, iter) )
    {
        total += get_intensity(refl);
        n++;
    }

    return total/n;
}


/* Reject really obvious outliers */
void early_rejection(struct crystal_refls *crystals, int n)
{
    int i;
    double m = 0.0;
    double mean_m;
    FILE *fh = fopen("reject.dat", "w");
    int n_flag = 0;

    if ( fh == NULL ) {
        ERROR("Failed to open reject.dat\n");
        return;
    }

    for ( i=0; i<n; i++ ) {
        double u;
        u = mean_intensity(crystals[i].refls);
        m += u;
        fprintf(fh, "%i %f\n", i, u);
    }
    mean_m = m/n;
    for ( i=0; i<n; i++ ) {
        double u;
        u = mean_intensity(crystals[i].refls);
        if ( u/mean_m < 0.2 ) {
            crystal_set_user_flag(crystals[i].cr, 5);
            n_flag++;
        }
    }
    fclose(fh);

    STATUS("Mean intensity/peak = %f ADU\n", m/n);
    STATUS("%i crystals flagged\n", n_flag);
}


static double calculate_cchalf(RefList *template, RefList *full,
                               int crystal_number, int do_exclude, int *pnref)
{
    Reflection *trefl;
    RefListIterator *iter;
    double sig2E, sig2Y;
    int n = 0;
    signed int oh = 0;
    signed int ok = 0;
    signed int ol = 0;

    double *means = NULL;
    double *variances = NULL;
    int max_n = 0;

    double *Ii = NULL;
    double *weights = NULL;
    int max_I = 0;

    /* "template" is the list of reflections to be included in CChalf */
    for ( trefl = first_refl(template, &iter);
          trefl != NULL;
          trefl = next_refl(trefl, iter) )
    {
        signed int h, k, l;
        struct reflection_contributions *c;
        Reflection *refl;
        int ninc = 0;
        int i;
        int n_would_inc = 0;

        /* The values we need are stored in the "full" list, not the
         * template list */
        get_indices(trefl, &h, &k, &l);

        if ( (h==oh) && (k==ok) && (l==ol) ) continue;
        oh = h;  ok = k;  ol = l;

        refl = find_refl(full, h, k, l);
        if ( refl == NULL ) continue;  /* There might not have
                                        * been enough measurements for
                                        * it to appear in "full" */

        if ( get_intensity(refl) < 3.0*get_esd_intensity(refl) ) continue;

        c = get_contributions(refl);
        assert(c != NULL);

        if ( c->n_contrib > max_I ) {
            Ii = realloc(Ii, c->n_contrib*sizeof(double));
            weights = realloc(weights, c->n_contrib*sizeof(double));
            max_I = c->n_contrib;
            if ( (Ii == NULL) || (weights == NULL) ) {
                ERROR("Failed to allocate memory for CChalf calculation\n");
                return NAN;
            }
        }

        if ( n >= max_n ) {
            max_n += 1024;
            means = realloc(means, max_n*sizeof(double));
            variances = realloc(variances, max_n*sizeof(double));
            if ( (means == NULL) || (variances == NULL) ) {
                ERROR("Failed to allocate memory for CChalf calculation\n");
                return NAN;
            }
        }

        for ( i=0; i<c->n_contrib; i++ ) {

            int is_crystal = (c->contrib_crystal_idxs[i] == crystal_number);

            if ( !do_exclude && !is_crystal ) n_would_inc++;
            if ( do_exclude && is_crystal ) continue;

            Ii[ninc] = c->contribs[i];
            weights[ninc] = c->contrib_legacy_weights[i];
            ninc++;

        }

        if ( ninc < 2 ) continue;

        /* Even if we are not excluding a crystal, we need to make sure that
         * the calculation uses the same reflections that would be considered
         * if we were excluding. */
        if ( !do_exclude && n_would_inc < 2 ) continue;

        means[n] = gsl_stats_wmean(weights, 1, Ii, 1, ninc);
        variances[n] = gsl_stats_wvariance_m(weights, 1, Ii, 1, ninc, means[n])/(ninc/2.0);
        n++;

    }

    sig2E = gsl_stats_mean(variances, 1, n);
    sig2Y = gsl_stats_variance(means, 1, n);

    free(Ii);
    free(weights);

    free(variances);
    free(means);

    if ( pnref != NULL ) {
        *pnref = n;
    }
    return (sig2Y - 0.5*sig2E) / (sig2Y + 0.5*sig2E);
}


struct deltacchalf_queue_args
{
    RefList *full;
    struct crystal_refls *crystals;
    int n_crystals;
    int n_done;
    int n_started;
    int n_non;
    int n_nan;
    double *vals;
};


struct deltacchalf_worker_args
{
    RefList *full;
    RefList *refls;
    int crystal_number;
    int non;
    int nan;
    double deltaCChalf;
};


static void *create_deltacchalf_job(void *vqargs)
{
    struct deltacchalf_worker_args *wargs;
    struct deltacchalf_queue_args *qargs = vqargs;

    wargs = malloc(sizeof(struct deltacchalf_worker_args));

    wargs->full = qargs->full;
    wargs->refls = qargs->crystals[qargs->n_started].refls;
    wargs->crystal_number = qargs->n_started;
    wargs->non = 0;
    wargs->nan = 0;

    qargs->n_started++;

    return wargs;
}


static void run_deltacchalf_job(void *vwargs, int cookie)
{
    double cchalf, cchalfi;
    struct deltacchalf_worker_args *wargs = vwargs;
    int nref1 = 0;
    int nref2 = 0;
    cchalf = calculate_cchalf(wargs->refls, wargs->full, wargs->crystal_number, 0, &nref1);
    cchalfi = calculate_cchalf(wargs->refls, wargs->full, wargs->crystal_number, 1, &nref2);
    assert(nref1 == nref2);
    //STATUS("Frame %i:", wargs->crystal_number);
    //STATUS("   With = %f  ", cchalf*100.0);
    //STATUS("Without = %f", cchalfi*100.0);
    //STATUS("  Delta = %f  ", (cchalf - cchalfi)*100.0);
    //STATUS("(nref = %i)\n", nref1);
    if ( nref1 == 0 ) {
        wargs->deltaCChalf = 0.0;
        wargs->non = 1;
    } else {
        wargs->deltaCChalf = cchalf - cchalfi;
        if ( isnan(wargs->deltaCChalf) || isinf(wargs->deltaCChalf) ) {
            wargs->deltaCChalf = 0.0;
            wargs->nan = 1;
        }
    }
}


static void finalise_deltacchalf_job(void *vqargs, void *vwargs)
{
    struct deltacchalf_queue_args *qargs = vqargs;
    struct deltacchalf_worker_args *wargs = vwargs;
    qargs->n_done++;
    if ( wargs->nan ) qargs->n_nan++;
    if ( wargs->non ) qargs->n_non++;
    qargs->vals[wargs->crystal_number] = wargs->deltaCChalf;
    progress_bar(qargs->n_done, qargs->n_crystals,
                 "Calculating deltaCChalf");
    free(vwargs);
}


static void check_deltacchalf(struct crystal_refls *crystals, int n,
                              RefList *full, int n_threads)
{
    double cchalf;
    int i;
    double *vals;
    double mean, sd;
    int nref = 0;
    struct deltacchalf_queue_args qargs;

    cchalf = calculate_cchalf(full, full, -1, 0, &nref);
    STATUS("Overall CChalf = %f %% (%i reflections with I/sigI>3)\n", cchalf*100.0, nref);

    vals = malloc(n*sizeof(double));
    if ( vals == NULL ) {
        ERROR("Not enough memory for deltaCChalf check\n");
        return;
    }

    qargs.full = full;
    qargs.crystals = crystals;
    qargs.n_started = 0;
    qargs.n_crystals = n;
    qargs.n_done = 0;
    qargs.n_nan = 0;
    qargs.n_non = 0;
    qargs.vals = vals;
    run_threads(n_threads, run_deltacchalf_job, create_deltacchalf_job,
                finalise_deltacchalf_job, &qargs, n, 0, 0, 0);

    if ( qargs.n_non > 0 ) {
        STATUS("WARNING: %i patterns had no reflections in deltaCChalf "
               "calculation (I set deltaCChalf=zero for them)\n",
               qargs.n_non);
    }
    if ( qargs.n_nan > 0 ) {
        STATUS("WARNING: %i NaN or inf deltaCChalf values were "
               "replaced with zero\n", qargs.n_nan);
    }

    mean = gsl_stats_mean(vals, 1, n);
    sd = gsl_stats_sd_m(vals, 1, n, mean);
    STATUS("deltaCChalf = %f ± %f %%\n", mean*100.0, sd*100.0);

    for ( i=0; i<n; i++ ) {
        if ( isnan(vals[i]) || isinf(vals[i])
          || ((vals[i]<0.0) && (vals[i] < mean-2.0*sd)) )
        {
            crystal_set_user_flag(crystals[i].cr, PRFLAG_DELTACCHALF);
        }
    }

    free(vals);
}


static void show_duds(struct crystal_refls *crystals, int n_crystals)
{
    int j;
    int bads[32];
    int any_bad = 0;

    for ( j=0; j<32; j++ ) bads[j] = 0;

    for ( j=0; j<n_crystals; j++ ) {
        int flag;
        flag = crystal_get_user_flag(crystals[j].cr);
        assert(flag < 32);
        bads[flag]++;
        if ( flag != PRFLAG_OK ) any_bad++;
    }

    if ( any_bad ) {
        STATUS("%i bad crystals:\n", any_bad);
        for ( j=0; j<32; j++ ) {
            if ( bads[j] ) {
                STATUS("  %i %s\n", bads[j], str_prflag(j));
            }
        }
    }
}


void check_rejection(struct crystal_refls *crystals, int n, RefList *full,
                     double max_B, int no_deltacchalf, int n_threads)
{
    int i;
    int n_acc = 0;

    /* Check according to delta CC½ */
    if ( !no_deltacchalf && (full != NULL) ) {
        check_deltacchalf(crystals, n, full, n_threads);
    }

    for ( i=0; i<n; i++ ) {
        if ( fabs(crystal_get_Bfac(crystals[i].cr)) > max_B ) {
            crystal_set_user_flag(crystals[i].cr, PRFLAG_BIGB);
        }
        if ( crystal_get_user_flag(crystals[i].cr) == 0 ) n_acc++;
    }

    show_duds(crystals, n);

    if ( n_acc < 2 ) {
        ERROR("Not enough crystals left to proceed (%i).  Sorry.\n",
              n_acc);
        exit(1);
    }
}
