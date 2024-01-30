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


static int calculate_refl_mean_var(RefList *full)
{
	Reflection *refl;
	RefListIterator *iter;
	signed int oh = 0;
	signed int ok = 0;
	signed int ol = 0;

	/* Iterate over all reflections */
	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		struct reflection_contributions *c;
		int j;
		signed int h, k, l;
		double res;
		double K;
		double Ex = 0.0;
		double Ex2 = 0.0;

		get_indices(refl, &h, &k, &l);

		/* next_refl() will iterate over multiple copies of the same
		 * unique reflection, but we are only interested in seeing each
		 * unique reflection once. */
		if ( (h==oh) && (k==ok) && (l==ol) ) continue;
		oh = h;  ok = k;  ol = l;

		/* We use the mean (merged) intensity as the reference point
		 * for shifting the data in the variance calculation */
		K = get_intensity(refl);

		c = get_contributions(refl);
		if ( c == NULL ) return 1;

		/* Calculate the resolution just once, using the cell from the
		 * first crystal to contribute, otherwise it takes too long */
		res = resolution(crystal_get_cell(c->contrib_crystals[0]),
		                 h, k, l);

		/* Mean of contributions */
		for ( j=0; j<c->n_contrib; j++ ) {

			double Ii, G, B;

			G = crystal_get_osf(c->contrib_crystals[j]);
			B = crystal_get_Bfac(c->contrib_crystals[j]);
			Ii = correct_reflection(get_intensity(c->contribs[j]), c->contribs[j], G, B, res);

			Ex += Ii - K;
			Ex2 += (Ii - K) * (Ii - K);

		}

		if ( c->n_contrib < 2 ) continue;

		set_temp1(refl, Ex);
		set_temp2(refl, Ex2);
	}

	return 0;
}


static double calculate_cchalf(RefList *template, RefList *full,
                               Crystal *exclude, int *pnref)
{
	Reflection *trefl;
	RefListIterator *iter;
	double sig2E, sig2Y;
	int n = 0;
	double wSum = 0.0;
	double mean = 0.0;
	double S = 0.0;
	double all_sum_var = 0.0;
	signed int oh = 0;
	signed int ok = 0;
	signed int ol = 0;

	/* "template" is the list of reflections to be included in CChalf */
	for ( trefl = first_refl(template, &iter);
	      trefl != NULL;
	      trefl = next_refl(trefl, iter) )
	{
		signed int h, k, l;
		double refl_mean, refl_var;
		double Ex, Ex2, K;
		int n_removed = 0;
		double w = 1.0;
		double meanOld;
		double res;
		struct reflection_contributions *c;
		Reflection *refl;
		Reflection *exrefl;

		/* The values we need are stored in the "full" list, not the
		 * template list */
		get_indices(trefl, &h, &k, &l);

		if ( (h==oh) && (k==ok) && (l==ol) ) continue;
		oh = h;  ok = k;  ol = l;

		refl = find_refl(full, h, k, l);
		if ( refl == NULL ) continue;  /* However, there might not have
		                                * been enough measurements for
		                                * it to appear in "full" */

		/* We use the mean (merged) intensity as the reference point
		 * for shifting the data in the variance calculation */
		K = get_intensity(refl);
		Ex = get_temp1(refl);
		Ex2 = get_temp2(refl);
		c = get_contributions(refl);
		assert(c != NULL);

		/* Resolution from first contributing crystal, like above */
		res = resolution(crystal_get_cell(c->contrib_crystals[0]),
		                 h, k, l);

		/* Remove contribution(s) from the excluded crystal.
		 * If the crystal is marked as bad, we should not remove it
		 * because it did not contribute in the first place.  */
		if ( exclude != NULL && !crystal_get_user_flag(exclude) ) {
			exrefl = find_refl(template, h, k, l);
		} else {
			exrefl = NULL;
		}

		while ( exrefl != NULL ) {

			double G, B;

			G = crystal_get_osf(exclude);
			B = crystal_get_Bfac(exclude);

			if ( get_partiality(exrefl) > MIN_PART_MERGE ) {

				double Ii = correct_reflection(get_intensity(exrefl), exrefl, G, B, res);

				/* Remove contribution of this reflection */
				Ex -= Ii - K;
				Ex2 -= (Ii - K)*(Ii - K);

				n_removed++;

			}

			exrefl = next_found_refl(exrefl);
		}

		if ( c->n_contrib - n_removed < 2 ) continue;

		refl_mean = K + (Ex / (c->n_contrib - n_removed));
		refl_var = (Ex2 - (Ex*Ex)/(c->n_contrib - n_removed)) / (c->n_contrib - n_removed - 1);
		refl_var /= (c->n_contrib - n_removed) / 2.0;

		all_sum_var += refl_var;
		n++;

		/* Running variance calculation to get sig2Y */
		wSum += w;
		meanOld = mean;
		mean = meanOld + (w/wSum) * (refl_mean - meanOld);
		S += w * (refl_mean - meanOld) * (refl_mean - mean);
	}

	sig2E = all_sum_var / n;
	sig2Y = S / (wSum - 1.0);

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
	Crystal *crystal;
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
	wargs->crystal = qargs->crystals[qargs->n_started].cr;
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
	int nref = 0;
	cchalf = calculate_cchalf(wargs->refls, wargs->full, NULL, &nref);
	cchalfi = calculate_cchalf(wargs->refls, wargs->full, wargs->crystal, &nref);
	//STATUS("Frame %i:", i);
	//STATUS("   With = %f  ", cchalf*100.0);
	//STATUS("Without = %f", cchalfi*100.0);
	//STATUS("  Delta = %f  ", (cchalf - cchalfi)*100.0);
	//STATUS("(nref = %i)\n", nref);
	if ( nref == 0 ) {
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

	if ( calculate_refl_mean_var(full) ) {
		STATUS("No reflection contributions for deltaCChalf "
		       "calculation (using reference reflections?)\n");
		return;
	}

	cchalf = calculate_cchalf(full, full, NULL, &nref);
	STATUS("Overall CChalf = %f %% (%i reflections)\n", cchalf*100.0, nref);

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
