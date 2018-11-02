/*
 * rejection.c
 *
 * Crystal rejection for scaling
 *
 * Copyright © 2012-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2018 Thomas White <taw@physics.org>
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
void early_rejection(Crystal **crystals, int n)
{
	int i;
	double m = 0.0;
	double mean_m;
	FILE *fh = fopen("reject.dat", "w");
	int n_flag = 0;

	for ( i=0; i<n; i++ ) {
		double u;
		RefList *list = crystal_get_reflections(crystals[i]);
		u = mean_intensity(list);
		m += u;
		fprintf(fh, "%i %f\n", i, u);
	}
	mean_m = m/n;
	for ( i=0; i<n; i++ ) {
		double u;
		RefList *list = crystal_get_reflections(crystals[i]);
		u = mean_intensity(list);
		if ( u/mean_m < 0.2 ) {
			crystal_set_user_flag(crystals[i], 5);
			n_flag++;
		}
	}
	fclose(fh);

	STATUS("Mean intensity/peak = %f ADU\n", m/n);
	STATUS("%i crystals flagged\n", n_flag);
}


static double calculate_refl_mean_var(RefList *full)
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
		assert(c != NULL);

		/* Calculate the resolution just once, using the cell from the
		 * first crystal to contribute, otherwise it takes too long */
		res = resolution(crystal_get_cell(c->contrib_crystals[0]),
		                 h, k, l);

		/* Mean of contributions */
		for ( j=0; j<c->n_contrib; j++ ) {

			double Ii, G, B;

			G = crystal_get_osf(c->contrib_crystals[j]);
			B = crystal_get_Bfac(c->contrib_crystals[j]);
			Ii = get_intensity(c->contribs[j]);

			/* Total (multiplicative) correction factor */
			Ii *= 1.0/G * exp(B*res*res) * get_lorentz(c->contribs[j])
			          / get_partiality(c->contribs[j]);

			Ex += Ii - K;
			Ex2 += (Ii - K) * (Ii - K);

		}

		if ( c->n_contrib < 2 ) continue;

		set_temp1(refl, Ex);
		set_temp2(refl, Ex2);
	}

}


static double calculate_cchalf(RefList *template, RefList *full,
                               Crystal *exclude, int *pnref)
{
	Reflection *trefl;
	RefListIterator *iter;
	double sig2E, sig2Y;
	int n = 0;
	double wSum = 0.0;
	double wSum2 = 0.0;
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

		/* Remove contribution(s) from the excluded crystal */
		if ( exclude != NULL ) {
			refl = find_refl(crystal_get_reflections(exclude), h, k, l);
		} else {
			refl = NULL;
		}

		while ( refl != NULL ) {

			double G, B;
			double Ii = get_intensity(refl);

			G = crystal_get_osf(exclude);
			B = crystal_get_Bfac(exclude);

			/* Total (multiplicative) correction factor */
			Ii *= 1.0/G * exp(B*res*res) * get_lorentz(refl)
			          / get_partiality(refl);

			/* Remove contribution of this reflection */
			Ex -= Ii - K;
			Ex2 -= (Ii - K)*(Ii - K);

			n_removed++;

			refl = next_found_refl(refl);
		}

		if ( c->n_contrib - n_removed < 2 ) continue;

		refl_mean = K + (Ex / (c->n_contrib - n_removed));
		refl_var = (Ex2 - (Ex*Ex)/(c->n_contrib - n_removed)) / (c->n_contrib - n_removed - 1);
		refl_var /= (c->n_contrib - n_removed) / 2.0;

		all_sum_var += refl_var;
		n++;

		/* Running variance calculation to get sig2Y */
		wSum += w;
		wSum2 += w*w;
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


static void check_deltacchalf(Crystal **crystals, int n, RefList *full)
{
	double cchalf;
	int i;
	double *vals;
	double mean, sd;
	int nref = 0;

	calculate_refl_mean_var(full);

	cchalf = calculate_cchalf(full, full, NULL, &nref);
	STATUS("Overall CChalf = %f %% (%i reflections)\n", cchalf*100.0, nref);

	vals = malloc(n*sizeof(double));
	if ( vals == NULL ) {
		ERROR("Not enough memory for deltaCChalf check\n");
		return;
	}

	for ( i=0; i<n; i++ ) {
		double cchalf, cchalfi;
		RefList *template = crystal_get_reflections(crystals[i]);
		cchalf = calculate_cchalf(template, full, NULL, &nref);
		cchalfi = calculate_cchalf(template, full, crystals[i], &nref);
		//STATUS("Frame %i:", i);
		//STATUS("   With = %f  ", cchalf*100.0);
		//STATUS("Without = %f", cchalfi*100.0);
		//STATUS("  Delta = %f  ", (cchalf - cchalfi)*100.0);
		//STATUS("(nref = %i)\n", nref);
		vals[i] = cchalf - cchalfi;
		progress_bar(i, n-1, "Calculating deltaCChalf");
	}

	mean = gsl_stats_mean(vals, 1, n);
	sd = gsl_stats_sd_m(vals, 1, n, mean);
	STATUS("deltaCChalf = %f ± %f %%\n", mean*100.0, sd*100.0);

	for ( i=0; i<n; i++ ) {
		if ( vals[i] < mean-2.0*sd ) {
			crystal_set_user_flag(crystals[i], PRFLAG_DELTACCHALF);
		}
	}

	free(vals);
}


static void show_duds(Crystal **crystals, int n_crystals)
{
	int j;
	int bads[32];
	int any_bad = 0;

	for ( j=0; j<32; j++ ) bads[j] = 0;

	for ( j=0; j<n_crystals; j++ ) {
		int flag;
		flag = crystal_get_user_flag(crystals[j]);
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


void check_rejection(Crystal **crystals, int n, RefList *full, double max_B,
                     int no_deltacchalf)
{
	int i;
	int n_acc = 0;

	/* Check according to delta CC½ */
	if ( !no_deltacchalf && (full != NULL) ) {
		 check_deltacchalf(crystals, n, full);
	}

	for ( i=0; i<n; i++ ) {
		if ( fabs(crystal_get_Bfac(crystals[i])) > max_B ) {
			crystal_set_user_flag(crystals[i], PRFLAG_BIGB);
		}
		if ( crystal_get_user_flag(crystals[i]) == 0 ) n_acc++;
	}

	show_duds(crystals, n);

	if ( n_acc < 2 ) {
		ERROR("Not enough crystals left to proceed (%i).  Sorry.\n",
		      n_acc);
		exit(1);
	}
}
