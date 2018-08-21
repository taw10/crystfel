/*
 * rejection.c
 *
 * Crystal rejection for scaling
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
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


static double calculate_cchalf(RefList *template, RefList *full,
                               Crystal *exclude)
{
	Reflection *trefl;
	RefListIterator *iter;
	double all_sum_mean = 0.0;
	double all_sumsq_mean = 0.0;
	double all_sum_var = 0.0;
	double sig2E, sig2Y;
	int n = 0;

	/* Iterate over all reflections */
	for ( trefl = first_refl(template, &iter);
	      trefl != NULL;
	      trefl = next_refl(trefl, iter) )
	{
		struct reflection_contributions *c;
		int j;
		double refl_sum = 0.0;
		double refl_sumsq = 0.0;
		double refl_mean, refl_var;
		signed int h, k, l;
		Reflection *refl;

		get_indices(trefl, &h, &k, &l);
		refl = find_refl(full, h, k, l);
		if ( refl == NULL ) continue;

		c = get_contributions(refl);
		assert(c != NULL);

		for ( j=0; j<c->n_contrib; j++ ) {

			double Ii;

			if ( c->contrib_crystals[j] == exclude ) {
				continue;
			}

			/* FIXME: Apply corrections */
			Ii = get_intensity(c->contribs[j]);

			refl_sum += Ii;
			refl_sumsq += Ii * Ii;

		}

		refl_mean = refl_sum / c->n_contrib;
		refl_var = refl_sumsq - refl_sum*refl_sum;

		all_sum_mean += refl_mean;
		all_sumsq_mean += refl_mean*refl_mean;
		all_sum_var += refl_var;
		n++;

	}

	sig2E = all_sum_var / n;
	sig2Y = all_sumsq_mean - all_sum_mean*all_sum_mean;

	return (sig2Y - 0.5*sig2E) / (sig2Y + 0.5*sig2E);
}


static void check_deltacchalf(Crystal **crystals, int n, RefList *full)
{
	double cchalf;
	int i;

	STATUS("Overall CChalf = %f\n", cchalf);

	for ( i=0; i<n; i++ ) {
		double cchalf, cchalfi;
		cchalf = calculate_cchalf(crystal_get_reflections(crystals[i]),
		                          full, NULL);
		cchalfi = calculate_cchalf(crystal_get_reflections(crystals[i]),
		                           full, crystals[i]);
		STATUS("DeltaCChalf_%i = %e\n", i, cchalfi - cchalf);
	}
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


void check_rejection(Crystal **crystals, int n, RefList *full, double max_B)
{
	int i;
	int n_acc = 0;

	/* Check according to delta CC½ */
	if ( full != NULL ) check_deltacchalf(crystals, n, full);

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
