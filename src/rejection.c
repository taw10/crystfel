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


struct contribs
{
	int          serial;  /* h,k,l */
	int          n_contrib;
	int          max_contrib;
	Reflection **contribs;
	Crystal    **contrib_crystals;
};


struct contributionlist
{
	struct contribs *contribs;
	int n;  /* Number of reflections */
};


static int alloc_contribs(struct contribs *c)
{
	c->contribs = realloc(c->contribs, c->max_contrib*sizeof(Reflection *));
	c->contrib_crystals = realloc(c->contrib_crystals,
	                              c->max_contrib*sizeof(Crystal *));
	if ( c->contribs == NULL ) return 1;
	if ( c->contrib_crystals == NULL ) return 1;
	return 0;
}


static void free_contributions(struct contributionlist *clist)
{
	int i;
	for ( i=0; i<clist->n; i++ ) {
		free(clist->contribs[i].contribs);
		free(clist->contribs[i].contrib_crystals);
	}
	free(clist->contribs);
	free(clist);
}


static struct contributionlist *find_all_contributions(Crystal **crystals,
                                                       int n, RefList *full)
{
	Reflection *refl;
	RefListIterator *iter;
	struct contributionlist *clist;
	int nc = 0;

	clist = malloc(sizeof(struct contributionlist));
	if ( clist == NULL ) return NULL;

	clist->n = num_reflections(full);
	clist->contribs = malloc(clist->n*sizeof(struct contribs));
	if ( clist->contribs == NULL ) {
		free(clist);
		return NULL;
	}

	for ( refl = first_refl(full, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		int i;
		signed int h, k, l;
		struct contribs *c;

		get_indices(refl, &h, &k, &l);
		c = &clist->contribs[nc++];

		c->serial = SERIAL(h, k, l);
		c->n_contrib = 0;
		c->max_contrib = 32;
		c->contribs = NULL;
		c->contrib_crystals = NULL;
		if ( alloc_contribs(c) ) return NULL;

		/* Find all versions of this reflection */
		for ( i=0; i<n; i++ ) {

			Reflection *refl2;
			RefList *list2 = crystal_get_reflections(crystals[i]);
			refl2 = find_refl(list2, h, k, l);

			if ( refl2 == NULL ) continue;
			do {
				c->contribs[c->n_contrib] = refl2;
				c->contrib_crystals[c->n_contrib] = crystals[i];
				c->n_contrib++;

				if ( c->n_contrib == c->max_contrib ) {
					c->max_contrib += 64;
					if ( alloc_contribs(c) ) return NULL;
				}

				refl2 = next_found_refl(refl2);

			} while ( refl2 != NULL );
		}

		progress_bar(nc, clist->n, "Collecting contributions");
	}

	return clist;
}


static struct contribs *lookup_contribs(struct contributionlist *clist,
                                        signed int h, signed int k, signed int l)
{
	int i, serial;

	serial = SERIAL(h, k, l);

	for ( i=0; i<clist->n; i++ ) {
		if ( clist->contribs[i].serial == serial ) {
			return &clist->contribs[i];
		}
	}
	return NULL;
}


static double calculate_cchalf(struct contributionlist *clist, Crystal *exclude)
{
	int i;
	double all_sum_mean = 0.0;
	double all_sumsq_mean = 0.0;
	double all_sum_var = 0.0;
	double sig2E, sig2Y;

	/* Iterate over all reflections */
	for ( i=0; i<clist->n; i++ ) {

		struct contribs *c;
		int j;
		double refl_sum = 0.0;
		double refl_sumsq = 0.0;
		double refl_mean, refl_var;

		c = &clist->contribs[i];

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

	}

	sig2E = all_sum_var / clist->n;
	sig2Y = all_sumsq_mean - all_sum_mean*all_sum_mean;

	return (sig2Y - 0.5*sig2E) / (sig2Y + 0.5*sig2E);
}


static void check_deltacchalf(Crystal **crystals, int n, RefList *full)
{
	struct contributionlist *clist;
	double cchalf;
	int i;

	clist = find_all_contributions(crystals, n, full);

	cchalf = calculate_cchalf(clist, NULL);
	STATUS("Overall CChalf = %f\n", cchalf);

	for ( i=0; i<n; i++ ) {
		double cchalfi;
		cchalfi = calculate_cchalf(clist, crystals[i]);
		STATUS("DeltaCChalf_%i = %e\n", i, cchalfi - cchalf);
	}

	free_contributions(clist);
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
