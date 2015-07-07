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


static void check_cc(Crystal *cr, RefList *full)
{
	RefList *list = crystal_get_reflections(cr);
	Reflection *refl;
	RefListIterator *iter;
	double G = crystal_get_osf(cr);
	double B = crystal_get_Bfac(cr);
	UnitCell *cell = crystal_get_cell(cr);
	double *vec1, *vec2;
	int n, i;
	double cc;

	n = num_reflections(list);
	vec1 = malloc(n*sizeof(double));
	vec2 = malloc(n*sizeof(double));
	if ( (vec1 == NULL) || (vec2 == NULL) ) {
		ERROR("Not enough memory to check CCs\n");
		return;
	}

	i = 0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double pobs, pcalc;
		double res, corr, Ipart;
		Reflection *match;

		if ( !get_flag(refl) ) continue;  /* Not free-flagged */

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);
		if ( 2.0*res > crystal_get_resolution_limit(cr) ) continue;

		match = find_refl(full, h, k, l);
		if ( match == NULL ) continue;

		/* Calculated partiality */
		pcalc = get_partiality(refl);

		/* Observed partiality */
		corr = exp(-G) * exp(B*res*res) * get_lorentz(refl);
		Ipart = get_intensity(refl) * corr;
		pobs = Ipart / get_intensity(match);

		vec1[i] = pobs;
		vec2[i] = pcalc;
		i++;
	}

	cc = gsl_stats_correlation(vec1, 1, vec2, 1, i);
	//printf("%f\n", cc);
	if ( cc < 0.5 ) crystal_set_user_flag(cr, 6);

	free(vec1);
	free(vec2);
}


static void check_ccs(Crystal **crystals, int n_crystals, RefList *full)
{
	int i;

	for ( i=0; i<n_crystals; i++ ) {
		check_cc(crystals[i], full);
	}
}


void check_rejection(Crystal **crystals, int n, RefList *full)
{
	int i;
	int n_acc = 0;

	/* Check according to CCs FIXME: Disabled */
	//if ( full != NULL ) check_ccs(crystals, n, full);

	for ( i=0; i<n; i++ ) {

		/* Reject if B factor modulus is very large */
		if ( fabs(crystal_get_Bfac(crystals[i])) > 1e-17 ) {
			crystal_set_user_flag(crystals[i], 1);
		}

		if ( crystal_get_user_flag(crystals[i]) == 0 ) n_acc++;

	}

	if ( n_acc < 2 ) {
		ERROR("Not enough crystals left to proceed (%i).  Sorry.\n",
		      n_acc);
		exit(1);
	}
}
