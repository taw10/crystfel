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

#include "crystal.h"
#include "reflist.h"
#include "rejection.h"


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
	check_rejection(crystals, n);
}


void check_rejection(Crystal **crystals, int n)
{
	int i;
	int n_acc = 0;

	for ( i=0; i<n; i++ ) {
		if ( crystal_get_user_flag(crystals[i]) == 0 ) {
			n_acc++;
			if ( n_acc >= 2 ) break;
		}
	}

	if ( n_acc < 2 ) {
		ERROR("Not enough crystals left to proceed (%i).  Sorry.\n",
		      n_acc);
		exit(1);
	}
}
