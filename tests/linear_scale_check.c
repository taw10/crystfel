/*
 * linear_scale_check.c
 *
 * Check that linear scaling works
 *
 * Copyright © 2017 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017 Thomas White <taw@physics.org>
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
#include <stdio.h>

#include <reflist.h>

#include "../src/scaling.h"

int main(int argc, char *argv[])
{
	int fail = 0;
	int i;
	gsl_rng *rng;
	RefList *list1;
	RefList *list2;
	double G;
	int r;

	list1 = reflist_new();
	list2 = reflist_new();

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	for ( i=0; i<50; i++ ) {
		signed int h, k, l;
		Reflection *refl1;
		Reflection *refl2;
		double intens;
		h = gsl_rng_uniform_int(rng, 20) - gsl_rng_uniform_int(rng, 40);
		k = gsl_rng_uniform_int(rng, 20) - gsl_rng_uniform_int(rng, 40);
		l = gsl_rng_uniform_int(rng, 20) - gsl_rng_uniform_int(rng, 40);
		refl1 = add_refl(list1, h, k, l);
		refl2 = add_refl(list2, h, k, l);
		intens =  gsl_rng_uniform(rng);  /* [0,1) */
		set_intensity(refl1, intens);
		set_partiality(refl1, 1.0);
		set_intensity(refl2, intens*2.0);
		set_partiality(refl2, 1.0);
	}

	r = linear_scale(list1, list2, &G);
	STATUS("Scaling result: %i, G = %f\n", r, G);

	return fail;
}
