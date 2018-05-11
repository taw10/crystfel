/*
 * scaling_check.c
 *
 * Check that scaling works
 *
 * Copyright © 2017-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2017-2018 Thomas White <taw@physics.org>
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
#include <cell-utils.h>

#include "../src/scaling.h"

int test_scaling(double G, double B, int scaleflags, int do_partials,
                 gsl_rng *rng)
{
	int i;
	Crystal *cr;
	RefList *list1;
	RefList *list2;
	int r;
	UnitCell *cell;

	list1 = reflist_new();
	list2 = reflist_new();

	cell = cell_new();
	cell_set_parameters(cell, 50e-10, 50e-10, 50e-10,
	                    deg2rad(90), deg2rad(90), deg2rad(90));

	for ( i=0; i<50; i++ ) {

		signed int h, k, l;
		Reflection *refl1;
		Reflection *refl2;
		double intens, p, s, L;

		h = gsl_rng_uniform_int(rng, 20) - gsl_rng_uniform_int(rng, 40);
		k = gsl_rng_uniform_int(rng, 20) - gsl_rng_uniform_int(rng, 40);
		l = gsl_rng_uniform_int(rng, 20) - gsl_rng_uniform_int(rng, 40);

		refl1 = add_refl(list1, h, k, l);
		refl2 = add_refl(list2, h, k, l);
		intens =  gsl_rng_uniform(rng);  /* [0,1) */
		p = do_partials ? gsl_rng_uniform(rng) : 1.0;
		L = gsl_rng_uniform(rng);

		s = resolution(cell, h, k, l);

		/* Reference */
		set_intensity(refl2, intens);
		set_partiality(refl2, 1.0);
		set_lorentz(refl2, 1.0);
		set_redundancy(refl2, 2);

		/* Crystal */
		set_intensity(refl1, intens * G * exp(-B*s*s) * p / L);
		set_partiality(refl1, p);
		set_lorentz(refl1, L);

	}

	cr = crystal_new();
	crystal_set_reflections(cr, list1);
	crystal_set_cell(cr, cell);

	crystal_set_osf(cr, 999.0);
	crystal_set_Bfac(cr, 999.0);

	r = scale_one_crystal(cr, list2, scaleflags | SCALE_VERBOSE_ERRORS);
	STATUS("Scaling result: %i, G = %8.4f, B = %8.4f A^2\n", r,
	       crystal_get_osf(cr), crystal_get_Bfac(cr)*1e20);

	if ( fabs(G - crystal_get_osf(cr)) > 0.001 ) r = 1;
	if ( fabs(B - crystal_get_Bfac(cr)) > 0.001e-20 ) r = 1;

	reflist_free(list1);
	reflist_free(list2);
	cell_free(cell);
	crystal_free(cr);

	if ( r ) {
		STATUS("  (should be: G = %8.4f, B = %8.4f A^2), %s partials\n",
		       G, B*1e20, do_partials ? "with" : "no");
	}

	return r;
}


int main(int argc, char *argv[])
{
	int fail = 0;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	fail += test_scaling(2.0, 0.0, SCALE_NO_B, 0, rng);
	fail += test_scaling(2.0, 0.0, SCALE_NONE, 0, rng);
	fail += test_scaling(2.0, 10.0e-20, SCALE_NONE, 0, rng);
	fail += test_scaling(5.0, 30.0e-20, SCALE_NONE, 0, rng);
	fail += test_scaling(2.0, 0.0, SCALE_NO_B, 1, rng);
	fail += test_scaling(2.0, 0.0, SCALE_NONE, 1, rng);
	fail += test_scaling(2.0, 10.0e-20, SCALE_NONE, 1, rng);
	fail += test_scaling(5.0, 30.0e-20, SCALE_NONE, 1, rng);

	gsl_rng_free(rng);

	return fail;
}
