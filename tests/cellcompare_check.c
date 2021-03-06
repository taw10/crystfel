/*
 * cellcompare_check.c
 *
 * Check that unit cell comparison works
 *
 * Copyright © 2019-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
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
#include <gsl/gsl_randist.h>

#include <cell.h>
#include <cell-utils.h>
#include <utils.h>


static void complain(UnitCell *cell, UnitCell *cref, const char *t, const char *c)
{
	STATUS("These cells should %sbe the same%s:\n", t, c);
	STATUS("Transformed: ----------------------------\n");
	cell_print_full(cell);
	STATUS("Original: ---------------------------\n");
	cell_print_full(cref);
}


static RationalMatrix *random_derivative(gsl_rng *rng)
{
	int i, j;
	RationalMatrix *tr;

	tr = rtnl_mtx_new(3, 3);
	if ( tr == NULL ) return NULL;

	do {
		for ( i=0; i<3; i++ ) {
			for ( j=0; j<3; j++ ) {
				/* 0..6 inclusive -> -3..3 inclusive */
				signed int a = gsl_rng_uniform_int(rng, 7) - 3;
				/* 0..2 inclusive */
				signed int b = gsl_rng_uniform_int(rng, 3);
				if ( b == 0 ) {
					a = 0;
					b = 1;
				}
				rtnl_mtx_set(tr, i, j, rtnl(a, b));
			}
		}

	} while ( rtnl_cmp(rtnl_mtx_det(tr), rtnl_zero()) == 0 );

	return tr;
}


static IntegerMatrix *random_permutation(gsl_rng *rng)
{
	int i;
	IntegerMatrix *tr;
	int v[] = {0, 1, 2};

	tr = intmat_new(3, 3);
	if ( tr == NULL ) return NULL;

	gsl_ran_shuffle(rng, v, 3, sizeof(int));
	for ( i=0; i<3; i++ ) {
		intmat_set(tr, i, v[i], 1);
	}

	return tr;
}


static int check_ccp(UnitCell *cell, UnitCell *cref, double *tols,
                     int should_match)
{
	const char *a;
	const char *b;

	a = should_match ? "" : "NOT ";
	b = " with compare_cell_parameters";

	if ( compare_cell_parameters(cell, cref, tols) != should_match )
	{
		complain(cell, cref, a, b);
		return 1;
	}
	return 0;
}


static int check_ccpao(UnitCell *cell, UnitCell *cref, double *tols,
                       int should_match)
{
	const char *a;
	const char *b;

	a = should_match ? "" : "NOT ";
	b = " with compare_cell_parameters_and_orientation";

	if ( compare_cell_parameters_and_orientation(cell, cref,
	                                             tols) != should_match )
	{
		complain(cell, cref, a, b);
		return 1;
	}
	return 0;
}


static int check_cpcpao(UnitCell *cell, UnitCell *cref, double *tols,
                        int should_match)
{
	IntegerMatrix *m = NULL;
	const char *a;
	const char *b;

	a = should_match ? "" : "NOT ";
	b = " with compare_permuted_cell_parameters_and_orientation";

	if ( compare_permuted_cell_parameters_and_orientation(cell, cref, tols, &m)
	                        != should_match )
	{
		complain(cell, cref, a, b);
		intmat_free(m);
		return 1;
	}
	intmat_free(m);
	return 0;
}


static int check_cdcp(UnitCell *cell, UnitCell *cref, double *tols,
                      RationalMatrix *tr, int should_match)
{
	RationalMatrix *m = NULL;
	const char *a;
	const char *b;

	a = should_match ? "" : "NOT ";
	b = " with compare_derivative_cell_parameters";

	if ( compare_derivative_cell_parameters(cref, cell, tols, 1, &m) != should_match )
	{
		complain(cell, cref, a, b);
		STATUS("The transformation matrix to create the cell was:\n");
		rtnl_mtx_print(tr);
		STATUS("Cell comparison says do this to the reference to "
		       "create the cell:\n");
		rtnl_mtx_print(m);
		rtnl_mtx_free(m);
		return 1;
	}
	rtnl_mtx_free(m);
	return 0;
}


static int check_crcp(UnitCell *cell, UnitCell *cref, double *tols,
                      RationalMatrix *tr, int should_match)
{
	RationalMatrix *m = NULL;
	const char *a;
	const char *b;
	UnitCell *match;

	a = should_match ? "" : "NOT ";
	b = " with compare_reindexed_cell_parameters";

	match = compare_reindexed_cell_parameters(cell, cref, tols, &m);

	if ( (match != NULL) != should_match )
	{
		complain(cell, cref, a, b);
		STATUS("The transformation matrix to create the cell was:\n");
		rtnl_mtx_print(tr);
		STATUS("Cell comparison says do this to the reference to "
		       "create the cell:\n");
		rtnl_mtx_print(m);
		rtnl_mtx_free(m);
		return 1;
	}
	rtnl_mtx_free(m);
	return 0;
}


static int test_plain_rotation(UnitCell *cref, gsl_rng *rng,
                               double *tols)
{
	UnitCell *cell;

	cell = cell_rotate(cref, random_quaternion(rng));
	if ( cell == NULL ) return 1;

	if ( check_ccp(cell, cref, tols, 1) ) return 1;
	if ( check_ccpao(cell, cref, tols, 0) ) return 1;
	if ( check_cpcpao(cell, cref, tols, 0) ) return 1;
	if ( check_cdcp(cell, cref, tols, NULL, 1) ) return 1;
	if ( check_crcp(cell, cref, tols, NULL, 1) ) return 1;

	cell_free(cell);

	return 0;
}


static int test_permutation_no_rotation(UnitCell *cref, gsl_rng *rng,
                                        double *tols)
{
	IntegerMatrix *tr;
	UnitCell *cell;

	tr = random_permutation(rng);
	cell = cell_transform_intmat(cref, tr);

	if ( check_ccp(cell, cref, tols, intmat_is_identity(tr)) ) return 1;
	if ( check_ccpao(cell, cref, tols, intmat_is_identity(tr)) ) return 1;
	if ( check_cpcpao(cell, cref, tols, 1) ) return 1;
	if ( check_cdcp(cell, cref, tols, NULL, 1) ) return 1;
	if ( check_crcp(cell, cref, tols, NULL, 1) ) return 1;

	cell_free(cell);
	intmat_free(tr);

	return 0;
}


static int test_rotation_and_permutation(UnitCell *cref, gsl_rng *rng,
                                         double *tols)
{
	IntegerMatrix *tr;
	UnitCell *cell;
	UnitCell *cell2;

	cell2 = cell_rotate(cref, random_quaternion(rng));
	if ( cell2 == NULL ) return 1;

	tr = random_permutation(rng);
	cell = cell_transform_intmat(cell2, tr);
	cell_free(cell2);

	if ( check_ccp(cell, cref, tols, intmat_is_identity(tr)) ) return 1;
	if ( check_ccpao(cell, cref, tols, 0) ) return 1;
	if ( check_cpcpao(cell, cref, tols, 0) ) return 1;
	if ( check_cdcp(cell, cref, tols, NULL, 1) ) return 1;
	if ( check_crcp(cell, cref, tols, NULL, 1) ) return 1;

	cell_free(cell);
	intmat_free(tr);

	return 0;
}


static int test_derivative_lattice(UnitCell *cref, gsl_rng *rng,
                                   double *tols)
{
	RationalMatrix *tr;
	UnitCell *cell;

	cell = NULL;
	tr = NULL;
	do {
		cell_free(cell);
		rtnl_mtx_free(tr);
		tr = random_derivative(rng);
		cell = cell_transform_rational(cref, tr);
	} while ( (cell_get_centering(cell) == '?')
	       || (cell_get_centering(cell) == 'H' ) );
	/* H centering is no good because it needs a unique axis to
	 * be specified in order for uncentering in c_r_c_p to work.
	 * cell_transform_rational doesn't set the unique axis (yet?) */

	if ( check_ccp(cell, cref, tols, rtnl_mtx_is_identity(tr)) ) return 1;
	if ( check_ccpao(cell, cref, tols, rtnl_mtx_is_identity(tr)) ) return 1;
	if ( check_cpcpao(cell, cref, tols, rtnl_mtx_is_perm(tr)) ) return 1;
	if ( check_cdcp(cell, cref, tols, tr, 1) ) return 1;
	/* check_crcp: Sometimes yes, hard to tell */

	cell_free(cell);
	rtnl_mtx_free(tr);

	return 0;
}


static int test_derivative_lattice_rotation(UnitCell *cref, gsl_rng *rng,
                                            double *tols)
{
	RationalMatrix *tr;
	UnitCell *cell;
	UnitCell *cell2;

	cell2 = cell_rotate(cref, random_quaternion(rng));
	if ( cell2 == NULL ) return 1;

	cell = NULL;
	tr = NULL;
	do {
		cell_free(cell);
		rtnl_mtx_free(tr);
		tr = random_derivative(rng);
		cell = cell_transform_rational(cell2, tr);
	} while ( (cell_get_centering(cell) == '?')
	       || (cell_get_centering(cell) == 'H' ) );  /* See above */
	cell_free(cell2);

	if ( check_ccp(cell, cref, tols, rtnl_mtx_is_identity(tr)) ) return 1;
	if ( check_ccpao(cell, cref, tols, 0) ) return 1;
	if ( check_cpcpao(cell, cref, tols, 0) ) return 1;
	if ( check_cdcp(cell, cref, tols, tr, 1) ) return 1;
	/* check_crcp: Sometimes yes, hard to tell */

	cell_free(cell);
	rtnl_mtx_free(tr);

	return 0;
}


extern IntegerMatrix *reduce_g6(struct g6 g, double epsrel);

int main(int argc, char *argv[])
{
	UnitCell *cref;
	gsl_rng *rng;
	int fail = 0;
	double tols[] = { 0.01, 0.01, 0.01,
	                  deg2rad(1.0), deg2rad(1.0), deg2rad(1.0) };

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	cref = cell_new_from_parameters(3e-10, 5.196e-10, 2e-10,
	                                deg2rad(103.9166666),
	                                deg2rad(109.4666666),
	                                deg2rad(134.8833333));
	cell_set_centering(cref, 'P');
	if ( cref == NULL ) return 1;

	fail += test_plain_rotation(cref, rng, tols);
	fail += test_permutation_no_rotation(cref, rng, tols);
	fail += test_rotation_and_permutation(cref, rng, tols);
	fail += test_derivative_lattice(cref, rng, tols);
	fail += test_derivative_lattice_rotation(cref, rng, tols);

	cell_free(cref);
	gsl_rng_free(rng);

	return fail;
}
