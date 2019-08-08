/*
 * cellcompare_check.c
 *
 * Check that unit cell comparison works
 *
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


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


static RationalMatrix *random_reindexing(gsl_rng *rng)
{
	int i;
	RationalMatrix *tr;
	int v[] = {0, 1, 2};

	tr = rtnl_mtx_new(3, 3);
	if ( tr == NULL ) return NULL;

	gsl_ran_shuffle(rng, v, 3, sizeof(int));
	for ( i=0; i<3; i++ ) {
		rtnl_mtx_set(tr, i, v[i], rtnl(1, 1));
	}

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


static int check_cpcp(UnitCell *cell, UnitCell *cref, double *tols,
                      int should_match)
{
	IntegerMatrix *m = NULL;
	const char *a;
	const char *b;

	a = should_match ? "" : "NOT ";
	b = " with compare_permuted_cell_parameters";

	if ( compare_permuted_cell_parameters(cell, cref, tols, &m) != should_match )
	{
		complain(cell, cref, a, b);
		STATUS("Matrix was:\n");
		intmat_print(m);
		intmat_free(m);
		return 1;
	}
	intmat_free(m);
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
	                                             tols[0], tols[3]) != should_match )
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

	if ( compare_permuted_cell_parameters_and_orientation(cell, cref,
	                                                      tols[0], tols[3], &m) != should_match )
	{
		complain(cell, cref, a, b);
		intmat_free(m);
		return 1;
	}
	intmat_free(m);
	return 0;
}


static int check_crcp(UnitCell *cell, UnitCell *cref, double *tols,
                      int should_match)
{
	RationalMatrix *m = NULL;
	const char *a;
	const char *b;

	a = should_match ? "" : "NOT ";
	b = " with compare_reindexed_cell_parameters";

	if ( compare_reindexed_cell_parameters(cell, cref, tols, 0, &m) != should_match )
	{
		complain(cell, cref, a, b);
		rtnl_mtx_free(m);
		return 1;
	}
	rtnl_mtx_free(m);
	return 0;
}


int main(int argc, char *argv[])
{
	UnitCell *cell, *cref;
	gsl_rng *rng;
	int i;
	double tols[] = { 0.01, 0.01, 0.01,
	                  deg2rad(1.0), deg2rad(1.0), deg2rad(1.0) };

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	cref = cell_new_from_parameters(50e-10, 55e-10, 70e-10,
	                                deg2rad(67.0),
	                                deg2rad(70.0),
	                                deg2rad(77.0));
	if ( cref == NULL ) return 1;

	/* Just rotate cell */
	STATUS("Testing plain rotation...\n");
	for ( i=0; i<100; i++ ) {

		cell = cell_rotate(cref, random_quaternion(rng));
		if ( cell == NULL ) return 1;

		if ( check_ccp(cell, cref, tols, 1) ) return 1;
		if ( check_cpcp(cell, cref, tols, 1) ) return 1;
		if ( check_ccpao(cell, cref, tols, 0) ) return 1;
		if ( check_cpcpao(cell, cref, tols, 0) ) return 1;
		if ( check_crcp(cell, cref, tols, 1) ) return 1;

		cell_free(cell);
	}

	/* Permute axes but don't rotate */
	STATUS("Testing axis permutation...\n");
	for ( i=0; i<100; i++ ) {

		IntegerMatrix *tr;

		tr = random_permutation(rng);
		cell = cell_transform_intmat(cref, tr);

		if ( check_ccp(cell, cref, tols, intmat_is_identity(tr)) ) return 1;
		if ( check_cpcp(cell, cref, tols, 1) ) return 1;
		if ( check_ccpao(cell, cref, tols, intmat_is_identity(tr)) ) return 1;
		if ( check_cpcpao(cell, cref, tols, 1) ) return 1;
		if ( check_crcp(cell, cref, tols, 1) ) return 1;

		cell_free(cell);
		intmat_free(tr);
	}

	/* Rotate cell and permute axes */
	STATUS("Testing rotation with axis permutation...\n");
	for ( i=0; i<100; i++ ) {

		IntegerMatrix *tr;
		UnitCell *cell2;

		cell2 = cell_rotate(cref, random_quaternion(rng));
		if ( cell2 == NULL ) return 1;

		tr = random_permutation(rng);
		cell = cell_transform_intmat(cell2, tr);
		cell_free(cell2);

		if ( check_ccp(cell, cref, tols, intmat_is_identity(tr)) ) return 1;
		if ( check_cpcp(cell, cref, tols, 1) ) return 1;
		if ( check_ccpao(cell, cref, tols, 0) ) return 1;
		if ( check_cpcpao(cell, cref, tols, 0) ) return 1;
		if ( check_crcp(cell, cref, tols, 1) ) return 1;

		cell_free(cell);
		intmat_free(tr);
	}

	/* Reindex */
	STATUS("Testing reindexing...\n");
	for ( i=0; i<100; i++ ) {

		RationalMatrix *tr;

		tr = random_reindexing(rng);
		cell = cell_transform_rational(cref, tr);

		if ( check_ccp(cell, cref, tols, rtnl_mtx_is_identity(tr)) ) return 1;
		if ( check_cpcp(cell, cref, tols, rtnl_mtx_is_perm(tr)) ) return 1;
		if ( check_ccpao(cell, cref, tols, rtnl_mtx_is_identity(tr)) ) return 1;
		if ( check_cpcpao(cell, cref, tols, rtnl_mtx_is_perm(tr)) ) return 1;
		if ( check_crcp(cell, cref, tols, 1) ) return 1;

		cell_free(cell);
		rtnl_mtx_free(tr);
	}

	/* Reindex and rotate */
	STATUS("Testing reindexing with rotation...\n");
	for ( i=0; i<100; i++ ) {

		RationalMatrix *tr;
		UnitCell *cell2;

		cell2 = cell_rotate(cref, random_quaternion(rng));
		if ( cell2 == NULL ) return 1;

		tr = random_reindexing(rng);
		cell = cell_transform_rational(cell2, tr);
		cell_free(cell2);

		if ( check_ccp(cell, cref, tols, rtnl_mtx_is_identity(tr)) ) return 1;
		if ( check_cpcp(cell, cref, tols, rtnl_mtx_is_perm(tr)) ) return 1;
		if ( check_ccpao(cell, cref, tols, 0) ) return 1;
		if ( check_cpcpao(cell, cref, tols, 0) ) return 1;
		if ( check_crcp(cell, cref, tols, 1) ) return 1;

		cell_free(cell);
		rtnl_mtx_free(tr);
	}

	/* NB There's no compare_reindexed_cell_parameters_and_orientation */

	/* Test using vectors */

	cell_free(cref);
	gsl_rng_free(rng);

	return 0;
}
