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


static double moduli_check(double ax, double ay, double az,
                           double bx, double by, double bz)
{
	double ma = modulus(ax, ay, az);
	double mb = modulus(bx, by, bz);
	return fabs(ma-mb)/ma;
}
static void complain(UnitCell *cell, UnitCell *cref, const char *t)
{
	STATUS("These cells should %s be the same:\n", t);
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


int compare_cell_parameters_and_orientation2(UnitCell *cell1, UnitCell *cell2,
                                            const double ltl, const double atl)
{
	double ax1, ay1, az1, bx1, by1, bz1, cx1, cy1, cz1;
	double ax2, ay2, az2, bx2, by2, bz2, cx2, cy2, cz2;

	if ( cell_get_centering(cell1) != cell_get_centering(cell2) ) return 0;

	cell_get_cartesian(cell1, &ax1, &ay1, &az1,
	                          &bx1, &by1, &bz1,
	                          &cx1, &cy1, &cz1);

	cell_get_cartesian(cell2, &ax2, &ay2, &az2,
	                          &bx2, &by2, &bz2,
	                          &cx2, &cy2, &cz2);

	cell_print_full(cell1);
	cell_print_full(cell2);

	STATUS("%f\n", rad2deg(atl));
	STATUS("%f\n", rad2deg(angle_between(ax1, ay1, az1, ax2, ay2, az2)));
	STATUS("%f\n", rad2deg(angle_between(bx1, by1, bz1, bx2, by2, bz2)));
	STATUS("%f\n", rad2deg(angle_between(cx1, cy1, cz1, cx2, cy2, cz2)));

	if ( angle_between(ax1, ay1, az1, ax2, ay2, az2) > atl ) return 0;
	if ( angle_between(bx1, by1, bz1, bx2, by2, bz2) > atl ) return 0;
	if ( angle_between(cx1, cy1, cz1, cx2, cy2, cz2) > atl ) return 0;

	if ( moduli_check(ax1, ay1, az1, ax2, ay2, az2) > ltl ) return 0;
	if ( moduli_check(bx1, by1, bz1, bx2, by2, bz2) > ltl ) return 0;
	if ( moduli_check(cx1, cy1, cz1, cx2, cy2, cz2) > ltl ) return 0;

	return 1;
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
		if ( !compare_cell_parameters(cell, cref, tols) ) {
			complain(cell, cref, "");
			return 1;
		}
		if ( compare_cell_parameters_and_orientation(cell, cref,
		                                             tols[0], tols[3]) )
		{
			complain(cell, cref, "NOT");
			return 1;
		}
		cell_free(cell);
	}

	/* Permute axes but don't rotate */
	STATUS("Testing axis permutation...\n");
	for ( i=0; i<100; i++ ) {

		IntegerMatrix *m;
		IntegerMatrix *tr;

		tr = random_permutation(rng);
		cell = cell_transform_intmat(cref, tr);

		if ( !compare_permuted_cell_parameters_and_orientation(cell, cref,
		                                                       tols[0], tols[3], &m) )
		{
			complain(cell, cref, "");
			return 1;
		}

		if ( compare_cell_parameters(cell, cref, tols)
		  && !intmat_is_identity(tr) )
		{
			complain(cell, cref, "NOT");
			return 1;
		}

		cell_free(cell);
		intmat_free(tr);
		intmat_free(m);
	}

	/* Rotate cell and permute axes */
	STATUS("Testing rotation with axis permutation...\n");
	for ( i=0; i<100; i++ ) {

		IntegerMatrix *m;
		IntegerMatrix *tr;
		UnitCell *cell2;

		cell = cell_rotate(cref, random_quaternion(rng));
		if ( cell == NULL ) return 1;

		tr = random_permutation(rng);
		cell2 = cell_transform_intmat(cell, tr);

		if ( !compare_permuted_cell_parameters(cell2, cref, tols, &m) ) {
			complain(cell2, cref, "");
			return 1;
		}

		if ( compare_permuted_cell_parameters_and_orientation(cell2, cref,
		                                                      tols[0], tols[3], &m) )
		{
			UnitCell *cc;
			complain(cell2, cref, "NOT, with just permutation,");
			STATUS("Matrix was (det=%i):\n", intmat_det(m));
			intmat_print(m);
			STATUS("Transformed version of cref:\n");
			cc = cell_transform_intmat(cref, m);
			cell_print_full(cc);
			cell_free(cc);
			return 1;
		}

		if ( compare_cell_parameters_and_orientation(cell2, cref,
		                                             tols[0], tols[3]) )
		{
			complain(cell2, cref, "NOT, without any change,");
			return 1;
		}

		cell_free(cell);
		cell_free(cell2);
		intmat_free(tr);
		intmat_free(m);
	}

	/* Reindex and rotate */
	STATUS("Testing rotation with reindexing...\n");
	for ( i=0; i<100; i++ ) {

		RationalMatrix *m;
		RationalMatrix *tr;
		UnitCell *cell2;

		cell = cell_rotate(cref, random_quaternion(rng));
		if ( cell == NULL ) return 1;

		tr = random_reindexing(rng);
		cell2 = cell_transform_rational(cell, tr);
		cell_free(cell);

		if ( !compare_reindexed_cell_parameters(cell2, cref, tols, 0, &m) ) {
			complain(cell2, cref, "");
			return 1;
		}
		cell_free(cell2);
		rtnl_mtx_free(tr);
		rtnl_mtx_free(m);
	}

	/* NB There's no compare_reindexed_cell_parameters_and_orientation */

	/* Test using vectors */

	cell_free(cref);
	gsl_rng_free(rng);

	return 0;
}
