/*
 * impose_bravais.c
 *
 * Test imposition of Bravais conditions
 *
 * Copyright © 2025 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2025 Thomas White <taw@physics.org>
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

#include <cell.h>
#include <cell-utils.h>


static int test(UnitCell *in, LatticeType latt, char ua,
                double a1, double b1, double c1,
                double al1, double be1, double ga1)
{
	UnitCell *out;
	double a2, b2, c2, al2, be2, ga2;
	int fail = 0;

	out = impose_bravais(in, latt, ua);
	cell_print_full(out);
	cell_get_parameters(out, &a2, &b2, &c2, &al2, &be2, &ga2);
	if ( fabs(a1*1e-10 - a2) > 1e-13 ) { ERROR("a axis length wrong\n"); fail = 1; }
	if ( fabs(b1*1e-10 - b2) > 1e-13 ) { ERROR("b axis length wrong\n"); fail = 1; }
	if ( fabs(c1*1e-10 - c2) > 1e-13 ) { ERROR("c axis length wrong\n"); fail = 1; }
	if ( fabs(deg2rad(al1) - al2) > deg2rad(0.01) ) { ERROR("al angle wrong\n"); fail = 1; }
	if ( fabs(deg2rad(be1) - be2) > deg2rad(0.01) ) { ERROR("be angle wrong\n"); fail = 1; }
	if ( fabs(deg2rad(ga1) - ga2) > deg2rad(0.01) ) { ERROR("ga angle wrong\n"); fail = 1; }
	cell_free(out);

	return fail;
}


int main(int argc, char *argv[])
{
	UnitCell *cell;
	int fail = 0;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	STATUS("Orthorhombic cell in reference orientation\n");
	cell = cell_new_from_parameters(123e-10, 45e-10, 67e-10,
	                                deg2rad(89.9), deg2rad(90.1), deg2rad(90.0));
	fail += test(cell, L_ORTHORHOMBIC, '*', 123, 45, 67, 90, 90, 90);

	STATUS("\n\nOrthorhombic cell in random orientation\n");
	cell = cell_rotate(cell, random_quaternion(rng));
	fail += test(cell, L_ORTHORHOMBIC, '*', 123, 45, 67, 90, 90, 90);

	STATUS("\n\nRhombohedral cell in reference orientation\n");
	cell = cell_new_from_parameters(120e-10, 121e-10, 122e-10,
	                                deg2rad(35.0), deg2rad(35.1), deg2rad(35.2));
	fail += test(cell, L_RHOMBOHEDRAL, '*', 121, 121, 121, 35.1, 35.1, 35.1);

	STATUS("\n\nRhombohedral cell in random orientation\n");
	cell = cell_rotate(cell, random_quaternion(rng));
	fail += test(cell, L_RHOMBOHEDRAL, '*', 121, 121, 121, 35.1, 35.1, 35.1);

	return fail;
}
