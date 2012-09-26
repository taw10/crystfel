/*
 * transformation_check.c
 *
 * Check that unit cell transformations work
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
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

#include <cell.h>
#include <cell-utils.h>


static int check_transformation(UnitCell *cell, UnitCellTransformation *tfn)
{
	UnitCell *cnew, *cback;
	UnitCellTransformation *inv;
	double a[9], b[9];
	int i;
	int fail = 0;

	cnew = cell_transform(cell, tfn);

	cback = cell_transform_inverse(cnew, tfn);

	cell_get_cartesian(cell, &a[0], &a[1], &a[2],
	                         &a[3], &a[4], &a[5],
	                         &a[6], &a[7], &a[8]);
	cell_get_cartesian(cback, &b[0], &b[1], &b[2],
	                          &b[3], &b[4], &b[5],
	                          &b[6], &b[7], &b[8]);
	for ( i=0; i<9; i++ ) {
		if ( !within_tolerance(a[i], b[i], 0.1) ) {
			fail = 1;
			STATUS("%e %e\n", a[i], b[i]);
		}
	}

	if ( fail ) {
		ERROR("Original cell not recovered after transformation:\n");
		cell_print(cell);
		tfn_print(tfn);
		inv = tfn_inverse(tfn);
		tfn_print(inv);
		cell_print(cback);
	}

	return fail;
}


int main(int argc, char *argv[])
{
	int fail = 0;
	UnitCell *cell, *cref;
	UnitCellTransformation *tfn;

	cref = cell_new_from_parameters(50e-10, 55e-10, 70e-10,
	                                deg2rad(67.0),
	                                deg2rad(70.0),
	                                deg2rad(77.0));
	if ( cref == NULL ) return 1;

	cell = cell_rotate(cref, random_quaternion());
	if ( cell == NULL ) return 1;
	cell_free(cref);

	/* Permutation of axes */
	tfn = tfn_identity();
	if ( tfn == NULL ) return 1;
	tfn_combine(tfn, tfn_vector(0,1,0),
	                 tfn_vector(0,0,1),
	                 tfn_vector(1,0,0));
	fail += check_transformation(cell, tfn);
	tfn_free(tfn);

	/* Doubling of cell in one direction */
	tfn = tfn_identity();
	if ( tfn == NULL ) return 1;
	tfn_combine(tfn, tfn_vector(2,0,0),
	                 tfn_vector(0,1,0),
	                 tfn_vector(0,0,1));
	fail += check_transformation(cell, tfn);
	tfn_free(tfn);

	/* Diagonal supercell */
	tfn = tfn_identity();
	if ( tfn == NULL ) return 1;
	tfn_combine(tfn, tfn_vector(1,1,0),
	                 tfn_vector(0,1,0),
	                 tfn_vector(0,0,1));
	fail += check_transformation(cell, tfn);
	tfn_free(tfn);

	/* Crazy */
	tfn = tfn_identity();
	if ( tfn == NULL ) return 1;
	tfn_combine(tfn, tfn_vector(1,1,0),
	                 tfn_vector(0,1,0),
	                 tfn_vector(0,1,1));
	fail += check_transformation(cell, tfn);
	tfn_free(tfn);

	cell_free(cell);

	return fail;
}
