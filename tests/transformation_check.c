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


int main(int argc, char *argv[])
{
	int fail = 0;
	UnitCell *cell, *cnew, *cback;
	UnitCellTransformation *tfn, *inv;

	cell = cell_new_from_parameters(50e-10, 55e-10, 70e-10,
	                                deg2rad(67.0),
	                                deg2rad(70.0),
	                                deg2rad(77.0));
	if ( cell == NULL ) return 1;

	tfn = tfn_identity();
	if ( tfn == NULL ) return 1;

	tfn_combine(tfn, tfn_vector(0,1,0),
	                 tfn_vector(1,0,0),
	                 tfn_vector(0,0,1));


	cell_print(cell);
	tfn_print(tfn);

	cnew = cell_transform(cell, tfn);
	cell_print(cnew);

	cback = cell_transform_inverse(cnew, tfn);
	inv = tfn_inverse(tfn);
	tfn_print(inv);
	cell_print(cback);

	cell_get_cartesian(cell, &ax1, &ay1, &az1,
	                         &by1, &by1, &bz1,
	                         &cx1, &cy1, &cz1);
	cell_get_cartesian(cback, &ax, &ay, &az, &by, &by, &bz, &cx, &cy, &cz);

	return fail;
}
