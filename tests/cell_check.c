/*
 * cell_check.c
 *
 * Check that unit cells work correctly
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2014 Thomas White <taw@physics.org>
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
#include <fenv.h>

#include <cell.h>
#include <cell-utils.h>


int main(int argc, char *argv[])
{
	int fail = 0;
	struct quaternion orientation;
	UnitCell *cell;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double asmag, bsmag, csmag;
	double als, bes, gas;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	cell = cell_new_from_parameters(27.155e-9, 28.155e-9, 10.987e-9,
	                                deg2rad(90.0),
	                                deg2rad(90.0),
	                                deg2rad(120.0));
	if ( cell == NULL ) return 1;

	orientation = random_quaternion(rng);
	cell = cell_rotate(cell, orientation);

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	asmag = sqrt(pow(asx, 2.0) + pow(asy, 2.0) + pow(asz, 2.0));
	bsmag = sqrt(pow(bsx, 2.0) + pow(bsy, 2.0) + pow(bsz, 2.0));
	csmag = sqrt(pow(csx, 2.0) + pow(csy, 2.0) + pow(csz, 2.0));

	als = angle_between(bsx, bsy, bsz, csx, csy, csz);
	bes = angle_between(asx, asy, asz, csx, csy, csz);
	gas = angle_between(asx, asy, asz, bsx, bsy, bsz);

	STATUS("Separation between (100) planes = %5.2f nm\n", 1e9/asmag);
	STATUS("Separation between (010) planes = %5.2f nm\n", 1e9/bsmag);
	STATUS("Separation between (001) planes = %5.2f nm\n", 1e9/csmag);
	STATUS("Angle between (100) and (010) planes = %5.2f deg\n",
	       rad2deg(gas));
	STATUS("Angle between (100) and (001) planes = %5.2f deg\n",
	       rad2deg(bes));
	STATUS("Angle between (010) and (001) planes = %5.2f deg\n",
	       rad2deg(als));

	cell_free(cell);

	cell = cell_new();
	if ( cell == NULL ) return 1;

	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);
	cell_print(cell);

	cell_get_cartesian(cell, &ax, &ay, &az,
	                         &bx, &by, &bz,
	                         &cx, &cy, &cz);

	STATUS("Cell choice 1:\n");
	cell_set_cartesian(cell, ax, ay, az,
	                         bx, by, bz,
	                         cx, cy, cz);
	cell_print(cell);

	STATUS("Cell choice 2:\n");
	cell_set_cartesian(cell, bx, by, bz,
	                         -ax-bx, -ay-by, -az-bz,
	                         cx, cy, cz);
	cell_print(cell);


	STATUS("Cell choice 3:\n");
	cell_set_cartesian(cell, -ax-bx, -ay-by, -az-bz,
	                         ax, ay, az,
	                         cx, cy, cz);
	cell_print(cell);

	gsl_rng_free(rng);

	return fail;
}
