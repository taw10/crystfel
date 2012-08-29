/*
 * centering_check.c
 *
 * Check that centering of cells works
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


static int check_cell(UnitCell *cell)
{
	int err = 0;

	if ( !cell_is_sensible(cell) ) {
		ERROR("Warning: Unit cell parameters are not sensible.\n");
		err = 1;
	}

	if ( !bravais_lattice(cell) ) {
		ERROR("Warning: Unit cell is not a conventional Bravais"
		      " lattice.\n");
		err = 1;
	}

	if ( !right_handed(cell) ) {
		ERROR("Warning: Unit cell is not right handed.\n");
		err = 1;
	}

	return err;
}


static int check_centering()
{
	UnitCell *cell;
	UnitCell *n;
	int fail = 0;

	cell = cell_new();

	cell_set_parameters(cell, 10e-10, 20e-10, 30e-10,
	                    deg2rad(90.0), deg2rad(90.0), deg2rad(100.0));
	cell_set_centering(cell, 'C');
	cell_set_lattice_type(cell, L_MONOCLINIC);
	cell_set_unique_axis(cell, 'c');
	if ( check_cell(cell) ) fail = 1;
	n = uncenter_cell(cell);
	if ( check_cell(n) ) fail = 1;

	return fail;
}



int main(int argc, char *argv[])
{
	int fail = 0;

	/* Triclinic P */
	check_centering(10e-10, 20e-10, 30e-10, 100.0, 120.0, 140.0,
	                L_TRICLINIC, 'P', '*',
	                );
	/* Monoclinic P */
	/* Monoclinic A */
	/* Monoclinic B */
	/* Monoclinic C */
	/* Orthorhombic P */
	/* Orthorhombic A */
	/* Orthorhombic B */
	/* Orthorhombic C */
	/* Orthorhombic I */
	/* Orthorhombic F */
	/* Tetragonal P */
	/* Tetragonal I */
	/* Rhombohedral R */
	/* Hexagonal P */
	/* Hexagonal H (PDB-speak for rhombohedral) */
	/* Cubic P */
	/* Cubic I */
	/* Cubic F */


	return fail;
}
