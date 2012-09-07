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


static int check_cell(UnitCell *cell, const char *text)
{
	int err = 0;

	if ( !cell_is_sensible(cell) ) {
		ERROR(" %s unit cell parameters are not sensible.\n", text);
		err = 1;
	}

	if ( !bravais_lattice(cell) ) {
		ERROR(" %s unit cell is not a conventional Bravais"
		      " lattice.\n", text);
		err = 1;
	}

	if ( !right_handed(cell) ) {
		ERROR(" %s unit cell is not right handed.\n", text);
		err = 1;
	}

	if ( err ) cell_print(cell);

	return err;
}


static int check_centering(double a, double b, double c,
                           double al, double be, double ga,
                           LatticeType latt, char cen, char ua)
{
	UnitCell *cell;
	UnitCell *n;
	UnitCellTransformation *t;
	int fail = 0;

	STATUS("Checking %s %c (ua %c) %5.2e %5.2e %5.2e %5.2f %5.2f %5.2f\n",
	       str_lattice(latt), cen, ua, a, b, c, al, be, ga);

	cell = cell_new_from_parameters(a, b, c,
	                                deg2rad(al), deg2rad(be), deg2rad(ga));
	cell_set_lattice_type(cell, latt);
	cell_set_centering(cell, cen);
	cell_set_unique_axis(cell, ua);

	if ( check_cell(cell, "Input") ) fail = 1;
	//cell_print(cell);
	n = uncenter_cell(cell, &t);
	if ( n != NULL ) {
		if ( check_cell(n, "Output") ) fail = 1;
	} else {
		fail = 1;
	}

	STATUS("Transformation was:\n");
	tfn_print(t);

	if ( fail ) ERROR("\n");

	return fail;
}



int main(int argc, char *argv[])
{
	int fail = 0;

	/* Triclinic P */
	fail += check_centering(50e-10, 55e-10, 70e-10, 67.0, 70.0, 77.0,
	                        L_TRICLINIC, 'P', '*');

	/* Monoclinic P */
	fail += check_centering(10e-10, 20e-10, 30e-10, 100.0, 90.0, 90.0,
	                        L_MONOCLINIC, 'P', 'a');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 100.0, 90.0,
	                        L_MONOCLINIC, 'P', 'b');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 100.0,
	                        L_MONOCLINIC, 'P', 'c');

	/* Monoclinic A */
	fail += check_centering(10e-10, 20e-10, 30e-10, 100.0, 90.0, 90.0,
	                        L_MONOCLINIC, 'A', 'a');

	/* Monoclinic B */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 100.0, 90.0,
	                        L_MONOCLINIC, 'B', 'b');

	/* Monoclinic C */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 100.0,
	                        L_MONOCLINIC, 'C', 'c');

	/* Orthorhombic P */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'P', '*');

	/* Orthorhombic A */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'A', '*');

	/* Orthorhombic B */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'B', '*');

	/* Orthorhombic C */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'C', '*');

	/* Orthorhombic I */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'I', '*');

	/* Orthorhombic F */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'F', '*');

	/* Tetragonal P */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'P', 'a');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'P', 'b');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'P', 'c');

	/* Tetragonal I */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'I', 'a');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'I', 'b');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'I', 'c');

	/* Rhombohedral R */
	fail += check_centering(10e-10, 10e-10, 10e-10, 60.0, 60.0, 60.0,
	                        L_RHOMBOHEDRAL, 'R', '*');

	/* Hexagonal P */
	fail += check_centering(10e-10, 20e-10, 30e-10, 120.0, 90.0, 90.0,
	                        L_HEXAGONAL, 'P', 'a');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 120.0, 90.0,
	                        L_HEXAGONAL, 'P', 'b');
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 120.0,
	                        L_HEXAGONAL, 'P', 'c');

	/* Hexagonal H (PDB-speak for rhombohedral) */
	fail += check_centering(40e-10, 20e-10, 20e-10, 120.0, 90.0, 90.0,
	                        L_HEXAGONAL, 'H', 'a');
	fail += check_centering(20e-10, 40e-10, 20e-10, 90.0, 120.0, 90.0,
	                        L_HEXAGONAL, 'H', 'b');
	fail += check_centering(20e-10, 20e-10, 40e-10, 90.0, 90.0, 120.0,
	                        L_HEXAGONAL, 'H', 'c');

	/* Cubic P */
	fail += check_centering(30e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_CUBIC, 'P', '*');

	/* Cubic I */
	fail += check_centering(30e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_CUBIC, 'I', '*');

	/* Cubic F */
	fail += check_centering(30e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_CUBIC, 'F', '*');

	return fail;
}
