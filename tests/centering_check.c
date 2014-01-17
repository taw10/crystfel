/*
 * centering_check.c
 *
 * Check that centering of cells works
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


static int check_cell(UnitCell *cell, const char *text)
{
	int err = 0;

	err = validate_cell(cell);

	if ( err ) {
		ERROR("%s cell:\n", text);
		cell_print(cell);
	}

	return err;
}


static int check_centering(double a, double b, double c,
                           double al, double be, double ga,
                           LatticeType latt, char cen, char ua, gsl_rng *rng)
{
	UnitCell *cell, *cref;
	UnitCell *n;
	UnitCellTransformation *t;
	int fail = 0;
	int i;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	STATUS("   --------------->  "
	       "Checking %s %c (ua %c) %5.2e %5.2e %5.2e %5.2f %5.2f %5.2f\n",
	       str_lattice(latt), cen, ua, a, b, c, al, be, ga);

	cref = cell_new_from_parameters(a, b, c,
	                                deg2rad(al), deg2rad(be), deg2rad(ga));
	cell_set_lattice_type(cref, latt);
	cell_set_centering(cref, cen);
	cell_set_unique_axis(cref, ua);

	cell = cell_rotate(cref, random_quaternion(rng));
	if ( cell == NULL ) return 1;
	cell_free(cref);

	check_cell(cell, "Input");
	n = uncenter_cell(cell, &t);
	if ( n != NULL ) {
		STATUS("Transformation was:\n");
		tfn_print(t);
		if ( check_cell(n, "Output") ) fail = 1;
		if ( !fail ) cell_print(n);
	} else {
		fail = 1;
	}

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	cell_get_cartesian(n, &ax, &ay, &az,
	                      &bx, &by, &bz,
	                      &cx, &cy, &cz);

	fesetround(1);  /* Round towards nearest */
	for ( i=0; i<100; i++ ) {

		signed int h, k, l;
		double x, y, z;
		double nh, nk, nl;
		double dh, dk, dl;
		int f = 0;

		do {

			h = gsl_rng_uniform_int(rng, 30);
			k = gsl_rng_uniform_int(rng, 30);
			l = gsl_rng_uniform_int(rng, 30);

		} while ( forbidden_reflection(cell, h, k, l) );

		x = h*asx + k*bsx + l*csx;
		y = h*asy + k*bsy + l*csy;
		z = h*asz + k*bsz + l*csz;

		nh = x*ax + y*ay + z*az;
		nk = x*bx + y*by + z*bz;
		nl = x*cx + y*cy + z*cz;

		dh = nh - lrint(nh);
		dk = nk - lrint(nk);
		dl = nl - lrint(nl);
		if ( fabs(dh) > 0.1 ) f++;
		if ( fabs(dk) > 0.1 ) f++;
		if ( fabs(dl) > 0.1 ) f++;

		if ( f ) {
			STATUS("Centered %3i %3i %3i -> "
			       "Primitive %7.2f %7.2f %7.2f\n",
			       h, k, l, nh, nk, nl);
			fail = 1;
		}

	}

	cell_get_reciprocal(n, &asx, &asy, &asz,
	                       &bsx, &bsy, &bsz,
	                       &csx, &csy, &csz);
	cell_get_cartesian(cell, &ax, &ay, &az,
	                         &bx, &by, &bz,
	                         &cx, &cy, &cz);

	for ( i=0; i<100; i++ ) {

		signed int h, k, l;
		double x, y, z;
		double nh, nk, nl;
		double dh, dk, dl;
		int f = 0;
		long int ih, ik, il;

		h = gsl_rng_uniform_int(rng, 30);
		k = gsl_rng_uniform_int(rng, 30);
		l = gsl_rng_uniform_int(rng, 30);

		x = h*asx + k*bsx + l*csx;
		y = h*asy + k*bsy + l*csy;
		z = h*asz + k*bsz + l*csz;

		nh = x*ax + y*ay + z*az;
		nk = x*bx + y*by + z*bz;
		nl = x*cx + y*cy + z*cz;

		dh = nh - lrint(nh);  dk = nk - lrint(nk);  dl = nl - lrint(nl);

		if ( fabs(dh) > 0.1 ) f++;
		if ( fabs(dk) > 0.1 ) f++;
		if ( fabs(dl) > 0.1 ) f++;

		ih = lrint(nh);  ik = lrint(nk);  il = lrint(nl);
		if ( forbidden_reflection(cell, ih, ik, il) ) {
			STATUS("Primitive %3i %3i %3i -> "
			       "Centered %3li %3li %3li, "
			       "which is forbidden\n",
			       h, k, l, ih, ik, il);
			fail = 1;
		}

		if ( f ) {
			STATUS("Primitive %3i %3i %3i -> "
			       "Centered %7.2f %7.2f %7.2f\n",
			       h, k, l, nh, nk, nl);
			fail = 1;
		}

	}

	return fail;
}


int main(int argc, char *argv[])
{
	int fail = 0;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	/* Triclinic P */
	fail += check_centering(50e-10, 55e-10, 70e-10, 67.0, 70.0, 77.0,
	                        L_TRICLINIC, 'P', '*', rng);

	/* Monoclinic P */
	fail += check_centering(10e-10, 20e-10, 30e-10, 100.0, 90.0, 90.0,
	                        L_MONOCLINIC, 'P', 'a', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 100.0, 90.0,
	                        L_MONOCLINIC, 'P', 'b', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 100.0,
	                        L_MONOCLINIC, 'P', 'c', rng);

	/* Monoclinic "C"-centered, unique axis a, three cell choices */
	fail += check_centering(10e-10, 20e-10, 30e-10, 100.0, 90.0, 90.0,
	                        L_MONOCLINIC, 'B', 'a', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 100.0, 90.0, 90.0,
	                        L_MONOCLINIC, 'C', 'a', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 100.0, 90.0, 90.0,
	                        L_MONOCLINIC, 'I', 'a', rng);

	/* Monoclinic "C"-centered, unique axis b, three cell choices */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 100.0, 90.0,
	                        L_MONOCLINIC, 'C', 'b', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 100.0, 90.0,
	                        L_MONOCLINIC, 'A', 'b', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 100.0, 90.0,
	                        L_MONOCLINIC, 'I', 'b', rng);

	/* Monoclinic "C"-centered, unique axis c, three cell choices */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 100.0,
	                        L_MONOCLINIC, 'A', 'c', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 100.0,
	                        L_MONOCLINIC, 'B', 'c', rng);
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 100.0,
	                        L_MONOCLINIC, 'I', 'c', rng);

	/* Orthorhombic P */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'P', '*', rng);

	/* Orthorhombic A */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'A', '*', rng);

	/* Orthorhombic B */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'B', '*', rng);

	/* Orthorhombic C */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'C', '*', rng);

	/* Orthorhombic I */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'I', '*', rng);

	/* Orthorhombic F */
	fail += check_centering(10e-10, 20e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_ORTHORHOMBIC, 'F', '*', rng);

	/* Tetragonal P */
	fail += check_centering(10e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'P', 'a', rng);
	fail += check_centering(30e-10, 10e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'P', 'b', rng);
	fail += check_centering(30e-10, 30e-10, 10e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'P', 'c', rng);

	/* Tetragonal I */
	fail += check_centering(10e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'I', 'a', rng);
	fail += check_centering(30e-10, 10e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'I', 'b', rng);
	fail += check_centering(30e-10, 30e-10, 10e-10, 90.0, 90.0, 90.0,
	                        L_TETRAGONAL, 'I', 'c', rng);

	/* Rhombohedral R */
	fail += check_centering(10e-10, 10e-10, 10e-10, 60.0, 60.0, 60.0,
	                        L_RHOMBOHEDRAL, 'R', '*', rng);

	/* Hexagonal P */
	fail += check_centering(30e-10, 10e-10, 10e-10, 120.0, 90.0, 90.0,
	                        L_HEXAGONAL, 'P', 'a', rng);
	fail += check_centering(10e-10, 30e-10, 10e-10, 90.0, 120.0, 90.0,
	                        L_HEXAGONAL, 'P', 'b', rng);
	fail += check_centering(10e-10, 10e-10, 30e-10, 90.0, 90.0, 120.0,
	                        L_HEXAGONAL, 'P', 'c', rng);

	/* Hexagonal H (PDB-speak for rhombohedral) */
	fail += check_centering(20e-10, 20e-10, 40e-10, 90.0, 90.0, 120.0,
	                        L_HEXAGONAL, 'H', 'c', rng);

	/* Cubic P */
	fail += check_centering(30e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_CUBIC, 'P', '*', rng);

	/* Cubic I */
	fail += check_centering(30e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_CUBIC, 'I', '*', rng);

	/* Cubic F */
	fail += check_centering(30e-10, 30e-10, 30e-10, 90.0, 90.0, 90.0,
	                        L_CUBIC, 'F', '*', rng);

	gsl_rng_free(rng);

	return fail;
}
