/*
 * cell-utils.c
 *
 * Unit Cell utility functions
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2021 Thomas White <taw@physics.org>
 *   2012      Lorenzo Galli
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

#include <libcrystfel-config.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "image.h"


/**
 * \file cell-utils.h
 **/


/* Weighting factor of lengths relative to angles */
#define LWEIGHT (10.0e-9)


/**
 * \param in: A UnitCell to rotate
 * \param quat: A quaternion
 *
 * Rotate a UnitCell using a quaternion.
 *
 * \returns a newly allocated rotated copy of \p in.
 *
 */
UnitCell *cell_rotate(UnitCell *in, struct quaternion quat)
{
	struct rvec a, b, c;
	struct rvec an, bn, cn;
	UnitCell *out = cell_new_from_cell(in);

	cell_get_cartesian(in, &a.u, &a.v, &a.w,
	                       &b.u, &b.v, &b.w,
	                       &c.u, &c.v, &c.w);

	an = quat_rot(a, quat);
	bn = quat_rot(b, quat);
	cn = quat_rot(c, quat);

	cell_set_cartesian(out, an.u, an.v, an.w,
	                        bn.u, bn.v, bn.w,
	                        cn.u, cn.v, cn.w);

	return out;
}


const char *str_lattice(LatticeType l)
{
	switch ( l )
	{
		case L_TRICLINIC :    return "triclinic";
		case L_MONOCLINIC :   return "monoclinic";
		case L_ORTHORHOMBIC : return "orthorhombic";
		case L_TETRAGONAL :   return "tetragonal";
		case L_RHOMBOHEDRAL : return "rhombohedral";
		case L_HEXAGONAL :    return "hexagonal";
		case L_CUBIC :        return "cubic";
	}

	return "unknown lattice";
}


LatticeType lattice_from_str(const char *s)
{
	if ( strcmp(s, "triclinic") == 0 ) return L_TRICLINIC;
	if ( strcmp(s, "monoclinic") == 0 ) return L_MONOCLINIC;
	if ( strcmp(s, "orthorhombic") == 0 ) return L_ORTHORHOMBIC;
	if ( strcmp(s, "tetragonal") == 0 ) return L_TETRAGONAL;
	if ( strcmp(s, "rhombohedral") == 0 ) return L_RHOMBOHEDRAL;
	if ( strcmp(s, "hexagonal") == 0 ) return L_HEXAGONAL;
	if ( strcmp(s, "cubic") == 0 ) return L_CUBIC;

	ERROR("Unrecognised lattice type '%s'\n", s);
	return L_TRICLINIC;
}


static int check_centering(char cen)
{
	switch ( cen ) {

		case 'P' :
		case 'A' :
		case 'B' :
		case 'C' :
		case 'I' :
		case 'F' :
		case 'R' :
		case 'H' :
		return 0;

		default:
		return 1;

	}
}


static int check_unique_axis(char ua)
{
	switch ( ua ) {

		case 'a' :
		case 'b' :
		case 'c' :
		return 0;

		default:
		return 1;

	}
}


int right_handed(UnitCell *cell)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	struct rvec aCb;
	double aCb_dot_c;
	int rh_reciprocal;
	int rh_direct;

	if ( cell_get_reciprocal(cell, &asx, &asy, &asz,
	                               &bsx, &bsy, &bsz,
	                               &csx, &csy, &csz) ) {
		ERROR("Couldn't get reciprocal cell.\n");
		return 0;
	}

	/* "a" cross "b" */
	aCb.u = asy*bsz - asz*bsy;
	aCb.v = - (asx*bsz - asz*bsx);
	aCb.w = asx*bsy - asy*bsx;

	/* "a cross b" dot "c" */
	aCb_dot_c = aCb.u*csx + aCb.v*csy + aCb.w*csz;

	rh_reciprocal = aCb_dot_c > 0.0;

	if ( cell_get_cartesian(cell, &asx, &asy, &asz,
	                              &bsx, &bsy, &bsz,
	                              &csx, &csy, &csz) ) {
		ERROR("Couldn't get direct cell.\n");
		return 0;
	}

	/* "a" cross "b" */
	aCb.u = asy*bsz - asz*bsy;
	aCb.v = - (asx*bsz - asz*bsx);
	aCb.w = asx*bsy - asy*bsx;

	/* "a cross b" dot "c" */
	aCb_dot_c = aCb.u*csx + aCb.v*csy + aCb.w*csz;

	rh_direct = aCb_dot_c > 0.0;

	if ( rh_reciprocal != rh_direct ) {
		ERROR("Whoops, reciprocal and real space handedness are "
		      "not the same!\n");
	}

	return rh_direct;
}


void cell_print(UnitCell *cell)
{
	LatticeType lt;
	char cen;

	if ( cell == NULL ) {
		STATUS("(NULL cell)\n");
		return;
	}

	lt = cell_get_lattice_type(cell);
	cen = cell_get_centering(cell);

	STATUS("%s %c", str_lattice(lt), cen);

	if ( (lt==L_MONOCLINIC) || (lt==L_TETRAGONAL) || ( lt==L_HEXAGONAL)
	  || ( (lt==L_ORTHORHOMBIC) && (cen=='A') )
	  || ( (lt==L_ORTHORHOMBIC) && (cen=='B') )
	  || ( (lt==L_ORTHORHOMBIC) && (cen=='C') ) )
	{
		STATUS(", unique axis %c", cell_get_unique_axis(cell));
	}

	if ( cell_has_parameters(cell) ) {
		if ( right_handed(cell) ) {
			STATUS(", right handed.\n");
		} else {
			STATUS(", left handed.\n");
		}
	} else {
		STATUS(".\n");
	}

	if ( cell_has_parameters(cell) ) {

		double a, b, c, alpha, beta, gamma;
		cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

		STATUS("a      b      c            alpha   beta  gamma\n");
		STATUS("%6.2f %6.2f %6.2f A    %6.2f %6.2f %6.2f deg\n",
		       a*1e10, b*1e10, c*1e10,
		       rad2deg(alpha), rad2deg(beta), rad2deg(gamma));
	} else {
		STATUS("Unit cell parameters are not specified.\n");
	}
}


void cell_print_oneline(UnitCell *cell)
{
	LatticeType lt;
	char cen;

	if ( cell == NULL ) {
		STATUS("(NULL cell)\n");
		return;
	}

	lt = cell_get_lattice_type(cell);
	cen = cell_get_centering(cell);

	STATUS("%s %c", str_lattice(lt), cen);

	if ( (lt==L_MONOCLINIC) || (lt==L_TETRAGONAL) || ( lt==L_HEXAGONAL)
	  || ( (lt==L_ORTHORHOMBIC) && (cen=='A') )
	  || ( (lt==L_ORTHORHOMBIC) && (cen=='B') )
	  || ( (lt==L_ORTHORHOMBIC) && (cen=='C') ) )
	{
		STATUS(", unique axis %c", cell_get_unique_axis(cell));
	}

	if ( cell_has_parameters(cell) ) {

		double a, b, c, alpha, beta, gamma;

		if ( !right_handed(cell) ) {
			STATUS(" (left handed)");
		}

		cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

		STATUS("  %.2f  %.2f  %.2f A,   %.2f  %.2f  %.2f deg\n",
		       a*1e10, b*1e10, c*1e10,
		       rad2deg(alpha), rad2deg(beta), rad2deg(gamma));
	} else {
		STATUS(", no cell parameters.\n");
	}
}


void cell_print_full(UnitCell *cell)
{
	cell_print(cell);
	if ( cell == NULL ) return;

	if ( cell_has_parameters(cell) ) {

		double asx, asy, asz;
		double bsx, bsy, bsz;
		double csx, csy, csz;
		double ax, ay, az, bx, by, bz, cx, cy, cz;

		cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

		STATUS("a = %10.3e %10.3e %10.3e m\n", ax, ay, az);
		STATUS("b = %10.3e %10.3e %10.3e m\n", bx, by, bz);
		STATUS("c = %10.3e %10.3e %10.3e m\n", cx, cy, cz);

		cell_get_reciprocal(cell, &asx, &asy, &asz,
			                  &bsx, &bsy, &bsz,
			                  &csx, &csy, &csz);

		STATUS("a* = %10.3e %10.3e %10.3e m^-1 (modulus %10.3e m^-1)\n",
		       asx, asy, asz, modulus(asx, asy, asz));
		STATUS("b* = %10.3e %10.3e %10.3e m^-1 (modulus %10.3e m^-1)\n",
		       bsx, bsy, bsz, modulus(bsx, bsy, bsz));
		STATUS("c* = %10.3e %10.3e %10.3e m^-1 (modulus %10.3e m^-1)\n",
		       csx, csy, csz, modulus(csx, csy, csz));

		STATUS("alpha* = %6.2f deg, beta* = %6.2f deg, "
		       "gamma* = %6.2f deg\n",
		       rad2deg(angle_between(bsx, bsy, bsz, csx, csy, csz)),
		       rad2deg(angle_between(asx, asy, asz, csx, csy, csz)),
		       rad2deg(angle_between(asx, asy, asz, bsx, bsy, bsz)));

	}
}


int bravais_lattice(UnitCell *cell)
{
	LatticeType lattice = cell_get_lattice_type(cell);
	char centering = cell_get_centering(cell);
	char ua = cell_get_unique_axis(cell);

	switch ( centering )
	{
		case 'P' :
		return 1;

		case 'A' :
		case 'B' :
		case 'C' :
		if ( lattice == L_MONOCLINIC ) {
			if ( (ua=='a') && (centering!='A') ) return 1;
			if ( (ua=='b') && (centering!='B') ) return 1;
			if ( (ua=='c') && (centering!='C') ) return 1;
		} else if ( lattice == L_ORTHORHOMBIC) {
			return 1;
		}
		return 0;

		case 'I' :
		/* We accept monoclinic I as "Bravais", even though it's
		 * unconventional */
		if ( (lattice == L_MONOCLINIC)
		  || (lattice == L_ORTHORHOMBIC)
		  || (lattice == L_TETRAGONAL)
		  || (lattice == L_CUBIC) )
		{
			return 1;
		}
		return 0;

		case 'F' :
		if ( (lattice == L_ORTHORHOMBIC) || (lattice == L_CUBIC) ) {
			return 1;
		}
		return 0;

		case 'H' :
		/* "Hexagonal H" is not a Bravais lattice, but rather something
		 * invented by the PDB to make life difficult for programmers.
		 * Accepting it as Bravais seems to be the least painful way to
		 * handle it correctly. Yuk. */
		if ( ua != 'c' ) return 0;
		if ( lattice == L_HEXAGONAL ) return 1;
		return 0;

		case 'R' :
		if ( lattice == L_RHOMBOHEDRAL ) return 1;
		return 0;

		default :
		return 0;
	}
}


static RationalMatrix *create_rtnl_mtx(signed int a1, signed int a2,
                                       signed int b1, signed int b2,
                                       signed int c1, signed int c2,
                                       signed int d1, signed int d2,
                                       signed int e1, signed int e2,
                                       signed int f1, signed int f2,
                                       signed int g1, signed int g2,
                                       signed int h1, signed int h2,
                                       signed int i1, signed int i2)
{
	RationalMatrix *m = rtnl_mtx_new(3, 3);
	if ( m == NULL ) return NULL;
	rtnl_mtx_set(m, 0, 0, rtnl(a1, a2));
	rtnl_mtx_set(m, 0, 1, rtnl(b1, b2));
	rtnl_mtx_set(m, 0, 2, rtnl(c1, c2));
	rtnl_mtx_set(m, 1, 0, rtnl(d1, d2));
	rtnl_mtx_set(m, 1, 1, rtnl(e1, e2));
	rtnl_mtx_set(m, 1, 2, rtnl(f1, f2));
	rtnl_mtx_set(m, 2, 0, rtnl(g1, g2));
	rtnl_mtx_set(m, 2, 1, rtnl(h1, h2));
	rtnl_mtx_set(m, 2, 2, rtnl(i1, i2));
	return m;
}


/* Given a centered cell "in", return the integer transformation matrix which
 * turns a primitive cell into "in". Set new_centering and new_latt to the
 * centering and lattice type of the primitive cell (usually aP, sometimes rR,
 * rarely mP).  Store the inverse matrix at pCi */
static IntegerMatrix *centering_transformation(UnitCell *in,
                                               char *new_centering,
                                               LatticeType *new_latt,
                                               char *new_ua,
                                               RationalMatrix **pCi)
{
	LatticeType lt;
	char ua, cen;
	IntegerMatrix *C = NULL;
	RationalMatrix *Ci = NULL;

	lt = cell_get_lattice_type(in);
	ua = cell_get_unique_axis(in);
	cen = cell_get_centering(in);

	/* Write the matrices exactly as they appear in ITA Table 5.1.3.1.
	 * C is "P", and Ci is "Q=P^-1".  Vice-versa if the transformation
	 * should go the opposite way to what's written in the first column. */

	if ( (cen=='P') || (cen=='R') ) {
		*new_centering = 'P';
		*new_latt = lt;
		*new_ua = ua;
		C = intmat_identity(3);
		Ci = rtnl_mtx_identity(3);
	}

	if ( cen == 'I' ) {
		C = intmat_create_3x3(0, 1, 1,
		                      1, 0, 1,
		                      1, 1, 0);
		Ci = create_rtnl_mtx(-1,2,  1,2,  1,2,
		                      1,2, -1,2,  1,2,
		                      1,2,  1,2, -1,2);
		if ( lt == L_CUBIC ) {
			*new_latt = L_RHOMBOHEDRAL;
			*new_centering = 'R';
			*new_ua = '*';
		} else {
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
			*new_ua = '*';
		}
	}

	if ( cen == 'F' ) {
		C = intmat_create_3x3(-1,  1,  1,
		                       1, -1,  1,
		                       1,  1, -1);
		Ci = create_rtnl_mtx( 0,1,  1,2,  1,2,
		                      1,2,  0,1,  1,2,
		                      1,2,  1,2,  0,1);
		if ( lt == L_CUBIC ) {
			*new_latt = L_RHOMBOHEDRAL;
			*new_centering = 'R';
			*new_ua = '*';
		} else {
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
			*new_ua = '*';
		}
	}

	if ( (lt == L_HEXAGONAL) && (cen == 'H') && (ua == 'c') ) {
		/* Obverse setting */
		C = intmat_create_3x3( 1,  0,  1,
		                      -1,  1,  1,
		                       0, -1,  1);
		Ci = create_rtnl_mtx( 2,3, -1,3, -1,3,
		                      1,3,  1,3, -2,3,
		                      1,3,  1,3,  1,3);
		assert(lt == L_HEXAGONAL);
		assert(ua == 'c');
		*new_latt = L_RHOMBOHEDRAL;
		*new_centering = 'R';
		*new_ua = '*';
	}

	if ( cen == 'A' ) {
		C = intmat_create_3x3( 1,  0,  0,
		                       0,  1,  1,
		                       0, -1,  1);
		Ci = create_rtnl_mtx( 1,1,  0,1,  0,1,
		                      0,1,  1,2, -1,2,
		                      0,1,  1,2,  1,2);
		if ( lt == L_ORTHORHOMBIC ) {
			*new_latt = L_MONOCLINIC;
			*new_centering = 'P';
			*new_ua = 'a';
		} else {
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
			*new_ua = '*';
		}
	}

	if ( cen == 'B' ) {
		C = intmat_create_3x3( 1,  0,  1,
		                       0,  1,  0,
		                      -1,  0,  1);
		Ci = create_rtnl_mtx( 1,2,  0,1, -1,2,
		                      0,1,  1,1,  0,1,
		                      1,2,  0,1,  1,2);
		if ( lt == L_ORTHORHOMBIC ) {
			*new_latt = L_MONOCLINIC;
			*new_centering = 'P';
			*new_ua = 'b';
		} else {
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
			*new_ua = '*';
		}
	}

	if ( cen == 'C' ) {
		C = intmat_create_3x3( 1,  1,  0,
		                      -1,  1,  0,
		                       0,  0,  1);
		Ci = create_rtnl_mtx( 1,2, -1,2,  0,1,
		                      1,2,  1,2,  0,1,
		                      0,1,  0,1,  1,1);
		if ( lt == L_ORTHORHOMBIC ) {
			*new_latt = L_MONOCLINIC;
			*new_centering = 'P';
			*new_ua = 'c';
		} else {
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
			*new_ua = '*';
		}
	}

	*pCi = Ci;
	return C;
}


/**
 * \param in: A %UnitCell
 * \param pC: Location at which to store the centering transformation
 * \param pCi: Location at which to store the inverse centering transformation
 *
 * Turns any cell into a primitive one, e.g. for comparison purposes.
 *
 * The transformation which was used is stored at \p Ci. The centering
 * transformation, which is the transformation you should apply if you want to
 * get back the original cell, will be stored at \p C.  Either or both of these
 * can be NULL if you don't need that information.
 *
 * \returns a primitive version of \p in.
 *
 */
UnitCell *uncenter_cell(UnitCell *in, IntegerMatrix **pC, RationalMatrix **pCi)
{
	IntegerMatrix *C;
	RationalMatrix *Ci;
	char new_centering;
	LatticeType new_latt;
	char new_ua;
	UnitCell *out;

	C = centering_transformation(in, &new_centering, &new_latt,
	                             &new_ua, &Ci);
	if ( C == NULL ) return NULL;

	out = cell_transform_rational(in, Ci);
	if ( out == NULL ) return NULL;

	cell_set_lattice_type(out, new_latt);
	cell_set_centering(out, new_centering);
	cell_set_unique_axis(out, new_ua);

	if ( pC != NULL ) {
		*pC = C;
	} else {
		intmat_free(C);
	}

	if ( pCi != NULL ) {
		*pCi = Ci;
	} else {
		rtnl_mtx_free(Ci);
	}

	return out;
}


/* Return sin(theta)/lambda = 1/2d.  Multiply by two if you want 1/d */
double resolution(UnitCell *cell, signed int h, signed int k, signed int l)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell_get_reciprocal(cell,
	                    &asx, &asy, &asz,
	                    &bsx, &bsy, &bsz,
	                    &csx, &csy, &csz);

	return modulus(h*asx + k*bsx + l*csx,
	               h*asy + k*bsy + l*csy,
	               h*asz + k*bsz + l*csz) / 2.0;
}


static void determine_lattice(UnitCell *cell,
                              const char *as, const char *bs, const char *cs,
                              const char *als, const char *bes, const char *gas)
{
	int n_right;

	/* Rhombohedral or cubic? */
	if ( (strcmp(as, bs) == 0) && (strcmp(as, cs) == 0) ) {

		if ( (strcmp(als, "  90.00") == 0)
		  && (strcmp(bes, "  90.00") == 0)
		  && (strcmp(gas, "  90.00") == 0) )
		{
			/* Cubic.  Unique axis irrelevant. */
			cell_set_lattice_type(cell, L_CUBIC);
			return;
		}

		if ( (strcmp(als, bes) == 0) && (strcmp(als, gas) == 0) ) {
			/* Rhombohedral.  Unique axis irrelevant. */
			cell_set_lattice_type(cell, L_RHOMBOHEDRAL);
			return;
		}

	}

	if ( (strcmp(als, "  90.00") == 0)
	  && (strcmp(bes, "  90.00") == 0)
	  && (strcmp(gas, "  90.00") == 0) )
	{
		if ( strcmp(bs, cs) == 0 ) {
			/* Tetragonal, unique axis a */
			cell_set_lattice_type(cell, L_TETRAGONAL);
			cell_set_unique_axis(cell, 'a');
			return;
		}

		if ( strcmp(as, cs) == 0 ) {
			/* Tetragonal, unique axis b */
			cell_set_lattice_type(cell, L_TETRAGONAL);
			cell_set_unique_axis(cell, 'b');
			return;
		}

		if ( strcmp(as, bs) == 0 ) {
			/* Tetragonal, unique axis c */
			cell_set_lattice_type(cell, L_TETRAGONAL);
			cell_set_unique_axis(cell, 'c');
			return;
		}

		/* Orthorhombic.  Unique axis irrelevant, but point group
		 * can have different orientations. */
		cell_set_lattice_type(cell, L_ORTHORHOMBIC);
		cell_set_unique_axis(cell, '*');
		return;
	}

	n_right = 0;
	if ( strcmp(als, "  90.00") == 0 ) n_right++;
	if ( strcmp(bes, "  90.00") == 0 ) n_right++;
	if ( strcmp(gas, "  90.00") == 0 ) n_right++;

	/* Hexgonal or monoclinic? */
	if ( n_right == 2 ) {

		if ( (strcmp(als, " 120.00") == 0)
		  && (strcmp(bs, cs) == 0) )
		{
			/* Hexagonal, unique axis a */
			cell_set_lattice_type(cell, L_HEXAGONAL);
			cell_set_unique_axis(cell, 'a');
			return;
		}

		if ( (strcmp(bes, " 120.00") == 0)
		  && (strcmp(as, cs) == 0) )
		{
			/* Hexagonal, unique axis b */
			cell_set_lattice_type(cell, L_HEXAGONAL);
			cell_set_unique_axis(cell, 'b');
			return;
		}

		if ( (strcmp(gas, " 120.00") == 0)
		  && (strcmp(as, bs) == 0) )
		{
			/* Hexagonal, unique axis c */
			cell_set_lattice_type(cell, L_HEXAGONAL);
			cell_set_unique_axis(cell, 'c');
			return;
		}

		if ( strcmp(als, "  90.00") != 0 ) {
			/* Monoclinic, unique axis a */
			cell_set_lattice_type(cell, L_MONOCLINIC);
			cell_set_unique_axis(cell, 'a');
			return;
		}

		if ( strcmp(bes, "  90.00") != 0 ) {
			/* Monoclinic, unique axis b */
			cell_set_lattice_type(cell, L_MONOCLINIC);
			cell_set_unique_axis(cell, 'b');
			return;
		}

		if ( strcmp(gas, "  90.00") != 0 ) {
			/* Monoclinic, unique axis c */
			cell_set_lattice_type(cell, L_MONOCLINIC);
			cell_set_unique_axis(cell, 'c');
			return;
		}
	}

	/* Triclinic, unique axis irrelevant. */
	cell_set_lattice_type(cell, L_TRICLINIC);
}


/**
 * \param filename: The filename from which to load the cell
 *
 * Loads a unit cell from the CRYST1 line of a PDB file.
 *
 * \returns a newly allocated %UnitCell.
 *
 */
UnitCell *load_cell_from_pdb(const char *filename)
{
	FILE *fh;
	char *rval;
	UnitCell *cell = NULL;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		return NULL;
	}

	do {

		char line[1024];

		rval = fgets(line, 1023, fh);

		if ( strncmp(line, "CRYST1", 6) == 0 ) {

			float a, b, c, al, be, ga;
			char as[10], bs[10], cs[10];
			char als[8], bes[8], gas[8];
			int r;

			memcpy(as, line+6, 9);    as[9] = '\0';
			memcpy(bs, line+15, 9);   bs[9] = '\0';
			memcpy(cs, line+24, 9);   cs[9] = '\0';
			memcpy(als, line+33, 7);  als[7] = '\0';
			memcpy(bes, line+40, 7);  bes[7] = '\0';
			memcpy(gas, line+47, 7);  gas[7] = '\0';

			r = sscanf(as, "%f", &a);
			r += sscanf(bs, "%f", &b);
			r += sscanf(cs, "%f", &c);
			r += sscanf(als, "%f", &al);
			r += sscanf(bes, "%f", &be);
			r += sscanf(gas, "%f", &ga);

			if ( r != 6 ) {
				STATUS("Couldn't understand CRYST1 line.\n");
				continue;
			}

			cell = cell_new_from_parameters(a*1e-10,
			                                b*1e-10, c*1e-10,
	                                                deg2rad(al),
	                                                deg2rad(be),
	                                                deg2rad(ga));

			determine_lattice(cell, as, bs, cs, als, bes, gas);

			if ( strlen(line) > 55 ) {
				cell_set_centering(cell, line[55]);
			} else {
				ERROR("CRYST1 line without centering.\n");
			}

			break;  /* Done */
		}

	} while ( rval != NULL );

	fclose(fh);

	if ( cell != NULL ) {
		validate_cell(cell);
	} else {
		ERROR("Failed to load cell from %s\n", filename);
	}


	return cell;
}


static int get_length_m(char **bits, int nbits, double *pl)
{
	char *rval;

	if ( nbits < 4 ) {
		ERROR("No units specified for '%s'\n", bits[0]);
		return 1;
	}

	*pl = strtod(bits[2], &rval);
	if ( *rval != '\0' ) {
		ERROR("Invalid value '%s'.\n", bits[2]);
		return 1;
	}

	if ( strcmp(bits[3], "nm") == 0 ) {
		*pl *= 1e-9;
	} else if ( strcmp(bits[3], "A") == 0 ) {
		*pl *= 1e-10;
	} else {
		ERROR("Unrecognised length units '%s'\n", bits[3]);
		return 1;
	}

	return 0;
}


static int get_angle_rad(char **bits, int nbits, double *pl)
{
	char *rval;

	if ( nbits < 4 ) {
		ERROR("No units specified for '%s'\n", bits[0]);
		return 1;
	}

	*pl = strtod(bits[2], &rval);
	if ( *rval != '\0' ) {
		ERROR("Invalid value '%s'.\n", bits[2]);
		return 1;
	}

	if ( strcmp(bits[3], "rad") == 0 ) {
		/* Do nothing, already in rad */
	} else if ( strcmp(bits[3], "deg") == 0 ) {
		*pl = deg2rad(*pl);
	} else {
		ERROR("Unrecognised angle units '%s'\n", bits[3]);
		return 1;
	}

	return 0;
}

/**
 * \param cell: a %UnitCell
 * \param fh: a file handle
 *
 * Writes \p cell to \p fh, in CrystFEL unit cell file format
 *
 */
void write_cell(UnitCell *cell, FILE *fh)
{
	double a, b, c, al, be, ga;
	LatticeType lt;

	fprintf(fh, "CrystFEL unit cell file version 1.0\n\n");
	lt = cell_get_lattice_type(cell);
	fprintf(fh, "lattice_type = %s\n", str_lattice(lt));
	if ( (lt == L_MONOCLINIC)
	  || (lt == L_TETRAGONAL)
	  || (lt == L_HEXAGONAL) )
	{
		fprintf(fh, "unique_axis = %c\n", cell_get_unique_axis(cell));
	}
	fprintf(fh, "centering = %c\n", cell_get_centering(cell));

	if ( cell_has_parameters(cell) ) {
		cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
		fprintf(fh, "a = %.2f A\n", a*1e10);
		fprintf(fh, "b = %.2f A\n", b*1e10);
		fprintf(fh, "c = %.2f A\n", c*1e10);
		fprintf(fh, "al = %.2f deg\n", rad2deg(al));
		fprintf(fh, "be = %.2f deg\n", rad2deg(be));
		fprintf(fh, "ga = %.2f deg\n", rad2deg(ga));
	}
}


/**
 * \param filename: The filename from which to load the cell
 *
 * Loads a unit cell from a file of any type (PDB or CrystFEL format)
 *
 * \returns a newly allocated %UnitCell.
 *
 */
UnitCell *load_cell_from_file(const char *filename)
{
	FILE *fh;
	char *rval;
	char line[1024];
	UnitCell *cell;
	int have_a = 0, have_b = 0, have_c = 0;
	int have_al = 0,  have_be = 0,  have_ga = 0;
	double a, b, c, al, be, ga;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Couldn't open '%s'\n", filename);
		return NULL;
	}

	rval = fgets(line, 1023, fh);
	chomp(line);

	if ( strcmp(line, "CrystFEL unit cell file version 1.0") != 0 ) {
		fclose(fh);
		return load_cell_from_pdb(filename);
	}

	cell = cell_new();

	do {

		char line[1024];
		int n1;
		int i;
		char **bits;

		rval = fgets(line, 1023, fh);
		chomp(line);

		n1 = assplode(line, " \t", &bits, ASSPLODE_NONE);
		if ( n1 < 3 ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		if ( bits[0][0] == ';' ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		if ( strcmp(bits[0], "lattice_type") == 0 ) {
			LatticeType lt = lattice_from_str(bits[2]);
			cell_set_lattice_type(cell, lt);

		} else if ( strcmp(bits[0], "centering") == 0 ) {
			char cen = bits[2][0];
			if ( !check_centering(cen) ) {
				cell_set_centering(cell, cen);
			} else {
				ERROR("Unrecognised centering '%c'\n", cen);
			}

		} else if ( strcmp(bits[0], "unique_axis") == 0 ) {
			char ua = bits[2][0];
			if ( !check_unique_axis(ua) ) {
				cell_set_unique_axis(cell, ua);
			} else {
				ERROR("Unrecognised unique axis '%c'\n", ua);
			}

		} else if ( strcmp(bits[0], "a") == 0 ) {
			if ( !get_length_m(bits, n1, &a) ) {
				have_a = 1;
			}

		} else if ( strcmp(bits[0], "b") == 0 ) {
			if ( !get_length_m(bits, n1, &b) ) {
				have_b = 1;
			}

		} else if ( strcmp(bits[0], "c") == 0 ) {
			if ( !get_length_m(bits, n1, &c) ) {
				have_c = 1;
			}

		} else if ( strcmp(bits[0], "al") == 0 ) {
			if ( !get_angle_rad(bits, n1, &al) ) {
				have_al = 1;
			}

		} else if ( strcmp(bits[0], "be") == 0 ) {
			if ( !get_angle_rad(bits, n1, &be) ) {
				have_be = 1;
			}

		} else if ( strcmp(bits[0], "ga") == 0 ) {
			if ( !get_angle_rad(bits, n1, &ga) ) {
				have_ga = 1;
			}

		} else {
			ERROR("Unrecognised field '%s'\n", bits[0]);
		}

		for ( i=0; i<n1; i++ ) free(bits[i]);
		free(bits);

	} while ( rval != NULL );

	fclose(fh);

	if ( have_a && have_b && have_c && have_al && have_be && have_ga ) {
		cell_set_parameters(cell, a, b, c, al, be, ga);
	}

	switch ( cell_get_lattice_type(cell) ) {

		case L_TRICLINIC :
		case L_ORTHORHOMBIC :
		case L_CUBIC :
		case L_RHOMBOHEDRAL :
		if ( (cell_get_unique_axis(cell) != '?')
		  && (cell_get_unique_axis(cell) != '*') ) {
			ERROR("WARNING: Unique axis '%c' doesn't make sense "
			      "for lattice type %s.\n",
			      cell_get_unique_axis(cell),
			      str_lattice(cell_get_lattice_type(cell)));
		}
		break;

		case L_MONOCLINIC :
		case L_TETRAGONAL :
		case L_HEXAGONAL :
		if ( (cell_get_unique_axis(cell) == '?')
		  || (cell_get_unique_axis(cell) == '*') ) {
			ERROR("You must specify the unique axis for lattice "
			      "type %s.\n",
			      str_lattice(cell_get_lattice_type(cell)));
			return NULL;
		}
		break;

		default :
		ERROR("Unrecognised lattice type %i\n",
		      cell_get_lattice_type(cell));
		break;
	}

	validate_cell(cell);

	return cell;
}


/* Force the linker to bring in CBLAS to make GSL happy */
void cell_fudge_gslcblas()
{
        STATUS("%p\n", cblas_sgemm);
}


/**
 * \param in: A %UnitCell to rotate
 * \param omega: Euler angle about +z
 * \param phi: Euler angle about +x
 * \param rot: Euler angle about new +z
 *
 * Rotate a %UnitCell using Euler angles
 *
 * \returns a newly allocated rotated copy of \p in.
 *
 */
UnitCell *rotate_cell(UnitCell *in, double omega, double phi, double rot)
{
	UnitCell *out;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double xnew, ynew, znew;

	cell_get_reciprocal(in, &asx, &asy, &asz, &bsx, &bsy,
	                        &bsz, &csx, &csy, &csz);

	/* Rotate by "omega" about +z (parallel to c* and c unless triclinic) */
	xnew = asx*cos(omega) + asy*sin(omega);
	ynew = -asx*sin(omega) + asy*cos(omega);
	znew = asz;
	asx = xnew;  asy = ynew;  asz = znew;
	xnew = bsx*cos(omega) + bsy*sin(omega);
	ynew = -bsx*sin(omega) + bsy*cos(omega);
	znew = bsz;
	bsx = xnew;  bsy = ynew;  bsz = znew;
	xnew = csx*cos(omega) + csy*sin(omega);
	ynew = -csx*sin(omega) + csy*cos(omega);
	znew = csz;
	csx = xnew;  csy = ynew;  csz = znew;

	/* Rotate by "phi" about +x (not parallel to anything specific) */
	xnew = asx;
	ynew = asy*cos(phi) + asz*sin(phi);
	znew = -asy*sin(phi) + asz*cos(phi);
	asx = xnew;  asy = ynew;  asz = znew;
	xnew = bsx;
	ynew = bsy*cos(phi) + bsz*sin(phi);
	znew = -bsy*sin(phi) + bsz*cos(phi);
	bsx = xnew;  bsy = ynew;  bsz = znew;
	xnew = csx;
	ynew = csy*cos(phi) + csz*sin(phi);
	znew = -csy*sin(phi) + csz*cos(phi);
	csx = xnew;  csy = ynew;  csz = znew;

	/* Rotate by "rot" about the new +z (in-plane rotation) */
	xnew = asx*cos(rot) + asy*sin(rot);
	ynew = -asx*sin(rot) + asy*cos(rot);
	znew = asz;
	asx = xnew;  asy = ynew;  asz = znew;
	xnew = bsx*cos(rot) + bsy*sin(rot);
	ynew = -bsx*sin(rot) + bsy*cos(rot);
	znew = bsz;
	bsx = xnew;  bsy = ynew;  bsz = znew;
	xnew = csx*cos(rot) + csy*sin(rot);
	ynew = -csx*sin(rot) + csy*cos(rot);
	znew = csz;
	csx = xnew;  csy = ynew;  csz = znew;

	out = cell_new_from_cell(in);
	cell_set_reciprocal(out, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	return out;
}


int cell_is_sensible(UnitCell *cell)
{
	double a, b, c, al, be, ga;

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	if (   al + be + ga >= 2.0*M_PI ) return 0;
	if (   al + be - ga >= 2.0*M_PI ) return 0;
	if (   al - be + ga >= 2.0*M_PI ) return 0;
	if ( - al + be + ga >= 2.0*M_PI ) return 0;
	if (   al + be + ga <= 0.0 ) return 0;
	if (   al + be - ga <= 0.0 ) return 0;
	if (   al - be + ga <= 0.0 ) return 0;
	if ( - al + be + ga <= 0.0 ) return 0;
	if ( isnan(al) ) return 0;
	if ( isnan(be) ) return 0;
	if ( isnan(ga) ) return 0;
	return 1;
}


/**
 * \param cell: A %UnitCell to validate
 *
 * Perform some checks for crystallographic validity \p cell, such as that the
 * lattice is a conventional Bravais lattice.
 * Warnings are printied if any of the checks are failed.
 *
 * \returns zero if the cell is fine, 1 if it is unconventional but otherwise
 *  OK (e.g. left-handed or not a Bravais lattice), and 2 if there is a serious
 *  problem such as the parameters being physically impossible.
 *
 */
int validate_cell(UnitCell *cell)
{
	int err = 0;
	char cen, ua;

	if ( cell_has_parameters(cell) && !cell_is_sensible(cell) ) {
		ERROR("WARNING: Unit cell parameters are not sensible.\n");
		err = 2;
	}

	if ( !bravais_lattice(cell) ) {
		ERROR("WARNING: Unit cell is not a conventional Bravais"
		      " lattice.\n");
		err = 1;
	}

	if ( cell_has_parameters(cell) && !right_handed(cell) ) {
		ERROR("WARNING: Unit cell is not right handed.\n");
		err = 1;
	}

	/* For monoclinic A, B or C centering, the unique axis must be something
	 * other than the centering. */
	if ( cell_get_lattice_type(cell) == L_MONOCLINIC ) {
		cen = cell_get_centering(cell);
		ua = cell_get_unique_axis(cell);
		if ( ((cen == 'A') && (ua == 'a'))
		  || ((cen == 'B') && (ua == 'b'))
		  || ((cen == 'C') && (ua == 'c')) ) {
			ERROR("WARNING: A, B or C centering matches unique"
			      " axis.\n");
			err = 2;
		}
	}

	return err;
}


/**
 * \param cell: A %UnitCell
 * \param h: h index to check
 * \param k: k index to check
 * \param l: l index to check
 *
 * \returns true if this reflection is forbidden.
 *
 */
int forbidden_reflection(UnitCell *cell,
                         signed int h, signed int k, signed int l)
{
	char cen;

	cen = cell_get_centering(cell);

	/* Reflection conditions here must match the transformation matrices
	 * in centering_transformation().  tests/centering_check verifies
	 * this (amongst other things). */

	if ( cen == 'P' ) return 0;
	if ( cen == 'R' ) return 0;

	if ( cen == 'A' ) return (k+l) % 2;
	if ( cen == 'B' ) return (h+l) % 2;
	if ( cen == 'C' ) return (h+k) % 2;

	if ( cen == 'I' ) return (h+k+l) % 2;
	if ( cen == 'F' ) return ((h+k) % 2) || ((h+l) % 2) || ((k+l) % 2);

	/* Obverse setting */
	if ( cen == 'H' ) return (-h+k+l) % 3;

	return 0;
}


/**
 * \returns cell volume in m^3
 */
double cell_get_volume(UnitCell *cell)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	struct rvec aCb;

	if ( cell_get_cartesian(cell, &ax, &ay, &az,
	                              &bx, &by, &bz,
	                              &cx, &cy, &cz) ) {
		ERROR("Couldn't get reciprocal cell.\n");
		return 0;
	}

	/* "a" cross "b" */
	aCb.u = ay*bz - az*by;
	aCb.v = - (ax*bz - az*bx);
	aCb.w = ax*by - ay*bx;

	/* "a cross b" dot "c" */
	return (aCb.u*cx + aCb.v*cy + aCb.w*cz);
}


/**
 * \param cell: A %UnitCell
 *
 * \returns the value of 1/d for the lowest order reflection
 * that is not systematically absent according to the centering.
 *
 */
double lowest_reflection(UnitCell *cell)
{
	signed int h, k, l;
	double lowres = INFINITY;

	/* FIXME: Inelegant and nasty.  Anyone want to work out
	 * all the possible cases? */
	for ( h=0; h<4; h++ ) {
		for ( k=0; k<4; k++ ) {
			for ( l=0; l<4; l++ ) {
				if ( (h==0) && (k==0) && (l==0) ) continue;
				if ( !forbidden_reflection(cell, h, k, l) ) {
					double r = resolution(cell, h, k, l);
					if ( r < lowres ) {
						lowres = r;
					}
				}
			}
		}
	}
	return lowres;
}


/* Return true if the two centering symbols are identical,
 * or if they are a pair of R/P, which should be considered the
 * same for the purposes of cell comparison */
static int centering_equivalent(char cen1, char cen2)
{
	if ( cen1 == cen2 ) return 1;
	if ( (cen1=='P') && (cen2=='R') ) return 1;
	if ( (cen1=='R') && (cen2=='P') ) return 1;
	return 0;
}


/**
 * \param cell: A UnitCell
 * \param reference: Another UnitCell
 * \param tols: Pointer to tolerances for a,b,c (fractional), al,be,ga (radians)
 *
 * Compare the two unit cells.  If the real space parameters match to within
 * the specified tolerances, and the centering matches, this function returns 1.
 * Otherwise 0.
 *
 * This function considers the cell parameters and centering, but ignores the
 * orientation of the cell.  If you want to compare the orientation as well,
 * use compare_cell_parameters_and_orientation() instead.
 *
 * \returns non-zero if the cells match.
 *
 */
int compare_cell_parameters(UnitCell *cell, UnitCell *reference,
                            const double *tols)
{
	double a1, b1, c1, al1, be1, ga1;
	double a2, b2, c2, al2, be2, ga2;

	/* Centering must match: we don't arbitrate primitive vs centered,
	 * different cell choices etc */
	if ( !centering_equivalent(cell_get_centering(cell),
	                           cell_get_centering(reference)) ) return 0;

	cell_get_parameters(cell, &a1, &b1, &c1, &al1, &be1, &ga1);
	cell_get_parameters(reference, &a2, &b2, &c2, &al2, &be2, &ga2);

	/* within_tolerance() takes a percentage */
	if ( !within_tolerance(a1, a2, tols[0]*100.0) ) return 0;
	if ( !within_tolerance(b1, b2, tols[1]*100.0) ) return 0;
	if ( !within_tolerance(c1, c2, tols[2]*100.0) ) return 0;
	if ( fabs(al1-al2) > tols[3] ) return 0;
	if ( fabs(be1-be2) > tols[4] ) return 0;
	if ( fabs(ga1-ga2) > tols[5] ) return 0;

	return 1;
}


static double moduli_check(double ax, double ay, double az,
                           double bx, double by, double bz)
{
	double ma = modulus(ax, ay, az);
	double mb = modulus(bx, by, bz);
	return fabs(ma-mb)/ma;
}


/**
 * \param cell: A UnitCell
 * \param reference: Another UnitCell
 * \param tols: Pointer to six tolerance values (see below)
 *
 * Compare the two unit cells.  If the axes match in length (to within
 * the specified tolerances), this function returns non-zero.
 *
 * This function compares the orientation of the cell as well as the parameters.
 * If you just want to see if the parameters are the same, use
 * compare_cell_parameters() instead.
 *
 * The comparison is done by checking that the lengths of the unit cell axes are
 * the same between the two cells, and that the axes have similar directions in
 * 3D space.  The first three tolerance values are the maximum allowable fractional
 * differences between the a,b,c axis lengths (respectively) of the two cells.
 * The last three tolerance values are the maximum allowable angles, in radians,
 * between the directions of the a,b,c axes of the two cells.
 *
 * \p cell and \p reference must have the same centering.  Otherwise, this
 * function always returns zero.
 *
 * \returns non-zero if the cells match.
 *
 */
/* 'tols' is in frac (not %) and radians */
int compare_cell_parameters_and_orientation(UnitCell *cell, UnitCell *reference,
                                            const double *tols)
{
	double ax1, ay1, az1, bx1, by1, bz1, cx1, cy1, cz1;
	double ax2, ay2, az2, bx2, by2, bz2, cx2, cy2, cz2;

	if ( cell_get_centering(cell) != cell_get_centering(reference) ) return 0;

	cell_get_cartesian(cell, &ax1, &ay1, &az1,
	                         &bx1, &by1, &bz1,
	                         &cx1, &cy1, &cz1);

	cell_get_cartesian(reference, &ax2, &ay2, &az2,
	                              &bx2, &by2, &bz2,
	                              &cx2, &cy2, &cz2);

	if ( angle_between(ax1, ay1, az1, ax2, ay2, az2) > tols[3] ) return 0;
	if ( angle_between(bx1, by1, bz1, bx2, by2, bz2) > tols[4] ) return 0;
	if ( angle_between(cx1, cy1, cz1, cx2, cy2, cz2) > tols[5] ) return 0;

	if ( moduli_check(ax1, ay1, az1, ax2, ay2, az2) > tols[0] ) return 0;
	if ( moduli_check(bx1, by1, bz1, bx2, by2, bz2) > tols[1] ) return 0;
	if ( moduli_check(cx1, cy1, cz1, cx2, cy2, cz2) > tols[2] ) return 0;

	return 1;
}


/**
 * \param cell: A UnitCell
 * \param reference: Another UnitCell
 * \param tols: Pointer to six tolerance values (see below)
 * \param pmb: Place to store pointer to matrix
 *
 * Compare the two unit cells.  If, using any permutation of the axes, the
 * axes can be made to match in length and the axes aligned in space, this
 * function returns non-zero and stores the transformation which must be applied
 * to \p cell at \p pmb.
 *
 * A "permutation" means a transformation represented by a matrix with all
 * elements equal to +1, 0 or -1, having determinant +1 or -1.  That means that
 * this function will find the relationship between a left-handed and a right-
 * handed basis.
 *
 * Note that the orientations of the cells must match, not just the parameters.
 * The comparison is done after reindexing using
 * compare_cell_parameters_and_orientation().  See that function for more
 * details.
 *
 * \p cell and \p reference must have the same centering.  Otherwise, this
 * function always returns zero.
 *
 * \returns non-zero if the cells match.
 *
 */
/* 'tols' is in frac (not %) and radians */
int compare_permuted_cell_parameters_and_orientation(UnitCell *cell,
                                                     UnitCell *reference,
                                                     const double *tols,
                                                     IntegerMatrix **pmb)
{
	IntegerMatrix *m;
	int i[9];

	if ( cell_get_centering(cell) != cell_get_centering(reference) ) return 0;

	m = intmat_new(3, 3);

	for ( i[0]=-1; i[0]<=+1; i[0]++ ) {
	for ( i[1]=-1; i[1]<=+1; i[1]++ ) {
	for ( i[2]=-1; i[2]<=+1; i[2]++ ) {
	for ( i[3]=-1; i[3]<=+1; i[3]++ ) {
	for ( i[4]=-1; i[4]<=+1; i[4]++ ) {
	for ( i[5]=-1; i[5]<=+1; i[5]++ ) {
	for ( i[6]=-1; i[6]<=+1; i[6]++ ) {
	for ( i[7]=-1; i[7]<=+1; i[7]++ ) {
	for ( i[8]=-1; i[8]<=+1; i[8]++ ) {

		UnitCell *nc;
		int j, k;
		int l = 0;
		signed int det;

		for ( j=0; j<3; j++ )
			for ( k=0; k<3; k++ )
				intmat_set(m, j, k, i[l++]);

		det = intmat_det(m);
		if ( (det != +1) && (det != -1) ) continue;

		nc = cell_transform_intmat(cell, m);

		if ( compare_cell_parameters_and_orientation(nc, reference,
		                                             tols) )
		{
			if ( pmb != NULL ) *pmb = m;
			cell_free(nc);
			return 1;
		}

		cell_free(nc);

	}
	}
	}
	}
	}
	}
	}
	}
	}

	intmat_free(m);
	return 0;
}


struct cand
{
	Rational abc[3];
	double fom;
};


static int cmpcand(const void *av, const void *bv)
{
	const struct cand *a = av;
	const struct cand *b = bv;
	return a->fom > b->fom;
}


static Rational *find_candidates(double len, double *a, double *b, double *c,
                                 double ltl, int csl, int *pncand)
{
	Rational *r;
	struct cand *cands;
	const int max_cand = 1024;
	int ncand = 0;
	Rational *rat;
	int nrat;
	int nrej = 0;
	int ia, ib, ic;
	int i;

	cands = malloc(max_cand * sizeof(struct cand));
	if ( cands == NULL ) return NULL;

	rat = rtnl_list(-5, 5, 1, csl ? 4 : 1, &nrat);
	if ( rat == NULL ) return NULL;

	for ( ia=0; ia<nrat; ia++ ) {
	for ( ib=0; ib<nrat; ib++ ) {
	for ( ic=0; ic<nrat; ic++ ) {
		double vec[3];
		double abc[3];
		double veclen;
		abc[0] = rtnl_as_double(rat[ia]);
		abc[1] = rtnl_as_double(rat[ib]);
		abc[2] = rtnl_as_double(rat[ic]);
		vec[0] = a[0]*abc[0] + b[0]*abc[1] + c[0]*abc[2];
		vec[1] = a[1]*abc[0] + b[1]*abc[1] + c[1]*abc[2];
		vec[2] = a[2]*abc[0] + b[2]*abc[1] + c[2]*abc[2];
		veclen = modulus(vec[0], vec[1], vec[2]);
		if ( within_tolerance(len, veclen, ltl*100.0) ) {
			if ( ncand == max_cand ) {
				nrej++;
			} else {
				cands[ncand].abc[0] = rat[ia];
				cands[ncand].abc[1] = rat[ib];
				cands[ncand].abc[2] = rat[ic];
				cands[ncand].fom = fabs(veclen - len);
				ncand++;
			}
		}
	}
	}
	}

	if ( nrej ) {
		ERROR("WARNING: Too many vector candidates (%i rejected)\n", nrej);
	}

	/* Sort by difference from reference vector length */
	qsort(cands, ncand, sizeof(struct cand), cmpcand);

	r = malloc(ncand * 3 * sizeof(Rational));
	if ( r == 0 ) return NULL;

	for ( i=0; i<ncand; i++ ) {
		r[3*i+0] = cands[i].abc[0];
		r[3*i+1] = cands[i].abc[1];
		r[3*i+2] = cands[i].abc[2];
	}
	free(cands);

	*pncand = ncand;
	return r;
}


static double g6_distance(UnitCell *cell1, UnitCell *cell2)
{
	struct g6 g1, g2;
	double total = 0.0;
	g1 = cell_get_G6(cell1);
	g2 = cell_get_G6(cell2);
	total += (g1.A-g2.A)*(g1.A-g2.A);
	total += (g1.B-g2.B)*(g1.B-g2.B);
	total += (g1.C-g2.C)*(g1.C-g2.C);
	total += (g1.D-g2.D)*(g1.D-g2.D);
	total += (g1.E-g2.E)*(g1.E-g2.E);
	total += (g1.F-g2.F)*(g1.F-g2.F);
	return sqrt(total);
}


/**
 * \param cell_in: A UnitCell
 * \param reference_in: Another UnitCell
 * \param tols: Pointer to tolerances for a,b,c (fractional), al,be,ga (radians)
 * \param csl: Non-zero to look for coincidence site lattice relationships
 * \param pmb: Place to store pointer to matrix
 *
 * Compare the \p cell_in with \p reference_in.  If \p cell is a derivative lattice
 * of \p reference, within fractional axis length differences \p tols[0..2]
 * and absolute angle difference \p tols[3..5] (in radians), this function returns
 * non-zero and stores the transformation which needs to be applied to
 * \p cell_in at \p pmb.
 *
 * Note that the tolerances will be applied to the primitive unit cell.  If
 * the reference cell is centered, a primitive unit cell will first be calculated.
 *
 * Subject to the tolerances, this function will find the transformation which
 * gives the best match to the reference cell, using the Euclidian norm in
 * G6 [see e.g. Andrews and Bernstein, Acta Cryst. A44 (1988) p1009].
 *
 * Only the cell parameters will be compared.  The relative orientations are
 * irrelevant.
 *
 * If \p csl is zero, the lattices must be derivatives of one another.  If
 * non-zero, a coincidence site lattice relationship will be searched for,
 * meaning that the lattice points of the transformed version of \p cell_in
 * might not coincide with lattice points of \p reference_in.
 *
 * This function is used by CrystFEL's cell_tool program to find non-obvious
 * relationships between crystal lattices.  For most routine comparisons, this
 * function is probably not the one you need!
 *
 * \returns non-zero if the cells match, zero for no match or error.
 *
 */
int compare_derivative_cell_parameters(UnitCell *cell_in, UnitCell *reference_in,
                                       const double *tols, int csl,
                                       RationalMatrix **pmb)
{
	UnitCell *cell;
	UnitCell *reference;
	IntegerMatrix *CBint;
	RationalMatrix *CiA;
	RationalMatrix *CB;
	RationalMatrix *M;
	double a, b, c, al, be, ga;
	double av[3], bv[3], cv[3];
	Rational *cand_a;
	Rational *cand_b;
	Rational *cand_c;
	int ncand_a, ncand_b, ncand_c;
	int ia, ib;
	RationalMatrix *CiAMCB = NULL;
	double min_dist = +INFINITY;

	/* Actually compare against primitive version of reference */
	reference = uncenter_cell(reference_in, &CBint, NULL);
	if ( reference == NULL ) return 0;
	CB = rtnl_mtx_from_intmat(CBint);
	intmat_free(CBint);

	/* Actually compare primitive version of cell */
	cell = uncenter_cell(cell_in, NULL, &CiA);
	if ( cell == NULL ) return 0;

	/* Get target parameters */
	cell_get_parameters(reference, &a, &b, &c, &al, &be, &ga);
	cell_get_cartesian(cell, &av[0], &av[1], &av[2],
	                         &bv[0], &bv[1], &bv[2],
	                         &cv[0], &cv[1], &cv[2]);

	/* Find vectors in 'cell' with lengths close to a, b and c */
	cand_a = find_candidates(a, av, bv, cv, tols[0], csl, &ncand_a);
	cand_b = find_candidates(b, av, bv, cv, tols[1], csl, &ncand_b);
	cand_c = find_candidates(c, av, bv, cv, tols[2], csl, &ncand_c);

	if ( (ncand_a==0) || (ncand_b==0) || (ncand_c==0) ) {
		*pmb = NULL;
		cell_free(cell);
		cell_free(reference);
		rtnl_mtx_free(CB);
		rtnl_mtx_free(CiA);
		return 0;
	}

	M = rtnl_mtx_new(3, 3);
	for ( ia=0; ia<ncand_a; ia++ ) {
		for ( ib=0; ib<ncand_b; ib++ ) {

			UnitCell *test;
			double at, bt, ct, alt, bet, gat;
			double dist;
			int ic = 0;
			RationalMatrix *MCB;

			/* Form the matrix using the first candidate for c */
			rtnl_mtx_set(M, 0, 0, cand_a[3*ia+0]);
			rtnl_mtx_set(M, 1, 0, cand_a[3*ia+1]);
			rtnl_mtx_set(M, 2, 0, cand_a[3*ia+2]);
			rtnl_mtx_set(M, 0, 1, cand_b[3*ib+0]);
			rtnl_mtx_set(M, 1, 1, cand_b[3*ib+1]);
			rtnl_mtx_set(M, 2, 1, cand_b[3*ib+2]);
			rtnl_mtx_set(M, 0, 2, cand_c[3*ic+0]);
			rtnl_mtx_set(M, 1, 2, cand_c[3*ic+1]);
			rtnl_mtx_set(M, 2, 2, cand_c[3*ic+2]);

			/* Check angle between a and b */
			test = cell_transform_rational(cell, M);
			cell_get_parameters(test, &at, &bt, &ct, &alt, &bet, &gat);
			cell_free(test);
			if ( fabs(gat - ga) > tols[5] ) continue;

			/* Gamma OK, now look for place for c axis */
			for ( ic=0; ic<ncand_c; ic++ ) {

				rtnl_mtx_set(M, 0, 0, cand_a[3*ia+0]);
				rtnl_mtx_set(M, 1, 0, cand_a[3*ia+1]);
				rtnl_mtx_set(M, 2, 0, cand_a[3*ia+2]);
				rtnl_mtx_set(M, 0, 1, cand_b[3*ib+0]);
				rtnl_mtx_set(M, 1, 1, cand_b[3*ib+1]);
				rtnl_mtx_set(M, 2, 1, cand_b[3*ib+2]);
				rtnl_mtx_set(M, 0, 2, cand_c[3*ic+0]);
				rtnl_mtx_set(M, 1, 2, cand_c[3*ic+1]);
				rtnl_mtx_set(M, 2, 2, cand_c[3*ic+2]);

				if ( rtnl_cmp(rtnl_mtx_det(M),rtnl_zero()) == 0 ) continue;

				test = cell_transform_rational(cell, M);

				if ( !csl && (cell_get_centering(test) != 'P') ) continue;

				cell_get_parameters(test, &at, &bt, &ct, &alt, &bet, &gat);
				if ( !right_handed(test) ) {
					cell_free(test);
					continue;
				}
				if ( fabs(alt - al) > tols[3] ) {
					cell_free(test);
					continue;
				}
				if ( fabs(bet - be) > tols[4] ) {
					cell_free(test);
					continue;
				}

				dist = g6_distance(test, reference);
				if ( dist < min_dist ) {
					min_dist = dist;
					MCB = rtnlmtx_times_rtnlmtx(M, CB);
					CiAMCB = rtnlmtx_times_rtnlmtx(CiA, MCB);
					rtnl_mtx_free(MCB);
				}

				cell_free(test);

			}
		}
	}

	rtnl_mtx_free(M);
	free(cand_a);
	free(cand_b);
	free(cand_c);

	if ( CiAMCB == NULL ) {
		*pmb = NULL;
		return 0;
	}

	/* Solution found */
	*pmb = CiAMCB;
	return 1;
}


/* Criteria from Grosse-Kunstleve et al., Acta Cryst A60 (2004) p1-6 */
#define GT(x,y) (y < x - eps)
#define LT(x,y) GT(y,x)
#define EQ(x,y) !(GT(x,y) || LT(x,y))
#define LTE(x,y) !(GT(x,y))


static int in_standard_presentation(struct g6 g, double eps)
{
	if ( GT(g.A, g.B) ) return 0;
	if ( GT(g.B, g.C) ) return 0;

	if ( EQ(g.A, g.B) && GT(fabs(g.D), fabs(g.E)) ) return 0;
	if ( EQ(g.B, g.C) && GT(fabs(g.E), fabs(g.F)) ) return 0;

	if ( ( GT(g.D, 0.0) && GT(g.E, 0.0) && GT(g.F, 0.0) )
	  || ( !GT(g.D, 0.0) && !GT(g.E, 0.0) && !GT(g.F, 0.0) ) ) return 1;

	return 0;
}


static int is_burger(struct g6 g, double eps)
{
	if ( !in_standard_presentation(g, eps) ) return 0;
	if ( GT(fabs(g.D), g.B) ) return 0;
	if ( GT(fabs(g.E), g.A) ) return 0;
	if ( GT(fabs(g.F), g.A) ) return 0;
	if ( LT(g.D + g.E + g.F + g.A + g.B, 0.0) ) return 0;
	return 1;
}


static int UNUSED is_niggli(struct g6 g, double eps)
{
	if ( !is_burger(g, eps) ) return 0;

	if ( EQ(g.D, g.B) && GT(g.F, 2.0*g.E) ) return 0;
	if ( EQ(g.E, g.A) && GT(g.F, 2.0*g.D) ) return 0;
	if ( EQ(g.F, g.A) && GT(g.E, 2.0*g.D) ) return 0;

	if ( EQ(g.D, -g.B) && !EQ(g.F, 0.0) ) return 0;
	if ( EQ(g.E, -g.A) && !EQ(g.F, 0.0) ) return 0;
	if ( EQ(g.F, -g.A) && !EQ(g.E, 0.0) ) return 0;

	if ( EQ(g.D + g.E + g.F + g.A + g.B, 0.0)
	  && GT(2.0*(g.A + g.E)+g.F, 0.0) ) return 0;

	return 1;
}


static int DEF_positive(struct g6 g, double eps)
{
	int n_zero = 0;
	int n_positive = 0;

	if ( LT(0.0, g.D) ) {
		n_positive++;
	} else {
		if ( !LT(g.D, 0.0) ) {
			n_zero++;
		}
	}

	if ( LT(0.0, g.E) ) {
		n_positive++;
	} else {
		if ( !LT(g.E, 0.0) ) {
			n_zero++;
		}
	}

	if ( LT(0.0, g.F) ) {
		n_positive++;
	} else {
		if ( !LT(g.F, 0.0) ) {
			n_zero++;
		}
	}

	return (n_positive==3) || ((n_zero==0) && (n_positive==1));
}


static void debug_lattice(struct g6 g, double eps, int step)
{
#if 0
	STATUS("After step %i:  %e %e %e %e %e %e --", step,
	       g.A/1e-0, g.B/1e-0, g.C/1e-0,
	       g.D/1e-0, g.E/1e-0, g.F/1e-0);
	if ( is_burger(g, eps) ) STATUS(" B");
	if ( is_niggli(g, eps) ) STATUS(" N");
	if ( eps_sign(g.D, eps)*eps_sign(g.E, eps)*eps_sign(g.F, eps) > 0 ) {
		STATUS(" I");
	} else {
		STATUS(" II");
	}
	STATUS("\n");
#endif
}


static void mult_in_place(IntegerMatrix *T, IntegerMatrix *M)
{
	int i, j;
	IntegerMatrix *tmp = intmat_times_intmat(T, M);
	assert(intmat_det(M) == 1);
	assert(intmat_det(T) == 1);
	assert(intmat_det(tmp) == 1);
	for ( i=0; i<3; i++ ) {
		for ( j=0; j<3; j++ ) {
			intmat_set(T, i, j, intmat_get(tmp, i, j));
		}
	}
	intmat_free(tmp);
}


/* Cell volume from G6 components */
static double g6_volume(struct g6 g)
{
	return sqrt(g.A*g.B*g.C
	       - 0.25*(g.A*g.D*g.D + g.B*g.E*g.E + g.C*g.F*g.F - g.D*g.E*g.F));
}


/* NB The G6 components are passed by value, not reference.
 * It's the caller's reponsibility to apply the matrix to the cell and
 * re-calculate the G6 vector, if required. */
IntegerMatrix *reduce_g6(struct g6 g, double epsrel)
{
	IntegerMatrix *T;
	IntegerMatrix *M;
	int finished;
	double eps;

	eps = pow(g6_volume(g), 1.0/3.0) * epsrel;
	eps = eps*eps;

	T = intmat_identity(3);
	M = intmat_new(3, 3);

	debug_lattice(g, eps, 0);

	do {

		int done_s1s2;

		do {
			done_s1s2 = 1;

			if ( GT(g.A, g.B)
			  || (EQ(g.A, g.B) && GT(fabs(g.D), fabs(g.E))) )
			{
				/* Swap a and b */
				double temp;
				temp = g.A; g.A = g.B;  g.B = temp;
				temp = g.D; g.D = g.E;  g.E = temp;

				intmat_zero(M);
				intmat_set(M, 1, 0, -1);
				intmat_set(M, 0, 1, -1);
				intmat_set(M, 2, 2, -1);
				mult_in_place(T, M);

				debug_lattice(g, eps, 1);

			}

			if ( GT(g.B, g.C)
			  || (EQ(g.B, g.C) && GT(fabs(g.E), fabs(g.F))) )
			{
				/* Swap b and c */
				double temp;
				temp = g.B; g.B = g.C;  g.C = temp;
				temp = g.E; g.E = g.F;  g.F = temp;

				intmat_zero(M);
				intmat_set(M, 0, 0, -1);
				intmat_set(M, 1, 2, -1);
				intmat_set(M, 2, 1, -1);
				mult_in_place(T, M);

				debug_lattice(g, eps, 2);

				/* ..."and go to 1." */
				done_s1s2 = 0;
			}

		} while ( !done_s1s2 );

		finished = 0;

		if ( DEF_positive(g, eps) ) {

			intmat_zero(M);
			intmat_set(M, 0, 0, LT(g.D, 0.0) ? -1 : 1);
			intmat_set(M, 1, 1, LT(g.E, 0.0) ? -1 : 1);
			intmat_set(M, 2, 2, LT(g.F, 0.0) ? -1 : 1);
			mult_in_place(T, M);

			assert(intmat_det(M) == 1);
			g.D = fabs(g.D);
			g.E = fabs(g.E);
			g.F = fabs(g.F);

			debug_lattice(g, eps, 3);

		} else  {

			int z = 4;
			signed int ijk[3] = { 1, 1, 1 };

			if ( GT(g.D, 0.0) ) {
				ijk[0] = -1;
			} else if ( !LT(g.D, 0.0) ) {
				z = 0;
			}
			if ( GT(g.E, 0.0) ) {
				ijk[1] = -1;
			} else if ( !LT(g.E, 0.0) ) {
				z = 1;
			}
			if ( GT(g.F, 0.0) ) {
				ijk[2] = -1;
			} else if ( !LT(g.F, 0.0) ) {
				z = 2;
			}
			if ( ijk[0]*ijk[1]*ijk[2] < 0 ) {
				if ( z == 4 ) {
					ERROR("No element in reduction step 4\n");
					abort();
				}
				ijk[z] = -1;
			}
			intmat_zero(M);
			intmat_set(M, 0, 0, ijk[0]);
			intmat_set(M, 1, 1, ijk[1]);
			intmat_set(M, 2, 2, ijk[2]);
			mult_in_place(T, M);

			g.D = -fabs(g.D);
			g.E = -fabs(g.E);
			g.F = -fabs(g.F);

			debug_lattice(g, eps, 4);

		}

		if ( GT(fabs(g.D), g.B)
		  || (EQ(g.B, g.D) && LT(2.0*g.E, g.F))
		  || (EQ(g.B, -g.D) && LT(g.F, 0.0)) )
		{
			signed int s = g.D > 0.0  ? 1 : -1;

			intmat_zero(M);
			intmat_set(M, 0, 0, 1);
			intmat_set(M, 1, 1, 1);
			intmat_set(M, 2, 2, 1);
			intmat_set(M, 1, 2, -s);
			mult_in_place(T, M);

			g.C = g.B + g.C -s*g.D;
			g.D = -2*s*g.B + g.D;
			g.E = g.E - s*g.F;

			debug_lattice(g, eps, 5);

		} else if ( GT(fabs(g.E), g.A)
		         || (EQ(g.A, g.E) && LT(2.0*g.D, g.F))
		         || (EQ(g.A, -g.E) && LT(g.F, 0.0)) )
		{
			signed int s = g.E > 0.0  ? 1 : -1;

			intmat_zero(M);
			intmat_set(M, 0, 0, 1);
			intmat_set(M, 1, 1, 1);
			intmat_set(M, 2, 2, 1);
			intmat_set(M, 0, 2, -s);
			mult_in_place(T, M);

			g.C = g.A + g.C -s*g.E;
			g.D = g.D - s*g.F;
			g.E = -2*s*g.A + g.E;

			debug_lattice(g, eps, 6);

		} else if ( GT(fabs(g.F), g.A)
		         || (EQ(g.A, g.F) && LT(2.0*g.D, g.E))
		         || (EQ(g.A, -g.F) && LT(g.E, 0.0)) )
		{
			signed int s = g.F > 0.0  ? 1 : -1;

			intmat_zero(M);
			intmat_set(M, 0, 0, 1);
			intmat_set(M, 1, 1, 1);
			intmat_set(M, 2, 2, 1);
			intmat_set(M, 0, 1, -s);
			mult_in_place(T, M);

			g.B = g.A + g.B -s*g.F;
			g.D = g.D - s*g.E;
			g.F = -2*s*g.A + g.F;

			debug_lattice(g, eps, 7);

		} else if ( LT(g.A+g.B+g.D+g.E+g.F, 0.0)  /* not g.C */
		         || ( (EQ(g.A+g.B+g.D+g.E+g.F, 0.0)
		               && GT(2.0*g.A + 2.0*g.E + g.F, 0.0)) ) )
		{
			intmat_zero(M);
			intmat_set(M, 0, 0, 1);
			intmat_set(M, 1, 1, 1);
			intmat_set(M, 2, 2, 1);
			intmat_set(M, 1, 2, 1);
			mult_in_place(T, M);

			g.C = g.A+g.B+g.C+g.D+g.E+g.F;
			g.D = 2.0*g.B + g.D + g.F;
			g.E = 2.0*g.A + g.E + g.F;

			debug_lattice(g, eps, 8);

		} else {
			finished = 1;
		}

	} while ( !finished );

	debug_lattice(g, eps, 99);

	assert(is_burger(g, eps));
	assert(is_niggli(g, eps));

	intmat_free(M);
	return T;
}


static double cell_diff(UnitCell *cell, double a, double b, double c,
                        double al, double be, double ga)
{
	double ta, tb, tc, tal, tbe, tga;
	double diff = 0.0;
	cell_get_parameters(cell, &ta, &tb, &tc, &tal, &tbe, &tga);
	diff += fabs(a - ta);
	diff += fabs(b - tb);
	diff += fabs(c - tc);
	return diff;
}


static int random_int(int max)
{
	int r;
	int limit = RAND_MAX;
	while ( limit % max ) limit--;
	do {
		r = rand();
	} while ( r > limit );
	return rand() % max;
}


static IntegerMatrix *check_permutations(UnitCell *cell_reduced, UnitCell *reference,
                                         RationalMatrix *CiARA, IntegerMatrix *RiBCB,
                                         const double *tols)
{
	IntegerMatrix *m;
	int i[9];
	double a, b, c, al, be, ga;
	double min_dist = +INFINITY;
	int s, sel;
	IntegerMatrix *best_m[24];
	int n_best = 0;

	m = intmat_new(3, 3);
	cell_get_parameters(reference, &a, &b, &c, &al, &be, &ga);

	for ( i[0]=-1; i[0]<=+1; i[0]++ ) {
	for ( i[1]=-1; i[1]<=+1; i[1]++ ) {
	for ( i[2]=-1; i[2]<=+1; i[2]++ ) {
	for ( i[3]=-1; i[3]<=+1; i[3]++ ) {
	for ( i[4]=-1; i[4]<=+1; i[4]++ ) {
	for ( i[5]=-1; i[5]<=+1; i[5]++ ) {
	for ( i[6]=-1; i[6]<=+1; i[6]++ ) {
	for ( i[7]=-1; i[7]<=+1; i[7]++ ) {
	for ( i[8]=-1; i[8]<=+1; i[8]++ ) {

		UnitCell *nc;
		int j, k;
		int l = 0;
		signed int det;
		UnitCell *tmp;

		for ( j=0; j<3; j++ )
			for ( k=0; k<3; k++ )
				intmat_set(m, j, k, i[l++]);

		det = intmat_det(m);
		if ( det != +1 ) continue;

		tmp = cell_transform_intmat(cell_reduced, m);
		nc = cell_transform_intmat(tmp, RiBCB);
		cell_free(tmp);

		if ( compare_cell_parameters(nc, reference, tols) ) {
			double dist = cell_diff(nc, a, b, c, al, be, ga);
			if ( dist < min_dist ) {

				/* If the new solution is significantly better,
				 * dump all the previous ones */
				for ( s=0; s<n_best; s++ ) {
					intmat_free(best_m[s]);
				}
				min_dist = dist;
				best_m[0] = intmat_copy(m);
				n_best = 1;

			} else if ( dist == min_dist ) {

				if ( n_best == 24 ) {
					ERROR("WARNING: Too many equivalent "
					      "reindexed lattices\n");
				} else {
					/* If the new solution is the same as the
					* previous one, add it to the list */
					best_m[n_best++] = intmat_copy(m);
				}

			}
		}

		cell_free(nc);

	}
	}
	}
	}
	}
	}
	}
	}
	}

	intmat_free(m);

	if ( n_best == 0 ) return NULL;

	sel = n_best;
	if ( n_best == 1 ) {

		/* If there's one solution, choose that one, of course */
		sel = 0;

	} else {

		/* If one of the solutions results in an identity applied to the
		 * original cell, choose that one */

		for ( s=0; s<n_best; s++ ) {
			RationalMatrix *tmp;
			RationalMatrix *comb;
			tmp = rtnlmtx_times_intmat(CiARA, best_m[s]);
			comb = rtnlmtx_times_intmat(tmp, RiBCB);
			if ( rtnl_mtx_is_identity(comb) ) {
				sel = s;
			}
			rtnl_mtx_free(tmp);
			rtnl_mtx_free(comb);
		}

	}

	/* Still undecided?  Choose randomly, to avoid weird distributions
	 * in the cell parameters */
	if ( sel == n_best ) {
		sel = random_int(n_best);
	}

	/* Free all the others */
	for ( s=0; s<n_best; s++ ) {
		if ( s != sel ) intmat_free(best_m[s]);
	}

	return best_m[sel];
}


/**
 * \param cell_in: A UnitCell
 * \param reference_in: Another UnitCell
 * \param tols: Pointer to tolerances for a,b,c (fractional), al,be,ga (radians)
 * \param pmb: Place to store pointer to matrix, or NULL if not needed
 *
 * Compare the \p cell_in with \p reference_in.  If they represent the same
 * lattice, this function returns a copy of \p cell_in transformed to look
 * similar to \p reference_in.  Otherwise, it returns NULL.
 *
 * If \pmb is non-NULL, the transformation which needs to be applied to
 * \p cell_in will be stored there.
 *
 * Only the cell parameters will be compared.  The relative orientations are
 * irrelevant.  The tolerances will be applied to the transformed copy of
 * \p cell_in, i.e. the version of the input cell which looks similar to
 * \p reference_in.  Subject to the tolerances, the cell will be chosen which
 * has the lowest total absolute error in unit cell axis lengths.
 *
 * There will usually be several transformation matrices which produce exactly
 * the same total absolute error.  If one of the matrices is an identity, that
 * one will be used.  Otherwise, the matrix will be selected at random from the
 * possibilities.  This avoids skewed distributions of unit cell parameters,
 * e.g. the angles always being greater than 90 degrees.
 *
 * This is the right function to use for deciding if an indexing solution
 * matches a reference cell or not.
 *
 * \returns A newly allocated UnitCell, or NULL.
 *
 */
UnitCell *compare_reindexed_cell_parameters(UnitCell *cell_in,
                                            UnitCell *reference_in,
                                            const double *tols,
                                            RationalMatrix **pmb)
{
	UnitCell *cell;
	UnitCell *reference;
	IntegerMatrix *CB;
	RationalMatrix *CiA;
	IntegerMatrix *RA;
	IntegerMatrix *RB;
	IntegerMatrix *RiB;
	IntegerMatrix *P;
	RationalMatrix *CiARA;
	IntegerMatrix *RiBCB;
	UnitCell *cell_reduced;
	UnitCell *reference_reduced;
	UnitCell *match;
	struct g6 g6cell;
	struct g6 g6ref;

	//STATUS("The input cell:\n");
	//cell_print(cell_in);

	//STATUS("The reference cell:\n");
	//cell_print(reference_in);

	/* Un-center both cells */
	reference = uncenter_cell(reference_in, &CB, NULL);
	if ( reference == NULL ) return NULL;

	cell = uncenter_cell(cell_in, NULL, &CiA);
	if ( cell == NULL ) return NULL;

	/* Calculate G6 components for both */
	g6cell = cell_get_G6(cell);
	g6ref = cell_get_G6(reference);

	/* Convert both to reduced basis (stably) */
	RA = reduce_g6(g6cell, 1e-5);
	//STATUS("------------------------------------------\n");
	RB = reduce_g6(g6ref, 1e-5);

	cell_reduced = cell_transform_intmat(cell, RA);

	//STATUS("Reduced cell:\n");
	//cell_print(cell_reduced);
	//STATUS("Reduced reference:\n");
	reference_reduced = cell_transform_intmat(reference, RB);
	//cell_print(reference_reduced);

	/* The primitive (non-reduced) cells are no longer needed */
	cell_free(reference);
	cell_free(cell);

	/* Within tolerance? */
	RiB = intmat_inverse(RB);
	RiBCB = intmat_times_intmat(RiB, CB);
	intmat_free(RiB);
	CiARA = rtnlmtx_times_intmat(CiA, RA);
	P = check_permutations(cell_reduced, reference_in, CiARA, RiBCB, tols);
	if ( P != NULL ) {

		RationalMatrix *tmp;
		RationalMatrix *comb;

		//STATUS("The best permutation matrix:\n");
		//intmat_print(P);

		/* Calculate combined matrix: CiA.RA.P.RiB.CB */
		tmp = rtnlmtx_times_intmat(CiARA, P);
		comb = rtnlmtx_times_intmat(tmp, RiBCB);
		rtnl_mtx_free(tmp);

		match = cell_transform_rational(cell_in, comb);
		//STATUS("Original cell transformed to look like reference:\n");
		//cell_print(match);

		if ( pmb != NULL ) *pmb = comb;

	} else {
		match = NULL;
	}

	rtnl_mtx_free(CiARA);
	rtnl_mtx_free(CiA);
	intmat_free(RiBCB);
	intmat_free(RA);
	intmat_free(RB);
	intmat_free(CB);
	intmat_free(P);
	cell_free(reference_reduced);
	cell_free(cell_reduced);

	return match;
}
