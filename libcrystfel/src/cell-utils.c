/*
 * cell-utils.c
 *
 * Unit Cell utility functions
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2012,2014 Thomas White <taw@physics.org>
 *   2012           Lorenzo Galli
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>

#include "cell.h"
#include "cell-utils.h"
#include "utils.h"
#include "image.h"


/**
 * SECTION:cell-utils
 * @short_description: Unit cell utilities
 * @title: Unit cell utilities
 * @section_id:
 * @see_also:
 * @include: "cell-utils.h"
 * @Image:
 *
 * There are some utility functions associated with the core %UnitCell.
 **/


/* Weighting factor of lengths relative to angles */
#define LWEIGHT (10.0e-9)


/**
 * cell_rotate:
 * @in: A %UnitCell to rotate
 * @quat: A %quaternion
 *
 * Rotate a %UnitCell using a %quaternion.
 *
 * Returns: a newly allocated rotated copy of @in.
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

	assert(rh_reciprocal == rh_direct);

	return rh_direct;
}


void cell_print(UnitCell *cell)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, alpha, beta, gamma;
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	LatticeType lt;
	char cen;

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
		cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

		STATUS("a      b      c            alpha   beta  gamma\n");
		STATUS("%6.2f %6.2f %6.2f A    %6.2f %6.2f %6.2f deg\n",
		       a*1e10, b*1e10, c*1e10,
		       rad2deg(alpha), rad2deg(beta), rad2deg(gamma));

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

		STATUS("Cell representation is %s.\n", cell_rep(cell));

	} else {
		STATUS("Unit cell parameters are not specified.\n");
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


static UnitCellTransformation *uncentering_transformation(UnitCell *in,
                                                          char *new_centering,
                                                          LatticeType *new_latt)
{
	UnitCellTransformation *t;
	const double OT = 1.0/3.0;
	const double TT = 2.0/3.0;
	const double H = 0.5;
	LatticeType lt;
	char ua, cen;

	lt = cell_get_lattice_type(in);
	ua = cell_get_unique_axis(in);
	cen = cell_get_centering(in);

	t = tfn_identity();
	if ( t == NULL ) return NULL;

	if ( ua == 'a' ) {
		tfn_combine(t, tfn_vector(0,1,0),
		               tfn_vector(0,0,1),
		               tfn_vector(1,0,0));
		if ( lt == L_MONOCLINIC ) {
			assert(cen != 'A');
			switch ( cen ) {
				case 'B' : cen = 'A'; break;
				case 'C' : cen = 'B'; break;
				case 'I' : cen = 'I'; break;
			}
		}
	}

	if ( ua == 'b' ) {
		tfn_combine(t, tfn_vector(0,0,1),
		               tfn_vector(1,0,0),
		               tfn_vector(0,1,0));
		if ( lt == L_MONOCLINIC ) {
			assert(cen != 'B');
			switch ( cen ) {
				case 'C' : cen = 'A'; break;
				case 'A' : cen = 'B'; break;
				case 'I' : cen = 'I'; break;
			}
		}
	}

	switch ( cen ) {

		case 'P' :
		*new_latt = lt;
		*new_centering = 'P';
		break;

		case 'R' :
		*new_latt = L_RHOMBOHEDRAL;
		*new_centering = 'R';
		break;

		case 'I' :
		tfn_combine(t, tfn_vector(-H,H,H),
		               tfn_vector(H,-H,H),
		               tfn_vector(H,H,-H));
		if ( lt == L_CUBIC ) {
			*new_latt = L_RHOMBOHEDRAL;
			*new_centering = 'R';
		} else {
			/* Tetragonal or orthorhombic */
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
		}
		break;

		case 'F' :
		tfn_combine(t, tfn_vector(0,H,H),
		               tfn_vector(H,0,H),
		               tfn_vector(H,H,0));
		if ( lt == L_CUBIC ) {
			*new_latt = L_RHOMBOHEDRAL;
			*new_centering = 'R';
		} else {
			assert(lt == L_ORTHORHOMBIC);
			*new_latt = L_TRICLINIC;
			*new_centering = 'P';
		}
		break;

		case 'A' :
		tfn_combine(t, tfn_vector( 1, 0, 0),
		               tfn_vector( 0, H, H),
		               tfn_vector( 0,-H, H));
		if ( lt == L_ORTHORHOMBIC ) {
			*new_latt = L_MONOCLINIC;
		} else {
			*new_latt = L_TRICLINIC;
		}
		*new_centering = 'P';
		break;

		case 'B' :
		tfn_combine(t, tfn_vector( H, 0, H),
		               tfn_vector( 0, 1, 0),
		               tfn_vector(-H, 0, H));
		if ( lt == L_ORTHORHOMBIC ) {
			*new_latt = L_MONOCLINIC;
		} else {
			*new_latt = L_TRICLINIC;
		}
		*new_centering = 'P';
		break;

		case 'C' :
		tfn_combine(t, tfn_vector( H, H, 0),
		               tfn_vector(-H, H, 0),
		               tfn_vector( 0, 0, 1));
		if ( lt == L_ORTHORHOMBIC ) {
			*new_latt = L_MONOCLINIC;
		} else {
			*new_latt = L_TRICLINIC;
		}
		*new_centering = 'P';
		break;

		case 'H' :
		/* Obverse setting */
		tfn_combine(t, tfn_vector(TT,OT,OT),
		               tfn_vector(-OT,OT,OT),
		               tfn_vector(-OT,-TT,OT));
		assert(lt == L_HEXAGONAL);
		*new_latt = L_RHOMBOHEDRAL;
		*new_centering = 'R';
		break;

		default :
		ERROR("Invalid centering '%c'\n", cell_get_centering(in));
		return NULL;

	}

	/* Reverse the axis permutation, but only if this was not an H->R
	 * transformation */
	if ( !((cen=='H') && (*new_latt == L_RHOMBOHEDRAL)) ) {
		if ( ua == 'a' ) {
			tfn_combine(t, tfn_vector(0,0,1),
				       tfn_vector(1,0,0),
				       tfn_vector(0,1,0));
		}

		if ( ua == 'b' ) {
			tfn_combine(t, tfn_vector(0,1,0),
				       tfn_vector(0,0,1),
				       tfn_vector(1,0,0));
		}
	}

	return t;
}


/**
 * uncenter_cell:
 * @in: A %UnitCell
 * @t: Location at which to store the transformation which was used.
 *
 * Turns any cell into a primitive one, e.g. for comparison purposes.  The
 * transformation which was used is stored at @t, which can be NULL if the
 * transformation is not required.
 *
 * Returns: a primitive version of @in in a conventional (unique axis c)
 * setting.
 *
 */
UnitCell *uncenter_cell(UnitCell *in, UnitCellTransformation **t)
{
	UnitCellTransformation *tt;
	char new_centering;
	LatticeType new_latt;
	UnitCell *out;

	if ( !bravais_lattice(in) ) {
		ERROR("Cannot uncenter: not a Bravais lattice.\n");
		cell_print(in);
		return NULL;
	}

	tt = uncentering_transformation(in, &new_centering, &new_latt);
	if ( tt == NULL ) return NULL;

	out = cell_transform(in, tt);
	if ( out == NULL ) return NULL;

	cell_set_lattice_type(out, new_latt);
	cell_set_centering(out, new_centering);

	if ( t != NULL ) {
		*t = tt;
	} else {
		tfn_free(tt);
	}

	return out;
}


#define MAX_CAND (1024)

static int right_handed_vec(struct rvec a, struct rvec b, struct rvec c)
{
	struct rvec aCb;
	double aCb_dot_c;

	/* "a" cross "b" */
	aCb.u = a.v*b.w - a.w*b.v;
	aCb.v = - (a.u*b.w - a.w*b.u);
	aCb.w = a.u*b.v - a.v*b.u;

	/* "a cross b" dot "c" */
	aCb_dot_c = aCb.u*c.u + aCb.v*c.v + aCb.w*c.w;

	if ( aCb_dot_c > 0.0 ) return 1;
	return 0;
}


struct cvec {
	struct rvec vec;
	float na;
	float nb;
	float nc;
	float fom;
};


static int same_vector(struct cvec a, struct cvec b)
{
	if ( a.na != b.na ) return 0;
	if ( a.nb != b.nb ) return 0;
	if ( a.nc != b.nc ) return 0;
	return 1;
}


/* Attempt to make 'cell' fit into 'template' somehow */
UnitCell *match_cell(UnitCell *cell_in, UnitCell *template_in, int verbose,
                     const float *tols, int reduce)
{
	signed int n1l, n2l, n3l;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	int i, j;
	double lengths[3];
	double angles[3];
	struct cvec *cand[3];
	UnitCell *new_cell = NULL;
	float best_fom = +999999999.9; /* Large number.. */
	int ncand[3] = {0,0,0};
	signed int ilow, ihigh;
	float angtol = deg2rad(tols[3]);
	UnitCell *cell;
	UnitCell *template;
	UnitCellTransformation *uncentering;
	UnitCell *new_cell_trans;

	/* "Un-center" the template unit cell to make the comparison easier */
	template = uncenter_cell(template_in, &uncentering);

	/* The candidate cell is also uncentered, because it might be centered
	 * if it came from (e.g.) MOSFLM */
	cell = uncenter_cell(cell_in, NULL);

	if ( cell_get_reciprocal(template, &asx, &asy, &asz,
	                         &bsx, &bsy, &bsz,
	                         &csx, &csy, &csz) ) {
		ERROR("Couldn't get reciprocal cell for template.\n");
		cell_free(template);
		cell_free(cell);
		tfn_free(uncentering);
		return NULL;
	}

	lengths[0] = modulus(asx, asy, asz);
	lengths[1] = modulus(bsx, bsy, bsz);
	lengths[2] = modulus(csx, csy, csz);

	angles[0] = angle_between(bsx, bsy, bsz, csx, csy, csz);
	angles[1] = angle_between(asx, asy, asz, csx, csy, csz);
	angles[2] = angle_between(asx, asy, asz, bsx, bsy, bsz);

	cand[0] = malloc(MAX_CAND*sizeof(struct cvec));
	cand[1] = malloc(MAX_CAND*sizeof(struct cvec));
	cand[2] = malloc(MAX_CAND*sizeof(struct cvec));

	if ( cell_get_reciprocal(cell, &asx, &asy, &asz,
	                         &bsx, &bsy, &bsz,
	                         &csx, &csy, &csz) ) {
		ERROR("Couldn't get reciprocal cell.\n");
		cell_free(template);
		cell_free(cell);
		tfn_free(uncentering);
		return NULL;
	}

	if ( reduce ) {
		ilow = -2;  ihigh = 4;
	} else {
		ilow = 0;  ihigh = 1;
	}

	/* Negative values mean 1/n, positive means n, zero means zero */
	for ( n1l=ilow; n1l<=ihigh; n1l++ ) {
	for ( n2l=ilow; n2l<=ihigh; n2l++ ) {
	for ( n3l=ilow; n3l<=ihigh; n3l++ ) {

		float n1, n2, n3;
		signed int b1, b2, b3;

		n1 = (n1l>=0) ? (n1l) : (1.0/n1l);
		n2 = (n2l>=0) ? (n2l) : (1.0/n2l);
		n3 = (n3l>=0) ? (n3l) : (1.0/n3l);

		if ( !reduce ) {
			if ( n1l + n2l + n3l > 1 ) continue;
		}

		/* 'bit' values can be +1 or -1 */
		for ( b1=-1; b1<=1; b1+=2 ) {
		for ( b2=-1; b2<=1; b2+=2 ) {
		for ( b3=-1; b3<=1; b3+=2 ) {

			double tx, ty, tz;
			double tlen;
			int i;

			n1 *= b1;  n2 *= b2;  n3 *= b3;

			tx = n1*asx + n2*bsx + n3*csx;
			ty = n1*asy + n2*bsy + n3*csy;
			tz = n1*asz + n2*bsz + n3*csz;
			tlen = modulus(tx, ty, tz);

			/* Test modulus for agreement with moduli of template */
			for ( i=0; i<3; i++ ) {

				if ( !within_tolerance(lengths[i], tlen,
				                       tols[i]) )
				{
					continue;
				}

				if ( ncand[i] == MAX_CAND ) {
					ERROR("Too many cell candidates - ");
					ERROR("consider tightening the unit ");
					ERROR("cell tolerances.\n");
				} else {

					double fom;

					fom = fabs(lengths[i] - tlen);

					cand[i][ncand[i]].vec.u = tx;
					cand[i][ncand[i]].vec.v = ty;
					cand[i][ncand[i]].vec.w = tz;
					cand[i][ncand[i]].na = n1;
					cand[i][ncand[i]].nb = n2;
					cand[i][ncand[i]].nc = n3;
					cand[i][ncand[i]].fom = fom;

					ncand[i]++;

				}
			}

		}
		}
		}

	}
	}
	}

	if ( verbose ) {
		STATUS("Candidates: %i %i %i\n", ncand[0], ncand[1], ncand[2]);
	}

	for ( i=0; i<ncand[0]; i++ ) {
	for ( j=0; j<ncand[1]; j++ ) {

		double ang;
		int k;
		float fom1;

		if ( same_vector(cand[0][i], cand[1][j]) ) continue;

		/* Measure the angle between the ith candidate for axis 0
		 * and the jth candidate for axis 1 */
		ang = angle_between(cand[0][i].vec.u, cand[0][i].vec.v,
		                    cand[0][i].vec.w, cand[1][j].vec.u,
		                    cand[1][j].vec.v, cand[1][j].vec.w);

		/* Angle between axes 0 and 1 should be angle 2 */
		if ( fabs(ang - angles[2]) > angtol ) continue;

		fom1 = fabs(ang - angles[2]);

		for ( k=0; k<ncand[2]; k++ ) {

			float fom2, fom3;

			if ( same_vector(cand[1][j], cand[2][k]) ) continue;

			/* Measure the angle between the current candidate for
			 * axis 0 and the kth candidate for axis 2 */
			ang = angle_between(cand[0][i].vec.u, cand[0][i].vec.v,
			                    cand[0][i].vec.w, cand[2][k].vec.u,
			                    cand[2][k].vec.v, cand[2][k].vec.w);

			/* ... it should be angle 1 ... */
			if ( fabs(ang - angles[1]) > angtol ) continue;

			fom2 = fom1 + fabs(ang - angles[1]);

			/* Finally, the angle between the current candidate for
			 * axis 1 and the kth candidate for axis 2 */
			ang = angle_between(cand[1][j].vec.u, cand[1][j].vec.v,
			                    cand[1][j].vec.w, cand[2][k].vec.u,
			                    cand[2][k].vec.v, cand[2][k].vec.w);

			/* ... it should be angle 0 ... */
			if ( fabs(ang - angles[0]) > angtol ) continue;

			/* Unit cell must be right-handed */
			if ( !right_handed_vec(cand[0][i].vec, cand[1][j].vec,
			                       cand[2][k].vec) ) continue;

			fom3 = fom2 + fabs(ang - angles[0]);
			fom3 += LWEIGHT * (cand[0][i].fom + cand[1][j].fom
			                   + cand[2][k].fom);

			if ( fom3 < best_fom ) {
				if ( new_cell != NULL ) free(new_cell);
				new_cell = cell_new_from_reciprocal_axes(
				                 cand[0][i].vec, cand[1][j].vec,
			                         cand[2][k].vec);
				best_fom = fom3;
			}

		}

	}
	}

	free(cand[0]);
	free(cand[1]);
	free(cand[2]);

	cell_free(cell);

	/* Reverse the de-centering transformation */
	if ( new_cell != NULL ) {

		new_cell_trans = cell_transform_inverse(new_cell, uncentering);
		cell_free(new_cell);
		cell_set_lattice_type(new_cell_trans,
		                      cell_get_lattice_type(template_in));
		cell_set_centering(new_cell_trans,
		                   cell_get_centering(template_in));
		cell_set_unique_axis(new_cell_trans,
		                     cell_get_unique_axis(template_in));

		cell_free(template);
		tfn_free(uncentering);

		return new_cell_trans;

	} else {
		cell_free(template);
		tfn_free(uncentering);
		return NULL;
	}
}


UnitCell *match_cell_ab(UnitCell *cell_in, UnitCell *template_in)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	int i;
	double lengths[3];
	int used[3];
	struct rvec real_a, real_b, real_c;
	struct rvec params[3];
	double alen, blen;
	float ltl = 5.0;     /* percent */
	int have_real_a;
	int have_real_b;
	int have_real_c;
	UnitCell *cell;
	UnitCell *template;
	UnitCellTransformation *to_given_cell;
	UnitCell *new_cell;
	UnitCell *new_cell_trans;

	/* "Un-center" the template unit cell to make the comparison easier */
	template = uncenter_cell(template_in, &to_given_cell);

	/* The candidate cell is also uncentered, because it might be centered
	 * if it came from (e.g.) MOSFLM */
	cell = uncenter_cell(cell_in, NULL);

	/* Get the lengths to match */
	if ( cell_get_cartesian(template, &ax, &ay, &az,
	                                  &bx, &by, &bz,
	                                  &cx, &cy, &cz) )
	{
		ERROR("Couldn't get cell for template.\n");
		return NULL;
	}
	alen = modulus(ax, ay, az);
	blen = modulus(bx, by, bz);

	/* Get the lengths from the cell and turn them into anonymous vectors */
	if ( cell_get_cartesian(cell, &ax, &ay, &az,
	                              &bx, &by, &bz,
	                              &cx, &cy, &cz) )
	{
		ERROR("Couldn't get cell.\n");
		return NULL;
	}
	lengths[0] = modulus(ax, ay, az);
	lengths[1] = modulus(bx, by, bz);
	lengths[2] = modulus(cx, cy, cz);
	used[0] = 0;  used[1] = 0;  used[2] = 0;
	params[0].u = ax;  params[0].v = ay;  params[0].w = az;
	params[1].u = bx;  params[1].v = by;  params[1].w = bz;
	params[2].u = cx;  params[2].v = cy;  params[2].w = cz;

	real_a.u = 0.0;  real_a.v = 0.0;  real_a.w = 0.0;
	real_b.u = 0.0;  real_b.v = 0.0;  real_b.w = 0.0;
	real_c.u = 0.0;  real_c.v = 0.0;  real_c.w = 0.0;

	/* Check each vector against a and b */
	have_real_a = 0;
	have_real_b = 0;
	for ( i=0; i<3; i++ ) {
		if ( within_tolerance(lengths[i], alen, ltl)
		     && !used[i] && !have_real_a )
		{
			used[i] = 1;
			memcpy(&real_a, &params[i], sizeof(struct rvec));
			have_real_a = 1;
		}
		if ( within_tolerance(lengths[i], blen, ltl)
		     && !used[i] && !have_real_b )
		{
			used[i] = 1;
			memcpy(&real_b, &params[i], sizeof(struct rvec));
			have_real_b = 1;
		}
	}

	/* Have we matched both a and b? */
	if ( !(have_real_a && have_real_b) ) return NULL;

	/* "c" is "the other one" */
	have_real_c = 0;
	for ( i=0; i<3; i++ ) {
		if ( !used[i] ) {
			memcpy(&real_c, &params[i], sizeof(struct rvec));
			have_real_c = 1;
		}
	}

	if ( !have_real_c ) {
		ERROR("Huh?  Couldn't find the third vector.\n");
		ERROR("Matches: %i %i %i\n", used[0], used[1], used[2]);
		return NULL;
	}

	/* Flip c if not right-handed */
	if ( !right_handed_vec(real_a, real_b, real_c) ) {
		real_c.u = -real_c.u;
		real_c.v = -real_c.v;
		real_c.w = -real_c.w;
	}

	new_cell = cell_new_from_direct_axes(real_a, real_b, real_c);

	 /* Reverse the de-centering transformation */
	new_cell_trans = cell_transform_inverse(new_cell, to_given_cell);
	cell_free(new_cell);
	cell_set_lattice_type(new_cell, cell_get_lattice_type(template_in));
	cell_set_centering(new_cell, cell_get_centering(template_in));
	cell_set_unique_axis(new_cell, cell_get_unique_axis(template_in));

	return new_cell_trans;
}


/* Return sin(theta)/lambda = 1/2d.  Multiply by two if you want 1/d */
double resolution(UnitCell *cell, signed int h, signed int k, signed int l)
{
	double a, b, c, alpha, beta, gamma;

	cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);

	const double Vsq = a*a*b*b*c*c*(1 - cos(alpha)*cos(alpha)
	                                  - cos(beta)*cos(beta)
	                                  - cos(gamma)*cos(gamma)
				          + 2*cos(alpha)*cos(beta)*cos(gamma) );

	const double S11 = b*b*c*c*sin(alpha)*sin(alpha);
	const double S22 = a*a*c*c*sin(beta)*sin(beta);
	const double S33 = a*a*b*b*sin(gamma)*sin(gamma);
	const double S12 = a*b*c*c*(cos(alpha)*cos(beta) - cos(gamma));
	const double S23 = a*a*b*c*(cos(beta)*cos(gamma) - cos(alpha));
	const double S13 = a*b*b*c*(cos(gamma)*cos(alpha) - cos(beta));

	const double brackets = S11*h*h + S22*k*k + S33*l*l
	                         + 2*S12*h*k + 2*S23*k*l + 2*S13*h*l;
	const double oneoverdsq = brackets / Vsq;
	const double oneoverd = sqrt(oneoverdsq);

	return oneoverd / 2;
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
 * load_cell_from_pdb:
 *
 * Loads a unit cell from a PDB file.
 *
 * Returns: a newly allocated %UnitCell.
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

	validate_cell(cell);

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
 * load_cell_from_file:
 *
 * Loads a unit cell from a file of any type (PDB or CrystFEL format)
 *
 * Returns: a newly allocated %UnitCell.
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
 * rotate_cell:
 * @in: A %UnitCell to rotate
 * @omega: Euler angle about +z
 * @phi: Euler angle about +x
 * @rot: Euler angle about new +z
 *
 * Rotate a %UnitCell using Euler angles
 *
 * Returns: a newly allocated rotated copy of @in.
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
 * validate_cell:
 * @cell: A %UnitCell to validate
 *
 * Perform some checks for crystallographic validity @cell, such as that the
 * lattice is a conventional Bravais lattice.
 * Warnings are printied if any of the checks are failed.
 *
 * Returns: true if cell is invalid.
 *
 */
int validate_cell(UnitCell *cell)
{
	int err = 0;
	char cen, ua;

	if ( cell_has_parameters(cell) && !cell_is_sensible(cell) ) {
		ERROR("WARNING: Unit cell parameters are not sensible.\n");
		err = 1;
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
			err = 1;
		}
	}

	return err;
}


/**
 * forbidden_reflection:
 * @cell: A %UnitCell
 * @h: h index to check
 * @k: k index to check
 * @l: l index to check
 *
 * Returns: true if this reflection is forbidden.
 *
 */
int forbidden_reflection(UnitCell *cell,
                         signed int h, signed int k, signed int l)
{
	char cen;

	cen = cell_get_centering(cell);

	/* Reflection conditions here must match the transformation matrices
	 * in uncentering_transformation().  tests/centering_check verifies
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
