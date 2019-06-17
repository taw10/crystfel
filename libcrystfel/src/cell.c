/*
 * cell.c
 *
 * A class representing a unit cell
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2012,2014,2017 Thomas White <taw@physics.org>
 *   2010                Richard Kirian
 *   2012                Lorenzo Galli
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

#include "cell.h"
#include "utils.h"
#include "image.h"
#include "integer_matrix.h"
#include "rational.h"


/**
 * \file cell.h
 */


typedef enum {
	CELL_REP_CRYST,
	CELL_REP_CART,
	CELL_REP_RECIP
} CellRepresentation;

struct _unitcell {

	CellRepresentation rep;

	int have_parameters;

	/* Crystallographic representation */
	double a;	/* m */
	double b;	/* m */
	double c;	/* m */
	double alpha;	/* Radians */
	double beta;	/* Radians */
	double gamma;	/* Radians */

	/* Cartesian representation */
	double ax;	double bx;	double cx;
	double ay;	double by;	double cy;
	double az;	double bz;	double cz;

	/* Cartesian representation of reciprocal axes */
	double axs;	double bxs;	double cxs;
	double ays;	double bys;	double cys;
	double azs;	double bzs;	double czs;

	LatticeType  lattice_type;
	char         centering;
	char         unique_axis;
};

typedef enum {
	CMASK_P = 1<<0,
	CMASK_A = 1<<1,
	CMASK_B = 1<<2,
	CMASK_C = 1<<3,
	CMASK_I = 1<<4,
	CMASK_F = 1<<5,
	CMASK_H = 1<<6,
	CMASK_R = 1<<7
} CenteringMask;

#define CMASK_ALL (CMASK_P | CMASK_A | CMASK_B | CMASK_C | CMASK_I \
                     | CMASK_F | CMASK_H | CMASK_R)


/************************** Setters and Constructors **************************/


/**
 * Create a new UnitCell.
 *
 * \returns the new unit cell, or NULL on failure.
 *
 */
UnitCell *cell_new()
{
	UnitCell *cell;

	cell = malloc(sizeof(UnitCell));
	if ( cell == NULL ) return NULL;

	cell->a = 1.0;
	cell->b = 1.0;
	cell->c = 1.0;
	cell->alpha = 0.0;
	cell->beta = 0.0;
	cell->gamma = 0.0;

	cell->rep = CELL_REP_CRYST;

	cell->lattice_type = L_TRICLINIC;
	cell->centering = 'P';
	cell->unique_axis = '?';
	cell->have_parameters = 0;

	return cell;
}


/**
 * \param cell: A %UnitCell to free.
 *
 * Frees a %UnitCell, and all internal resources concerning that cell.
 *
 */
void cell_free(UnitCell *cell)
{
	if ( cell == NULL ) return;
	free(cell);
}


/**
 * \param cell: A %UnitCell
 *
 * \returns True if cell has its parameters specified.
 *
 */
int cell_has_parameters(UnitCell *cell)
{
	if ( cell == NULL ) return 0;
	return cell->have_parameters;
}


void cell_set_parameters(UnitCell *cell, double a, double b, double c,
                         double alpha, double beta, double gamma)
{
	if ( cell == NULL ) return;

	cell->a = a;
	cell->b = b;
	cell->c = c;
	cell->alpha = alpha;
	cell->beta = beta;
	cell->gamma = gamma;

	cell->rep = CELL_REP_CRYST;
	cell->have_parameters = 1;
}


void cell_set_cartesian(UnitCell *cell,
			double ax, double ay, double az,
			double bx, double by, double bz,
			double cx, double cy, double cz)
{
	if ( cell == NULL ) return;

	cell->ax = ax;  cell->ay = ay;  cell->az = az;
	cell->bx = bx;  cell->by = by;  cell->bz = bz;
	cell->cx = cx;  cell->cy = cy;  cell->cz = cz;

	cell->rep = CELL_REP_CART;
	cell->have_parameters = 1;
}


UnitCell *cell_new_from_parameters(double a, double b, double c,
                                   double alpha, double beta, double gamma)
{
	UnitCell *cell;

	cell = cell_new();
	if ( cell == NULL ) return NULL;

	cell_set_parameters(cell, a, b, c, alpha, beta, gamma);

	return cell;
}


UnitCell *cell_new_from_reciprocal_axes(struct rvec as, struct rvec bs,
                                        struct rvec cs)
{
	UnitCell *cell;

	cell = cell_new();
	if ( cell == NULL ) return NULL;

	cell->axs = as.u;  cell->ays = as.v;  cell->azs = as.w;
	cell->bxs = bs.u;  cell->bys = bs.v;  cell->bzs = bs.w;
	cell->cxs = cs.u;  cell->cys = cs.v;  cell->czs = cs.w;

	cell->rep = CELL_REP_RECIP;
	cell->have_parameters = 1;

	return cell;
}


UnitCell *cell_new_from_direct_axes(struct rvec a, struct rvec b, struct rvec c)
{
	UnitCell *cell;

	cell = cell_new();
	if ( cell == NULL ) return NULL;

	cell->ax = a.u;  cell->ay = a.v;  cell->az = a.w;
	cell->bx = b.u;  cell->by = b.v;  cell->bz = b.w;
	cell->cx = c.u;  cell->cy = c.v;  cell->cz = c.w;

	cell->rep = CELL_REP_CART;
	cell->have_parameters = 1;

	return cell;
}


UnitCell *cell_new_from_cell(const UnitCell *orig)
{
	UnitCell *new;
	new = cell_new();
	*new = *orig;
	return new;
}


void cell_set_reciprocal(UnitCell *cell,
                        double asx, double asy, double asz,
                        double bsx, double bsy, double bsz,
                        double csx, double csy, double csz)
{
	if ( cell == NULL ) return;

	cell->axs = asx;  cell->ays = asy;  cell->azs = asz;
	cell->bxs = bsx;  cell->bys = bsy;  cell->bzs = bsz;
	cell->cxs = csx;  cell->cys = csy;  cell->czs = csz;

	cell->rep = CELL_REP_RECIP;
	cell->have_parameters = 1;
}


void cell_set_centering(UnitCell *cell, char centering)
{
	cell->centering = centering;
}


void cell_set_lattice_type(UnitCell *cell, LatticeType lattice_type)
{
	cell->lattice_type = lattice_type;
}


void cell_set_unique_axis(UnitCell *cell, char unique_axis)
{
	cell->unique_axis = unique_axis;
}


/************************* Getter helper functions ****************************/

static int cell_crystallographic_to_cartesian(const UnitCell *cell,
                                              double *ax, double *ay, double *az,
                                              double *bx, double *by, double *bz,
                                              double *cx, double *cy, double *cz)
{
	double tmp, V, cosalphastar, cstar;

	if ( !cell->have_parameters ) {
		ERROR("Unit cell has unspecified parameters.\n");
		return 1;
	}

	/* Firstly: Get a in terms of x, y and z
	 * +a (cryst) is defined to lie along +x (cart) */
	*ax = cell->a;
	*ay = 0.0;
	*az = 0.0;

	/* b in terms of x, y and z
	 * b (cryst) is defined to lie in the xy (cart) plane */
	*bx = cell->b*cos(cell->gamma);
	*by = cell->b*sin(cell->gamma);
	*bz = 0.0;

	tmp = cos(cell->alpha)*cos(cell->alpha)
		+ cos(cell->beta)*cos(cell->beta)
		+ cos(cell->gamma)*cos(cell->gamma)
		- 2.0*cos(cell->alpha)*cos(cell->beta)*cos(cell->gamma);
	V = cell->a * cell->b * cell->c * sqrt(1.0 - tmp);

	cosalphastar = cos(cell->beta)*cos(cell->gamma) - cos(cell->alpha);
	cosalphastar /= sin(cell->beta)*sin(cell->gamma);

	cstar = (cell->a * cell->b * sin(cell->gamma))/V;

	/* c in terms of x, y and z */
	*cx = cell->c*cos(cell->beta);
	*cy = -cell->c*sin(cell->beta)*cosalphastar;
	*cz = 1.0/cstar;

	return 0;
}


/* Why yes, I do enjoy long argument lists...! */
static int cell_invert(double ax, double ay, double az,
                       double bx, double by, double bz,
                       double cx, double cy, double cz,
                       double *asx, double *asy, double *asz,
                       double *bsx, double *bsy, double *bsz,
                       double *csx, double *csy, double *csz)
{
	int s;
	gsl_matrix *m;
	gsl_matrix *inv;
	gsl_permutation *perm;

	m = gsl_matrix_alloc(3, 3);
	if ( m == NULL ) {
		ERROR("Couldn't allocate memory for matrix\n");
		return 1;
	}
	gsl_matrix_set(m, 0, 0, ax);
	gsl_matrix_set(m, 1, 0, ay);
	gsl_matrix_set(m, 2, 0, az);
	gsl_matrix_set(m, 0, 1, bx);
	gsl_matrix_set(m, 1, 1, by);
	gsl_matrix_set(m, 2, 1, bz);
	gsl_matrix_set(m, 0, 2, cx);
	gsl_matrix_set(m, 1, 2, cy);
	gsl_matrix_set(m, 2, 2, cz);

	/* Invert */
	perm = gsl_permutation_alloc(m->size1);
	if ( perm == NULL ) {
		ERROR("Couldn't allocate permutation\n");
		gsl_matrix_free(m);
		return 1;
	}
	inv = gsl_matrix_alloc(m->size1, m->size2);
	if ( inv == NULL ) {
		ERROR("Couldn't allocate inverse\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	if ( gsl_linalg_LU_decomp(m, perm, &s) ) {
		ERROR("Couldn't decompose matrix\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	if ( gsl_linalg_LU_invert(m, perm, inv)  ) {
		ERROR("Couldn't invert cell matrix:\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	gsl_permutation_free(perm);
	gsl_matrix_free(m);

	/* Transpose */
	gsl_matrix_transpose(inv);

	*asx = gsl_matrix_get(inv, 0, 0);
	*asy = gsl_matrix_get(inv, 1, 0);
	*asz = gsl_matrix_get(inv, 2, 0);
	*bsx = gsl_matrix_get(inv, 0, 1);
	*bsy = gsl_matrix_get(inv, 1, 1);
	*bsz = gsl_matrix_get(inv, 2, 1);
	*csx = gsl_matrix_get(inv, 0, 2);
	*csy = gsl_matrix_get(inv, 1, 2);
	*csz = gsl_matrix_get(inv, 2, 2);

	gsl_matrix_free(inv);

	return 0;
}


/********************************** Getters ***********************************/

int cell_get_parameters(const UnitCell *cell, double *a, double *b, double *c,
                        double *alpha, double *beta, double *gamma)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	if ( cell == NULL ) return 1;

	if ( !cell->have_parameters ) {
		ERROR("Unit cell has unspecified parameters.\n");
		return 1;
	}

	switch ( cell->rep ) {

		case CELL_REP_CRYST:
		/* Direct response */
		*a = cell->a;
		*b = cell->b;
		*c = cell->c;
		*alpha = cell->alpha;
		*beta = cell->beta;
		*gamma = cell->gamma;
		return 0;

		case CELL_REP_CART:
		/* Convert cartesian -> crystallographic */
		*a = modulus(cell->ax, cell->ay, cell->az);
		*b = modulus(cell->bx, cell->by, cell->bz);
		*c = modulus(cell->cx, cell->cy, cell->cz);

		*alpha = angle_between(cell->bx, cell->by, cell->bz,
		                       cell->cx, cell->cy, cell->cz);
		*beta = angle_between(cell->ax, cell->ay, cell->az,
		                      cell->cx, cell->cy, cell->cz);
		*gamma = angle_between(cell->ax, cell->ay, cell->az,
		                       cell->bx, cell->by, cell->bz);
		return 0;

		case CELL_REP_RECIP:
		/* Convert reciprocal -> crystallographic.
                 * Start by converting reciprocal -> cartesian */
		if ( cell_invert(cell->axs, cell->ays, cell->azs,
		                 cell->bxs, cell->bys, cell->bzs,
		                         cell->cxs, cell->cys, cell->czs,
		                         &ax, &ay, &az,
		                         &bx, &by, &bz,
		                         &cx, &cy, &cz) ) return 1;

		/* Now convert cartesian -> crystallographic */
		*a = modulus(ax, ay, az);
		*b = modulus(bx, by, bz);
		*c = modulus(cx, cy, cz);

		*alpha = angle_between(bx, by, bz, cx, cy, cz);
		*beta = angle_between(ax, ay, az, cx, cy, cz);
		*gamma = angle_between(ax, ay, az, bx, by, bz);
		return 0;
	}

	return 1;
}


int cell_get_cartesian(const UnitCell *cell,
                       double *ax, double *ay, double *az,
                       double *bx, double *by, double *bz,
                       double *cx, double *cy, double *cz)
{
	if ( cell == NULL ) return 1;

	if ( !cell->have_parameters ) {
		ERROR("Unit cell has unspecified parameters.\n");
		return 1;
	}

	switch ( cell->rep ) {

		case CELL_REP_CRYST:
		/* Convert crystallographic -> cartesian. */
		return cell_crystallographic_to_cartesian(cell,
		                                          ax, ay, az,
		                                          bx, by, bz,
		                                          cx, cy, cz);

		case CELL_REP_CART:
		/* Direct response */
		*ax = cell->ax;  *ay = cell->ay;  *az = cell->az;
		*bx = cell->bx;  *by = cell->by;  *bz = cell->bz;
		*cx = cell->cx;  *cy = cell->cy;  *cz = cell->cz;
		return 0;

		case CELL_REP_RECIP:
		/* Convert reciprocal -> cartesian */
		return cell_invert(cell->axs, cell->ays, cell->azs,
		                   cell->bxs, cell->bys, cell->bzs,
		                   cell->cxs, cell->cys, cell->czs,
		                   ax, ay, az, bx, by, bz, cx, cy, cz);

	}

	return 1;
}


int cell_get_reciprocal(const UnitCell *cell,
                        double *asx, double *asy, double *asz,
                        double *bsx, double *bsy, double *bsz,
                        double *csx, double *csy, double *csz)
{
	int r;
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	if ( cell == NULL ) return 1;

	if ( !cell->have_parameters ) {
		ERROR("Unit cell has unspecified parameters.\n");
		return 1;
	}

	switch ( cell->rep ) {

		case CELL_REP_CRYST:
		/* Convert crystallographic -> reciprocal */
		r = cell_crystallographic_to_cartesian(cell,
		                                       &ax, &ay, &az,
		                                       &bx, &by, &bz,
		                                       &cx, &cy, &cz);
		if ( r ) return r;
		return cell_invert(ax, ay, az,bx, by, bz, cx, cy, cz,
		                   asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

		case CELL_REP_CART:
		/* Convert cartesian -> reciprocal */
		cell_invert(cell->ax, cell->ay, cell->az,
		            cell->bx, cell->by, cell->bz,
		            cell->cx, cell->cy, cell->cz,
		            asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);
		return 0;

		case CELL_REP_RECIP:
		/* Direct response */
		*asx = cell->axs;  *asy = cell->ays;  *asz = cell->azs;
		*bsx = cell->bxs;  *bsy = cell->bys;  *bsz = cell->bzs;
		*csx = cell->cxs;  *csy = cell->cys;  *csz = cell->czs;
		return 0;

	}

	return 1;
}


char cell_get_centering(UnitCell *cell)
{
	return cell->centering;
}


LatticeType cell_get_lattice_type(UnitCell *cell)
{
	return cell->lattice_type;
}


char cell_get_unique_axis(UnitCell *cell)
{
	return cell->unique_axis;
}


const char *cell_rep(UnitCell *cell)
{
	switch ( cell->rep ) {

		case CELL_REP_CRYST:
		return "crystallographic, direct space";

		case CELL_REP_CART:
		return "cartesian, direct space";

		case CELL_REP_RECIP:
		return "cartesian, reciprocal space";

	}

	return "unknown";
}


UnitCell *cell_transform_gsl_direct(UnitCell *in, gsl_matrix *m)
{
	gsl_matrix *c;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	gsl_matrix *res;
	UnitCell *out;

	cell_get_cartesian(in, &asx, &asy, &asz, &bsx, &bsy,
	                       &bsz, &csx, &csy, &csz);

	c = gsl_matrix_alloc(3, 3);
	gsl_matrix_set(c, 0, 0, asx);
	gsl_matrix_set(c, 1, 0, asy);
	gsl_matrix_set(c, 2, 0, asz);
	gsl_matrix_set(c, 0, 1, bsx);
	gsl_matrix_set(c, 1, 1, bsy);
	gsl_matrix_set(c, 2, 1, bsz);
	gsl_matrix_set(c, 0, 2, csx);
	gsl_matrix_set(c, 1, 2, csy);
	gsl_matrix_set(c, 2, 2, csz);

	res = gsl_matrix_calloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, c, m, 0.0, res);

	out = cell_new_from_cell(in);
	cell_set_cartesian(out, gsl_matrix_get(res, 0, 0),
	                        gsl_matrix_get(res, 1, 0),
	                        gsl_matrix_get(res, 2, 0),
	                        gsl_matrix_get(res, 0, 1),
	                        gsl_matrix_get(res, 1, 1),
	                        gsl_matrix_get(res, 2, 1),
	                        gsl_matrix_get(res, 0, 2),
	                        gsl_matrix_get(res, 1, 2),
	                        gsl_matrix_get(res, 2, 2));

	gsl_matrix_free(res);
	gsl_matrix_free(c);
	return out;
}


static int centering_has_point(char cen, Rational *p)
{
	/* First, put the point into the range 0..1 */
	while ( rtnl_cmp(p[0], rtnl_zero()) < 0 ) p[0] = rtnl_add(p[0], rtnl(1, 1));
	while ( rtnl_cmp(p[1], rtnl_zero()) < 0 ) p[1] = rtnl_add(p[1], rtnl(1, 1));
	while ( rtnl_cmp(p[2], rtnl_zero()) < 0 ) p[2] = rtnl_add(p[2], rtnl(1, 1));
	while ( rtnl_cmp(p[0], rtnl(1, 1)) >= 0 ) p[0] = rtnl_sub(p[0], rtnl(1, 1));
	while ( rtnl_cmp(p[1], rtnl(1, 1)) >= 0 ) p[1] = rtnl_sub(p[1], rtnl(1, 1));
	while ( rtnl_cmp(p[2], rtnl(1, 1)) >= 0 ) p[2] = rtnl_sub(p[2], rtnl(1, 1));

	/* 0,0,0 is present in all centerings */
	if ( (rtnl_cmp(p[0], rtnl_zero()) == 0)
	  && (rtnl_cmp(p[1], rtnl_zero()) == 0)
	  && (rtnl_cmp(p[2], rtnl_zero()) == 0) ) return 1;

	/* Only I has 1/2 , 1/2, 1/2 */
	if ( (rtnl_cmp(p[0], rtnl(1,2)) == 0)
	  && (rtnl_cmp(p[1], rtnl(1,2)) == 0)
	  && (rtnl_cmp(p[2], rtnl(1,2)) == 0)
	  && (cen == 'I') ) return 1;

	/* A or F has 0 , 1/2, 1/2 */
	if ( (rtnl_cmp(p[0], rtnl_zero()) == 0)
	  && (rtnl_cmp(p[1], rtnl(1,2)) == 0)
	  && (rtnl_cmp(p[2], rtnl(1,2)) == 0)
	  && ((cen == 'A') || (cen == 'F')) ) return 1;

	/* B or F has 1/2 , 0 , 1/2 */
	if ( (rtnl_cmp(p[0], rtnl(1,2)) == 0)
	  && (rtnl_cmp(p[1], rtnl_zero()) == 0)
	  && (rtnl_cmp(p[2], rtnl(1,2)) == 0)
	  && ((cen == 'B') || (cen == 'F')) ) return 1;

	/* C or F has 1/2 , 1/2 , 0 */
	if ( (rtnl_cmp(p[0], rtnl(1,2)) == 0)
	  && (rtnl_cmp(p[1], rtnl(1,2)) == 0)
	  && (rtnl_cmp(p[2], rtnl_zero()) == 0)
	  && ((cen == 'C') || (cen == 'F')) ) return 1;

	/* H has 2/3 , 1/3 , 1/3 */
	if ( (rtnl_cmp(p[0], rtnl(2,3)) == 0)
	  && (rtnl_cmp(p[1], rtnl(1,3)) == 0)
	  && (rtnl_cmp(p[2], rtnl(1,3)) == 0)
	  && (cen == 'H') ) return 1;

	/* H has 1/3 , 2/3 , 2/3 */
	if ( (rtnl_cmp(p[0], rtnl(1,3)) == 0)
	  && (rtnl_cmp(p[1], rtnl(2,3)) == 0)
	  && (rtnl_cmp(p[2], rtnl(2,3)) == 0)
	  && (cen == 'H') ) return 1;

	return 0;
}


static void maybe_eliminate(CenteringMask c, CenteringMask *cmask, Rational *nc,
                            char cen)
{
	/* Skip test if this centering isn't even a candidate */
	if ( !(*cmask & c) ) return;

	if ( !centering_has_point(cen, nc) ) {
		*cmask |= c;
		*cmask ^= c;
	}
}


/* Check if the point x,y,z in the original cell matches any lattice point
 * in the transformed cell */
static void check_point_fwd(RationalMatrix *P, CenteringMask *cmask,
                            Rational x, Rational y, Rational z)
{
	Rational c[3] = {x, y, z};
	Rational nc[3];

	/* Transform the lattice point */
	transform_fractional_coords_rtnl(P, c, nc);

	/* Eliminate any centerings which don't include the transformed point */
	maybe_eliminate(CMASK_P, cmask, nc, 'P');
	maybe_eliminate(CMASK_R, cmask, nc, 'R');
	maybe_eliminate(CMASK_A, cmask, nc, 'A');
	maybe_eliminate(CMASK_B, cmask, nc, 'B');
	maybe_eliminate(CMASK_C, cmask, nc, 'C');
	maybe_eliminate(CMASK_I, cmask, nc, 'I');
	maybe_eliminate(CMASK_F, cmask, nc, 'F');
	maybe_eliminate(CMASK_H, cmask, nc, 'H');
}


/* Check if the point x,y,z in the transformed cell matches any lattice point
 * in the original cell.  If not, eliminate "exclude" from "*mask". */
static void check_point_bwd(RationalMatrix *P, CenteringMask *mask,
                            char cen, CenteringMask exclude,
                            Rational x, Rational y, Rational z)
{
	Rational nc[3];
	Rational c[3] = {x, y, z};

	transform_fractional_coords_rtnl_inverse(P, c, nc);

	if ( !centering_has_point(cen, nc) ) {
		*mask |= exclude;
		*mask ^= exclude;  /* Unset bits */
	}
}


static char cmask_decode(CenteringMask mask)
{
	char res[32];

	res[0] = '\0';

	if ( mask & CMASK_H ) strcat(res, "H");
	if ( mask & CMASK_F ) strcat(res, "F");
	if ( mask & CMASK_I ) strcat(res, "I");
	if ( mask & CMASK_A ) strcat(res, "A");
	if ( mask & CMASK_B ) strcat(res, "B");
	if ( mask & CMASK_C ) strcat(res, "C");
	if ( mask & CMASK_P ) strcat(res, "P");
	if ( mask & CMASK_R ) strcat(res, "R");

	if ( strlen(res) == 0 ) return '?';
	return res[0];
}


static char determine_centering(RationalMatrix *P, char cen)
{
	CenteringMask cmask = CMASK_ALL;

	/* Check whether the current centering can provide all the lattice
	 * points for the transformed cell.  Eliminate any centerings for which
	 * it can't. */
	check_point_bwd(P, &cmask, cen, CMASK_A | CMASK_F, rtnl_zero(), rtnl(1,2), rtnl(1,2));
	check_point_bwd(P, &cmask, cen, CMASK_B | CMASK_F, rtnl(1,2), rtnl_zero(), rtnl(1,2));
	check_point_bwd(P, &cmask, cen, CMASK_C | CMASK_F, rtnl(1,2), rtnl(1,2), rtnl_zero());
	check_point_bwd(P, &cmask, cen, CMASK_I, rtnl(1,2), rtnl(1,2), rtnl(1,2));
	check_point_bwd(P, &cmask, cen, CMASK_H, rtnl(2,3), rtnl(1,3), rtnl(1,3));
	check_point_bwd(P, &cmask, cen, CMASK_H, rtnl(1,3), rtnl(2,3), rtnl(2,3));
	check_point_bwd(P, &cmask, cen, CMASK_ALL, rtnl(1,1), rtnl(1,1), rtnl(1,1));

	/* Check whether the current centering's lattice points will all
	 * coincide with lattice points in the new centering.  Eliminate any
	 * centerings for which they don't (they give "excess lattice points"). */
	switch ( cen ) {

		case 'P' :
		case 'R' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		break;

		case 'A' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		check_point_fwd(P, &cmask, rtnl_zero(), rtnl(1,2), rtnl(1,2));
		break;

		case 'B' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		check_point_fwd(P, &cmask, rtnl(1,2), rtnl_zero(), rtnl(1,2));
		break;

		case 'C' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		check_point_fwd(P, &cmask, rtnl(1,2), rtnl(1,2), rtnl_zero());
		break;

		case 'I' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		check_point_fwd(P, &cmask, rtnl(1,2), rtnl(1,2), rtnl(1,2));
		break;

		case 'F' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		check_point_fwd(P, &cmask, rtnl_zero(), rtnl(1,2), rtnl(1,2));
		check_point_fwd(P, &cmask, rtnl(1,2), rtnl_zero(), rtnl(1,2));
		check_point_fwd(P, &cmask, rtnl(1,2), rtnl(1,2), rtnl_zero());
		break;

		case 'H' :
		check_point_fwd(P, &cmask, rtnl(1,1), rtnl(1,1), rtnl(1,1));
		check_point_fwd(P, &cmask, rtnl(2,3), rtnl(1,3), rtnl(1,3));
		check_point_fwd(P, &cmask, rtnl(1,3), rtnl(2,3), rtnl(2,3));
		break;

	}

	return cmask_decode(cmask);
}


/**
 * \param cell: A %UnitCell.
 * \param m: A %RationalMatrix.
 *
 * Applies \p m to \p cell.
 *
 * This function will determine the centering of the resulting unit cell,
 * producing '?' if any lattice points cannot be accounted for.  Note that if
 * there are 'excess' lattice points in the transformed cell, the centering
 * will still be '?' even if the lattice points for another centering are
 * all present.
 *
 * The lattice type will be set to triclinic, and the unique axis to '?'.
 *
 * \returns Transformed copy of \p cell.
 *
 */
UnitCell *cell_transform_rational(UnitCell *cell, RationalMatrix *m)
{
	UnitCell *out;
	gsl_matrix *tm;
	char ncen;
	int i, j;
	Rational det;

	if ( m == NULL ) return NULL;

	det = rtnl_mtx_det(m);
	if ( rtnl_cmp(det, rtnl_zero()) == 0 ) return NULL;

	tm = gsl_matrix_alloc(3,3);
	if ( tm == NULL ) {
		return NULL;
	}

	for ( i=0; i<3; i++ ) {
		for ( j=0; j<3; j++ ) {
			gsl_matrix_set(tm, i, j,
			               rtnl_as_double(rtnl_mtx_get(m, i, j)));
		}
	}

	out = cell_transform_gsl_direct(cell, tm);
	gsl_matrix_free(tm);

	ncen = determine_centering(m, cell_get_centering(cell));
	cell_set_centering(out, ncen);

	/* FIXME: Update unique axis, lattice type */
	cell_set_lattice_type(out, L_TRICLINIC);
	cell_set_unique_axis(out, '?');

	return out;
}


/**
 * \param cell: A %UnitCell.
 * \param m: An %IntegerMatrix.
 *
 * Applies \p m to \p cell.
 *
 * This is just a convenience function which turns \p m into a %RationalMatrix
 * and then calls cell_transform_rational().  See the documentation for that
 * function for some important information.
 *
 * \returns Transformed copy of \p cell.
 *
 */
UnitCell *cell_transform_intmat(UnitCell *cell, IntegerMatrix *m)
{
	UnitCell *ans;
	RationalMatrix *mtx = rtnl_mtx_from_intmat(m);
	ans = cell_transform_rational(cell, mtx);
	rtnl_mtx_free(mtx);
	return ans;
}


/**
 * \param cell: A %UnitCell.
 * \param m: A %RationalMatrix
 *
 * Applies the inverse of \p m to \p cell.
 *
 * \returns Transformed copy of \p cell.
 *
 */
UnitCell *cell_transform_rational_inverse(UnitCell *cell, RationalMatrix *m)
{
	UnitCell *out;
	gsl_matrix *tm;
	gsl_matrix *inv;
	gsl_permutation *perm;
	int s;
	int i, j;

	if ( m == NULL ) return NULL;

	tm = gsl_matrix_alloc(3,3);
	if ( tm == NULL ) {
		return NULL;
	}

	for ( i=0; i<3; i++ ) {
		for ( j=0; j<3; j++ ) {
			gsl_matrix_set(tm, i, j,
			               rtnl_as_double(rtnl_mtx_get(m, i, j)));
		}
	}

	perm = gsl_permutation_alloc(3);
	if ( perm == NULL ) {
		ERROR("Couldn't allocate permutation\n");
		return NULL;
	}
	inv = gsl_matrix_alloc(3, 3);
	if ( inv == NULL ) {
		ERROR("Couldn't allocate inverse\n");
		gsl_permutation_free(perm);
		return NULL;
	}
	if ( gsl_linalg_LU_decomp(tm, perm, &s) ) {
		ERROR("Couldn't decompose matrix\n");
		gsl_permutation_free(perm);
		return NULL;
	}
	if ( gsl_linalg_LU_invert(tm, perm, inv)  ) {
		ERROR("Couldn't invert transformation matrix\n");
		gsl_permutation_free(perm);
		return NULL;
	}
	gsl_permutation_free(perm);

	out = cell_transform_gsl_direct(cell, inv);

	/* FIXME: Update centering, unique axis, lattice type */

	gsl_matrix_free(tm);
	gsl_matrix_free(inv);

	return out;
}


/**
 * \param cell: A %UnitCell.
 * \param m: An %IntegerMatrix
 *
 * Applies the inverse of \p m to \p cell.
 *
 * \returns Transformed copy of \p cell.
 *
 */
UnitCell *cell_transform_intmat_inverse(UnitCell *cell, IntegerMatrix *m)
{
	UnitCell *ans;
	RationalMatrix *mtx = rtnl_mtx_from_intmat(m);
	ans = cell_transform_rational_inverse(cell, mtx);
	rtnl_mtx_free(mtx);
	return ans;
}
