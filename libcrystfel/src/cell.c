/*
 * cell.c
 *
 * A class representing a unit cell
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *   2010      Richard Kirian
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


/**
 * SECTION:unitcell
 * @short_description: Unit cell
 * @title: UnitCell
 * @section_id:
 * @see_also:
 * @include: "cell.h"
 * @Image:
 *
 * This structure represents a unit cell.
 */


typedef enum {
	CELL_REP_CRYST,
	CELL_REP_CART,
	CELL_REP_RECIP
} CellRepresentation;

struct _unitcell {

	CellRepresentation rep;

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

	char        *pointgroup;
	LatticeType  lattice_type;
	char         centering;
	char         unique_axis;
};


/************************** Setters and Constructors **************************/


/**
 * cell_new:
 *
 * Create a new %UnitCell.
 *
 * Returns: the new unit cell, or NULL on failure.
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
	cell->alpha = M_PI_2;
	cell->beta = M_PI_2;
	cell->gamma = M_PI_2;

	cell->rep = CELL_REP_CRYST;

	cell->pointgroup = strdup("1");
	cell->lattice_type = L_TRICLINIC;
	cell->centering = 'P';
	cell->unique_axis = 'c';

	return cell;
}


/**
 * cell_free:
 * @cell: A %UnitCell to free.
 *
 * Frees a %UnitCell, and all internal resources concerning that cell.
 *
 */
void cell_free(UnitCell *cell)
{
	if ( cell == NULL ) return;
	free(cell->pointgroup);
	free(cell);
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
}


void cell_set_cartesian_a(UnitCell *cell, double ax, double ay, double az)
{
	if ( cell == NULL ) return;
	cell->ax = ax;  cell->ay = ay;  cell->az = az;
	cell->rep = CELL_REP_CART;
}


void cell_set_cartesian_b(UnitCell *cell, double bx, double by, double bz)
{
	if ( cell == NULL ) return;
	cell->bx = bx;  cell->by = by;  cell->bz = bz;
	cell->rep = CELL_REP_CART;
}


void cell_set_cartesian_c(UnitCell *cell, double cx, double cy, double cz)
{
	if ( cell == NULL ) return;
	cell->cx = cx;  cell->cy = cy;  cell->cz = cz;
	cell->rep = CELL_REP_CART;
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

	return cell;
}


UnitCell *cell_new_from_cell(UnitCell *orig)
{
	UnitCell *new;
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	new = cell_new();

	cell_get_cartesian(orig, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	cell_set_cartesian(new, ax, ay, az, bx, by, bz, cx, cy, cz);
	cell_set_pointgroup(new, orig->pointgroup);

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
}


void cell_set_pointgroup(UnitCell *cell, const char *sym)
{
	free(cell->pointgroup);
	cell->pointgroup = strdup(sym);
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

static int cell_crystallographic_to_cartesian(UnitCell *cell,
                                             double *ax, double *ay, double *az,
                                             double *bx, double *by, double *bz,
                                             double *cx, double *cy, double *cz)
{
	double tmp, V, cosalphastar, cstar;

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
	gsl_matrix_set(m, 0, 1, bx);
	gsl_matrix_set(m, 0, 2, cx);
	gsl_matrix_set(m, 1, 0, ay);
	gsl_matrix_set(m, 1, 1, by);
	gsl_matrix_set(m, 1, 2, cy);
	gsl_matrix_set(m, 2, 0, az);
	gsl_matrix_set(m, 2, 1, bz);
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
		ERROR("Couldn't invert matrix\n");
		gsl_matrix_free(m);
		gsl_permutation_free(perm);
		return 1;
	}
	gsl_permutation_free(perm);
	gsl_matrix_free(m);

	/* Transpose */
	gsl_matrix_transpose(inv);

	*asx = gsl_matrix_get(inv, 0, 0);
	*bsx = gsl_matrix_get(inv, 0, 1);
	*csx = gsl_matrix_get(inv, 0, 2);
	*asy = gsl_matrix_get(inv, 1, 0);
	*bsy = gsl_matrix_get(inv, 1, 1);
	*csy = gsl_matrix_get(inv, 1, 2);
	*asz = gsl_matrix_get(inv, 2, 0);
	*bsz = gsl_matrix_get(inv, 2, 1);
	*csz = gsl_matrix_get(inv, 2, 2);

	gsl_matrix_free(inv);

	return 0;
}


/********************************** Getters ***********************************/

int cell_get_parameters(UnitCell *cell, double *a, double *b, double *c,
                        double *alpha, double *beta, double *gamma)
{
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	if ( cell == NULL ) return 1;

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
		cell_invert(cell->axs, cell->ays, cell->azs,
		            cell->bxs, cell->bys, cell->bzs,
		            cell->cxs, cell->cys, cell->czs,
		            &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

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


int cell_get_cartesian(UnitCell *cell,
                       double *ax, double *ay, double *az,
                       double *bx, double *by, double *bz,
                       double *cx, double *cy, double *cz)
{
	if ( cell == NULL ) return 1;

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


int cell_get_reciprocal(UnitCell *cell,
                        double *asx, double *asy, double *asz,
                        double *bsx, double *bsy, double *bsz,
                        double *csx, double *csy, double *csz)
{
	int r;
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	if ( cell == NULL ) return 1;

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


const char *cell_get_pointgroup(UnitCell *cell)
{
	return cell->pointgroup;
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


static UnitCellTransformation *inverse_transformation(UnitCellTransformation *t)
{
	/* FIXME: Implementation */
	return NULL;
}


/**
 * cell_transform:
 * @cell: A %UnitCell.
 * @t: A %UnitCellTransformation.
 *
 * Applies @t to @cell.
 *
 * Returns: Transformed copy of @cell.
 *
 */
UnitCell *cell_transform(UnitCell *cell, UnitCellTransformation *t)
{
	UnitCell *out;

	if ( t == NULL ) return NULL;

	out = cell_new();
	if ( out == NULL ) return NULL;

	/* FIXME: Implementation */

	return out;
}


/**
 * cell_transform_inverse:
 * @cell: A %UnitCell.
 * @t: A %UnitCellTransformation.
 *
 * Applies the inverse of @t to @cell.
 *
 * Returns: Transformed copy of @cell.
 *
 */
UnitCell *cell_transform_inverse(UnitCell *cell, UnitCellTransformation *t)
{
	return cell_transform(cell, inverse_transformation(t));
}


/**
 * tfn_identity:
 *
 * Returns: A %UnitCellTransformation corresponding to an identity operation.
 *
 */
static UnitCellTransformation *tfn_identity()
{
}


/**
 * tfn_combine:
 * @t: A %UnitCellTransformation
 * @na: Pointer to three doubles representing naa, nab, nac
 * @nb: Pointer to three doubles representing nba, nbb, nbc
 * @nc: Pointer to three doubles representing nca, ncb, ncc
 *
 * Updates @t such that it represents its previous transformation followed by
 * a new transformation, corresponding to letting a = naa*a + nab*b + nac*c.
 * Likewise, a = nba*a + nbb*b + nbc*c and c = nca*a + ncb*b + ncc*c.
 *
 */
static void tfn_combine(UnitCellTransformation *t,
                        double *na, double *nb, double *nc)
{
}


static double *v(double a, double b, double c)
{
	double *vec = malloc(3*sizeof(double));
	if ( vec == NULL ) return NULL;
	vec[0] = a;  vec[1] = b;  vec[2] = c;
	return vec;
}


void cell_transformation_print(UnitCellTransformation *t)
{

}
