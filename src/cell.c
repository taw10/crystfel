/*
 * cell.c
 *
 * Unit Cell Calculations
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "cell.h"
#include "utils.h"


/* Update the cartesian representation from the crystallographic one */
static void cell_update_cartesian(UnitCell *cell)
{
	double tmp, V, cosalphastar, cstar;

	if ( !cell ) return;

	/* a in terms of x, y and z
	 * +a (cryst) is defined to lie along +x (cart) */
	cell->ax = cell->a;
	cell->ay = 0.0;
	cell->az = 0.0;

	/* b in terms of x, y and z
	 * b (cryst) is defined to lie in the xy (cart) plane */
	cell->bx = cell->b*cos(cell->gamma);
	cell->by = cell->b*sin(cell->gamma);
	cell->bz = 0.0;

	tmp = cos(cell->alpha)*cos(cell->alpha)
		- cos(cell->beta)*cos(cell->beta)
		- cos(cell->gamma)*cos(cell->gamma)
		+ 2.0*cos(cell->alpha)*cos(cell->beta)*cos(cell->gamma);
	V = cell->a * cell->b * cell->c * sqrt(1.0 - tmp);

	cosalphastar = cos(cell->beta)*cos(cell->gamma) - cos(cell->alpha);
	cosalphastar /= sin(cell->beta)*sin(cell->gamma);

	cstar = (cell->a * cell->b * sin(cell->gamma))/V;

	/* c in terms of x, y and z */
	cell->cx = cell->c*cos(cell->beta);
	cell->cy = -cell->c*sin(cell->beta)*cosalphastar;
	cell->cz = 1.0/cstar;
}


/* Update the crystallographic representation from the cartesian one */
static void cell_update_crystallographic(UnitCell *cell)
{
	if ( !cell ) return;

	cell->a = modulus(cell->ax, cell->ay, cell->az);
	cell->b = modulus(cell->bx, cell->by, cell->bz);
	cell->c = modulus(cell->cx, cell->cy, cell->cz);

	cell->alpha = angle_between(cell->bx, cell->by, cell->bz,
					cell->cx, cell->cy, cell->cz);
	cell->beta = angle_between(cell->ax, cell->ay, cell->az,
					cell->cx, cell->cy, cell->cz);
	cell->gamma = angle_between(cell->ax, cell->ay, cell->az,
					cell->bx, cell->by, cell->bz);

	printf("a=%f nm\n", cell->a/1e9);
	printf("b=%f nm\n", cell->b/1e9);
	printf("c=%f nm\n", cell->c/1e9);
	printf("alpha = %f deg\n", rad2deg(cell->alpha));
	printf(" beta = %f deg\n", rad2deg(cell->beta));
	printf("gamma = %f deg\n", rad2deg(cell->gamma));
}


UnitCell *cell_new()
{
	UnitCell *cell;

	cell = malloc(sizeof(UnitCell));
	if ( !cell ) return NULL;

	cell->a = 1.0;
	cell->b = 1.0;
	cell->c = 1.0;
	cell->alpha = M_PI_2;
	cell->beta = M_PI_2;
	cell->gamma = M_PI_2;

	cell_update_cartesian(cell);

	return cell;
}


void cell_set_parameters(UnitCell *cell, double a, double b, double c,
                         double alpha, double beta, double gamma)
{
	if ( !cell ) return;

	cell->a = a;
	cell->b = b;
	cell->c = c;
	cell->alpha = alpha;
	cell->beta = beta;
	cell->gamma = gamma;

	cell_update_cartesian(cell);
}


void cell_set_cartesian(UnitCell *cell,
			double ax, double ay, double az,
			double bx, double by, double bz,
			double cx, double cy, double cz)
{
	if ( !cell ) return;

	cell->ax = ax;  cell->ay = ay;  cell->az = az;
	cell->bx = bx;  cell->by = by;  cell->bz = bz;
	cell->cx = cx;  cell->cy = cy;  cell->cz = cz;

	cell_update_crystallographic(cell);
}


UnitCell *cell_new_from_parameters(double a, double b, double c,
                                   double alpha, double beta, double gamma)
{
	UnitCell *cell;

	cell = cell_new();
	if ( !cell ) return NULL;

	cell_set_parameters(cell, a, b, c, alpha, beta, gamma);

	return cell;
}


void cell_get_cartesian(UnitCell *cell,
                        double *ax, double *ay, double *az,
                        double *bx, double *by, double *bz,
                        double *cx, double *cy, double *cz)
{
	if ( !cell ) return;

	*ax = cell->ax;  *ay = cell->ay;  *az = cell->az;
	*bx = cell->bx;  *by = cell->by;  *bz = cell->bz;
	*cx = cell->cx;  *cy = cell->cy;  *cz = cell->cz;
}


void cell_get_reciprocal(UnitCell *cell,
                         double *asx, double *asy, double *asz,
                         double *bsx, double *bsy, double *bsz,
                         double *csx, double *csy, double *csz)
{
	int s;
	gsl_matrix *m;
	gsl_matrix *inv;
	gsl_permutation *perm;

	m = gsl_matrix_alloc(3, 3);
	gsl_matrix_set(m, 0, 0, cell->ax);
	gsl_matrix_set(m, 0, 1, cell->bx);
	gsl_matrix_set(m, 0, 2, cell->cx);
	gsl_matrix_set(m, 1, 0, cell->ay);
	gsl_matrix_set(m, 1, 1, cell->by);
	gsl_matrix_set(m, 1, 2, cell->cy);
	gsl_matrix_set(m, 2, 0, cell->az);
	gsl_matrix_set(m, 2, 1, cell->bz);
	gsl_matrix_set(m, 2, 2, cell->cz);

	/* Invert */
	perm = gsl_permutation_alloc(m->size1);
	inv = gsl_matrix_alloc(m->size1, m->size2);
	gsl_linalg_LU_decomp(m, perm, &s);
	gsl_linalg_LU_invert(m, perm, inv);
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

	printf("a* = %+8.5f %+8.5f %+8.5f nm^-1\n",
	       *asx/1e9, *asy/1e9, *asz/1e9);
	printf("b* = %+8.5f %+8.5f %+8.5f nm^-1\n",
	       *bsx/1e9, *bsy/1e9, *bsz/1e9);
	printf("c* = %+8.5f %+8.5f %+8.5f nm^-1\n",
	       *csx/1e9, *csy/1e9, *csz/1e9);
}
