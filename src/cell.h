/*
 * cell.h
 *
 * Unit Cell Calculations
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef CELL_H
#define CELL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* A 3D vector in reciprocal space */
struct rvec
{
	double   u;
	double   v;
	double   w;
};

typedef struct _unitcell UnitCell;

extern UnitCell *cell_new(void);
extern UnitCell *cell_new_from_cell(UnitCell *orig);
extern void cell_free(UnitCell *cell);

/* Lengths in m, angles in radians */
extern UnitCell *cell_new_from_parameters(double a, double b, double c,
				double alpha, double beta, double gamma);

extern UnitCell *cell_new_from_axes(struct rvec as, struct rvec bs,
                                    struct rvec cs);

extern void cell_set_cartesian(UnitCell *cell,
			double ax, double ay, double az,
			double bx, double by, double bz,
			double cx, double cy, double cz);

extern void cell_set_parameters(UnitCell *cell, double a, double b, double c,
				double alpha, double beta, double gamma);

extern void cell_set_cartesian_a(UnitCell *cell, double ax, double ay, double az);
extern void cell_set_cartesian_b(UnitCell *cell, double bx, double by, double bz);
extern void cell_set_cartesian_c(UnitCell *cell, double cx, double cy, double cz);


extern int cell_get_parameters(UnitCell *cell, double *a, double *b, double *c,
                                double *alpha, double *beta, double *gamma);

extern int cell_get_cartesian(UnitCell *cell,
                               double *ax, double *ay, double *az,
                               double *bx, double *by, double *bz,
                               double *cx, double *cy, double *cz);

extern int cell_get_reciprocal(UnitCell *cell,
                               double *asx, double *asy, double *asz,
                               double *bsx, double *bsy, double *bsz,
                               double *csx, double *csy, double *csz);

extern void cell_set_reciprocal(UnitCell *cell,
                               double asx, double asy, double asz,
                               double bsx, double bsy, double bsz,
                               double csx, double csy, double csz);

extern const char *cell_get_pointgroup(UnitCell *cell);

extern double resolution(UnitCell *cell,
                         signed int h, signed int k, signed int l);

extern void cell_print(UnitCell *cell);

extern UnitCell *match_cell(UnitCell *cell, UnitCell *template, int verbose);
extern int cells_similar(UnitCell *c1, UnitCell *c2);

extern UnitCell *load_cell_from_pdb(const char *filename);

#endif	/* CELL_H */
