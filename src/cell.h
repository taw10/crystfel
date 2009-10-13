/*
 * cell.h
 *
 * Unit Cell Calculations
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * template_index - Indexing diffraction patterns by template matching
 *
 */

#ifndef CELL_H
#define CELL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {

	/* Crystallographic representation */
	double a;	/* nm */
	double b;	/* nm */
	double c;	/* nm */
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

} UnitCell;

extern UnitCell *cell_new(void);

/* Lengths in nm, angles in radians */
extern UnitCell *cell_new_from_parameters(double a, double b, double c,
				double alpha, double beta, double gamma);

extern void cell_set_cartesian(UnitCell *cell,
			double ax, double ay, double az,
			double bx, double by, double bz,
			double cx, double cy, double cz);

extern void cell_set_parameters(UnitCell *cell, double a, double b, double c,
				double alpha, double beta, double gamma);

#endif	/* CELL_H */
