/*
 * fom.h
 *
 * Figure of merit calculation
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
 *   2013      Lorenzo Galli <lorenzo.galli@desy.de>
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

#ifndef FOM_H
#define FOM_H

/**
 * \file fom.h
 * Figure of merit calculation
 */

#include <reflist.h>
#include <symmetry.h>

enum fom_type
{
	FOM_R1I,
	FOM_R1F,
	FOM_R2,
	FOM_RSPLIT,
	FOM_CC,
	FOM_CCSTAR,
	FOM_CCANO,
	FOM_CRDANO,
	FOM_RANO,
	FOM_RANORSPLIT,
	FOM_D1SIG,
	FOM_D2SIG
};

struct fom_shells
{
	int config_intshells;
	int nshells;
	double *rmins;
	double *rmaxs;
};

struct fom_context
{
	enum fom_type fom;
	int nshells;
	int *cts;

	/* For R-factors */
	double *num;
	double *den;

	/* For "double" R-factors */
	double *num2;
	double *den2;

	/* For CCs */
	double **vec1;
	double **vec2;
	int *n;
	int nmax;

	/* For "counting" things e.g. d1sig or d2sig */
	int *n_within;
};


extern struct fom_context *fom_calculate(RefList *list1, RefList *list2,
                                         UnitCell *cell,
                                         struct fom_shells *shells,
                                         enum fom_type fom, int noscale,
                                         SymOpList *sym);

extern struct fom_shells *fom_make_resolution_shells(double rmin, double rmax,
                                                     int nshells);

extern struct fom_shells *fom_make_intensity_shells(double min_I, double max_I,
                                                    int nshells);

extern double fom_shell_label(struct fom_shells *s, int i);

extern double fom_shell(struct fom_context *fctx, int i);

extern double fom_overall(struct fom_context *fctx);

extern enum fom_type fom_type_from_string(const char *s);

#endif	/* FOM */
