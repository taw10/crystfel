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

struct fom_rejections
{
	int common;
	int low_snr;
	int negative_deleted;
	int negative_zeroed;
	int few_measurements;
	int outside_resolution_range;
	int no_bijvoet;
	int centric;
	int nan_inf_value;
};

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
	FOM_D2SIG,
	FOM_NUM_MEASUREMENTS,
	FOM_REDUNDANCY,
	FOM_SNR,
	FOM_MEAN_INTENSITY,
	FOM_STDDEV_INTENSITY,
	FOM_COMPLETENESS,
};

struct fom_shells
{
	int nshells;
	double *rmins;
	double *rmaxs;
};

struct fom_context;

extern struct fom_rejections fom_select_reflection_pairs(RefList *list1,
                                                         RefList *list2,
                                                         RefList **plist1_acc,
                                                         RefList **plist2_acc,
                                                         UnitCell *cell,
                                                         SymOpList *sym,
                                                         int anom,
                                                         double rmin_fix,
                                                         double rmax_fix,
                                                         double sigma_cutoff,
                                                         int ignore_negs,
                                                         int zero_negs,
                                                         int mul_cutoff);

extern struct fom_rejections fom_select_reflections(RefList *list,
                                                    RefList **plist_acc,
                                                    UnitCell *cell,
                                                    SymOpList *sym,
                                                    double rmin_fix,
                                                    double rmax_fix,
                                                    double sigma_cutoff,
                                                    int ignore_negs,
                                                    int zero_negs,
                                                    int mul_cutoff);


extern struct fom_context *fom_calculate(RefList *list1, RefList *list2,
                                         UnitCell *cell,
                                         struct fom_shells *shells,
                                         enum fom_type fom, int noscale,
                                         const SymOpList *sym);

extern struct fom_shells *fom_make_resolution_shells(double rmin, double rmax,
                                                     int nshells);

extern double fom_shell_centre(struct fom_shells *s, int i);

extern double fom_overall_value(struct fom_context *fctx);
extern double fom_shell_value(struct fom_context *fctx, int i);

extern int fom_overall_num_reflections(struct fom_context *fctx);
extern int fom_shell_num_reflections(struct fom_context *fctx, int i);

extern int fom_overall_num_possible(struct fom_context *fctx);
extern int fom_shell_num_possible(struct fom_context *fctx, int i);

extern enum fom_type fom_type_from_string(const char *s);

#endif	/* FOM */
