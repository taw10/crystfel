/*
 * statistics.h
 *
 * Structure-factor statistics
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
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

#ifndef STATISTICS_H
#define STATISTICS_H


#include "reflist.h"

extern double stat_scale_intensity(RefList *list1, RefList *list2);

extern double stat_r1_zero(RefList *list1, RefList *list2,
                           double *scalep, int u);
extern double stat_r1_ignore(RefList *list1, RefList *list2,
                             double *scalep, int u);

extern double stat_r2(RefList *list1, RefList *list2, double *scalep, int u);

extern double stat_r1_i(RefList *list1, RefList *list2, double *scalep, int u);

extern double stat_rdiff_zero(RefList *list1, RefList *list2,
                              double *scalep, int u);
extern double stat_rdiff_ignore(RefList *list1, RefList *list2,
                                double *scalep, int u);
extern double stat_rdiff_intensity(RefList *list1, RefList *list2,
                                   double *scalep, int u);

extern double stat_pearson_i(RefList *list1, RefList *list2);
extern double stat_pearson_f_zero(RefList *list1, RefList *list2);
extern double stat_pearson_f_ignore(RefList *list1, RefList *list2);


#endif	/* STATISTICS_H */
