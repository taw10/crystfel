/*
 * statistics.h
 *
 * Structure-factor statistics
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef STATISTICS_H
#define STATISTICS_H


#include "reflist.h"

extern double stat_scale_intensity(RefList *list1, double *arr2);

extern double stat_r1_zero(RefList *list1, double *arr2, double *scalep, int u);
extern double stat_r1_ignore(RefList *list1, double *arr2,
                             double *scalep, int u);

extern double stat_r2(RefList *list1, double *arr2, double *scalep, int u);

extern double stat_r1_i(RefList *list1, double *arr2, double *scalep, int u);

extern double stat_rdiff_zero(RefList *list1, double *arr2,
                              double *scalep, int u);
extern double stat_rdiff_ignore(RefList *list1, double *arr2,
                                double *scalep, int u);
extern double stat_rdiff_intensity(RefList *list1, double *arr2,
                                   double *scalep, int u);

extern double stat_pearson_i(RefList *list1, double *arr2);
extern double stat_pearson_f_zero(RefList *list1, double *arr2);
extern double stat_pearson_f_ignore(RefList *list1, double *arr2);


#endif	/* STATISTICS_H */
