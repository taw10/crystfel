/*
 * histogram.h
 *
 * Quick histogram functions
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2013-2014 Thomas White <taw@physics.org>
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

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _histogram Histogram;

extern Histogram *histogram_init();
extern void histogram_free(Histogram *hi);
extern void histogram_add_value(Histogram *hi, double val);
extern void histogram_show(Histogram *hi);

extern int *histogram_get_data(Histogram *hi, int *n);
extern double histogram_get_min(Histogram *hi);
extern double histogram_get_max(Histogram *hi);
extern int histogram_get_num_bins(Histogram *hi);
extern void histogram_set_min(Histogram *hi, double min);
extern void histogram_set_max(Histogram *hi, double max);
extern void histogram_set_num_bins(Histogram *hi, int n);

#ifdef __cplusplus
}
#endif

#endif	/* HISTOGRAM_H */
