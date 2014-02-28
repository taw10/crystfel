/*
 * multihistogram.h
 *
 * Histogram with categories
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2013-2014 Thomas White <taw@physics.org>
 *
 * This program is free software: you can redistribute it and/or modify
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

#ifndef MULTIHISTOGRAM_H
#define MULTIHISTOGRAM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


typedef struct _multihistogram MultiHistogram;

extern MultiHistogram *multihistogram_new();
extern void multihistogram_free(MultiHistogram *hi);

extern void multihistogram_delete_all_values(MultiHistogram *hi);
extern void multihistogram_add_value(MultiHistogram *hi, double val, int cat);

extern void multihistogram_set_min(MultiHistogram *hi, double min);
extern void multihistogram_set_max(MultiHistogram *hi, double max);
extern void multihistogram_set_num_bins(MultiHistogram *hi, int n);

extern int *multihistogram_get_data(MultiHistogram *hi, int cat);


#endif	/* HISTOGRAM_H */
