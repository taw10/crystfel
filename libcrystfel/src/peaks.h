/*
 * peaks.h
 *
 * Peak search and other image analysis
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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

#ifndef PEAKS_H
#define PEAKS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <pthread.h>

#include "reflist.h"
#include "crystal.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int *make_BgMask(struct image *image, struct panel *p, double ir_inn);

extern void search_peaks(struct image *image, float threshold,
                         float min_gradient, float min_snr,
                         double ir_inn, double ir_mid, double ir_out,
                         int use_saturated);

extern int peak_sanity_check(struct image *image, Crystal **crystals,
                             int n_cryst);

extern void validate_peaks(struct image *image, double min_snr,
                           int ir_inn, int ir_mid, int ir_out,
                           int use_saturated, int check_snr);

#ifdef __cplusplus
}
#endif

#endif	/* PEAKS_H */
