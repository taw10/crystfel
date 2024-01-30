/*
 * merge.h
 *
 * Parallel weighted merging of intensities
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2019 Thomas White <taw@physics.org>
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

#ifndef MERGE_H
#define MERGE_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "crystal.h"
#include "reflist.h"
#include "geometry.h"

/* Minimum partiality of a reflection for it to be merged */
#define MIN_PART_MERGE (0.3)


extern RefList *merge_intensities(struct crystal_refls *crystals, int n, int n_threads,
                                  int min_meas, double push_res, int use_weak,
                                  int ln_merge);

extern double correct_reflection_nopart(double val, Reflection *refl,
                                        double osf, double Bfac, double res);

extern double correct_reflection(double val, Reflection *refl, double osf,
                                 double Bfac, double res);

extern double residual(RefList *list, Crystal *cr, const RefList *full, int free,
                       int *pn_used, const char *filename);

extern double log_residual(RefList *list, Crystal *cr, const RefList *full, int free,
                           int *pn_used, const char *filename);

extern void write_unmerged(const char *fn, struct crystal_refls *crystals, int n_crystals);

#endif	/* MERGE */
