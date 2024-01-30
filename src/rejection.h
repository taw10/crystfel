/*
 * rejection.h
 *
 * Crystal rejection for scaling
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

#ifndef REJECTION_H
#define REJECTION_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "image.h"

extern void early_rejection(struct crystal_refls *crystals, int n);
extern void check_rejection(struct crystal_refls *crystals, int n, RefList *full,
                            double max_B, int no_deltacchalf, int n_threads);

#endif	/* REJECTION_H */
