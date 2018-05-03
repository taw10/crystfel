/*
 * scaling.h
 *
 * Scaling
 *
 * Copyright © 2012-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2018 Thomas White <taw@physics.org>
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

#ifndef SCALING_H
#define SCALING_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "crystal.h"
#include "geometry.h"

enum ScaleFlags
{
	SCALE_NONE = 0,
	SCALE_NO_B = 1<<0,  /* Don't use Debye-Waller part */
};

extern int scale_one_crystal(Crystal *cr, const RefList *reference, int flags);

extern void scale_all(Crystal **crystals, int n_crystals, int nthreads,
                      int flags);

#endif	/* SCALING_H */
