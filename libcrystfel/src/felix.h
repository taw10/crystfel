/*
 * felix.h
 *
 * Invoke Felix for multi-crystal autoindexing
 *
 * Copyright © 2013-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013,2017 Thomas White <taw@physics.org>
 *   2013-2014      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
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

#ifndef FELIX_H
#define FELIX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cell.h"

struct felix_options
{
	double ttmin;  /* radians */
	double ttmax;  /* radians */
	int min_measurements;
	double min_completeness;
	double min_uniqueness;
	int n_voxels;
	double test_fraction;
	double sigma;
};

extern void *felix_prepare(IndexingMethod *indm, UnitCell *cell,
                           struct felix_options *opts);

extern const char *felix_probe(UnitCell *cell);

extern void felix_cleanup(IndexingPrivate *pp);

extern int felix_index(struct image *image, IndexingPrivate *p);


#endif	/* FELIX_H */
