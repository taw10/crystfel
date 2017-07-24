/*
 * taketwo.h
 *
 * Rewrite of TakeTwo algorithm (Acta D72 (8) 956-965) for CrystFEL
 *
 * Copyright © 2016-2017 Helen Ginn
 * Copyright © 2016-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2016      Helen Ginn <helen@strubi.ox.ac.uk>
 *   2016-2017 Thomas White <taw@physics.org>
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

#ifndef TAKETWO_H
#define TAKETWO_H

#include "cell.h"
#include "index.h"

struct taketwo_options
{
	int member_thresh;
	double len_tol;
	double angle_tol;
	double trace_tol;
};


extern void *taketwo_prepare(IndexingMethod *indm, UnitCell *cell,
                             struct detector *det, float *ltl);
extern int taketwo_index(struct image *image,
                         const struct taketwo_options *opts, void *priv);
extern void taketwo_cleanup(IndexingPrivate *pp);

#endif /* TAKETWO_H */
