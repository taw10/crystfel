/*
 * reax.h
 *
 * A new auto-indexer
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
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

#ifndef REAX_H
#define REAX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cell.h"

#ifdef HAVE_FFTW

extern IndexingPrivate *reax_prepare(void);
extern void reax_cleanup(IndexingPrivate *pp);
extern void reax_index(IndexingPrivate *p, struct image *image, UnitCell *cell);

#else /* HAVE_FFTW */

static IndexingPrivate *reax_prepare()
{
	return NULL;
}

static void reax_cleanup(IndexingPrivate *pp)
{
}

static void reax_index(IndexingPrivate *p, struct image *image, UnitCell *cell)
{
}


#endif /* HAVE_FFTW */

#endif	/* REAX_H */
