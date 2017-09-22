/*
 * asdf.h
 *
 * Alexandra's Superior Direction Finder, or
 * Algorithm Similar to DirAx, FFT-based
 *
 * Copyright © 2014-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2014-2015 Alexandra Tolstikova <alexandra.tolstikova@desy.de>
 *   2015,2017 Thomas White <taw@physics.org>
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

#ifndef ASDF_H
#define ASDF_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "index.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_FFTW

extern int run_asdf(struct image *image, void *ipriv);

extern void *asdf_prepare(IndexingMethod *indm, UnitCell *cell);

extern void asdf_cleanup(void *pp);

#else /* HAVE_FFTW */

int run_asdf(struct image *image, void *ipriv)
{
	ERROR("This copy of CrystFEL was compiled without FFTW support.\n");
	return 0;
}


void *asdf_prepare(IndexingMethod *indm, UnitCell *cell)
{
	ERROR("This copy of CrystFEL was compiled without FFTW support.\n");
	ERROR("To use asdf indexing, recompile with FFTW.\n");
	return NULL;
}

void asdf_cleanup(void *pp)
{
}


#endif /* HAVE_FFTW */


#ifdef __cplusplus
}
#endif

#endif	/* ASDF_H */
