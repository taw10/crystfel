/*
 * reax.h
 *
 * A new auto-indexer
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
