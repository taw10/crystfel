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

extern IndexingPrivate *reax_prepare(void);

extern void reax_index(IndexingPrivate *p, struct image *image, UnitCell *cell);


#endif	/* REAX_H */
