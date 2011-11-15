/*
 * mosflm.h
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * (c) 2010 Richard Kirian <rkirian@asu.edu>
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef MOSFLM_H
#define MOSFLM_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "utils.h"


extern void run_mosflm(struct image *image, UnitCell *cell);


#endif	/* MOSFLM_H */
