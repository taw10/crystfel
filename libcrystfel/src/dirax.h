/*
 * dirax.h
 *
 * Invoke the DirAx auto-indexing program
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef DIRAX_H
#define DIRAX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "utils.h"


extern void run_dirax(struct image *image);


#endif	/* DIRAX_H */
