/*
 * dirax.h
 *
 * Invoke the DirAx auto-indexing program
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef DIRAX_H
#define DIRAX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


extern void index_pattern(struct image *image, int no_index);


#endif	/* DIRAX_H */
