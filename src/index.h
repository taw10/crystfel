/*
 * index.h
 *
 * Perform indexing (somehow)
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef INDEX_H
#define INDEX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


extern void index_pattern(struct image *image, int no_index, int dump_peaks,
                          int use_dirax);


#endif	/* INDEX_H */
