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


typedef enum {
	INDEXING_NONE,
	INDEXING_DIRAX,
	INDEXING_MATCH
} IndexingMethod;


extern void index_pattern(struct image *image, IndexingMethod indm,
                          int no_match, int verbose);

#endif	/* INDEX_H */
