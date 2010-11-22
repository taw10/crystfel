/*
 * post-refinement.h
 *
 * Post refinement
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef POST_REFINEMENT_H
#define POST_REFINEMENT_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>

#include "image.h"


/* Refineable parameters */
enum {
	REF_SCALE,
	REF_DIV,
	REF_R,
	REF_ASX,
	REF_BSX,
	REF_CSX,
	REF_ASY,
	REF_BSY,
	REF_CSY,
	REF_ASZ,
	REF_BSZ,
	REF_CSZ,
	NUM_PARAMS
};

/* Apply the given shift to the 'k'th parameter of 'image'. */
void apply_shift(struct image *image, int k, double shift);


double mean_partial_dev(struct image *image, struct cpeak *spots, int n,
                        const char *sym, double *i_full, FILE *graph);


double pr_iterate(struct image *image, double *i_full, const char *sym,
                  struct cpeak **pspots, int *n);


#endif	/* POST_REFINEMENT_H */
