/*
 * post-refinement.h
 *
 * Post refinement
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
#include "utils.h"


/* Refineable parameters */
enum {
	REF_ASX,
	REF_ASY,
	REF_ASZ,
	REF_BSX,
	REF_BSY,
	REF_BSZ,
	REF_CSX,
	REF_CSY,
	REF_CSZ,
	NUM_PARAMS,
	REF_DIV,
	REF_R,
};


extern void pr_refine(struct image *image, const RefList *full,
                      const char *sym);

/* Exported so it can be poked by tests/pr_gradient_check */
extern double gradient(struct image *image, int k, Reflection *refl, double r);


#endif	/* POST_REFINEMENT_H */
