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


extern void pr_refine(struct image *image, const RefList *full,
                      const char *sym);


#endif	/* POST_REFINEMENT_H */
