/*
 * templates.h
 *
 * Indexing by template matching
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef TEMPLATES_H
#define TEMPLATES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "index.h"
#include "image.h"
#include "cell.h"

extern IndexingPrivate *generate_templates(UnitCell *cell, const char *filename,
                                           struct detector *det);


extern void match_templates(struct image *image, IndexingPrivate *ipriv);

extern void free_templates(IndexingPrivate *priv);

#endif	/* TEMPLATES_H */
