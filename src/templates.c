/*
 * templates.c
 *
 * Indexing by template matching
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "index.h"
#include "index-priv.h"


/* Private data for template indexing */
struct _indexingprivate_template
{
	struct _indexingprivate base;
};


IndexingPrivate *generate_templates(UnitCell *cell, const char *filename)
{
	struct _indexingprivate_template *priv;

	priv = calloc(1, sizeof(struct _indexingprivate_template));
	priv->base.indm = INDEXING_TEMPLATE;

	return (struct _indexingprivate *)priv;
}


void match_templates(struct image *image, IndexingPrivate *ipriv)
{
}
