/*
 * reax.c
 *
 * A new auto-indexer
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell.h"
#include "index.h"
#include "index-priv.h"


struct dvec
{
	double x;
	double y;
	double z;
};


struct reax_private
{
	IndexingPrivate base;
};


IndexingPrivate *reax_prepare()
{
	struct reax_private *priv;

	priv = calloc(1, sizeof(*priv));
	if ( priv == NULL ) return NULL;

	priv->base.indm = INDEXING_REAX;

	return (IndexingPrivate *)priv;
}


void reax_index(struct image *image, UnitCell *cell)
{

}
