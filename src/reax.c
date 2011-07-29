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
	struct dvec *directions;
};


IndexingPrivate *reax_prepare()
{
	struct reax_private *priv;
	int ui, vi;
	int samp;
	double angular_inc;

	priv = calloc(1, sizeof(*priv));
	if ( priv == NULL ) return NULL;

	priv->base.indm = INDEXING_REAX;

	/* Decide on sampling interval */
	angular_inc = 0.03;  /* From Steller (1997) */
	samp = (2.0 * M_PI) / angular_inc;

	priv->directions = malloc(samp*samp*sizeof(struct dvec));
	if ( priv == NULL) {
		free(priv);
		return NULL;
	}

	for ( ui=0; ui<samp; ui++ ) {
	for ( vi=0; vi<samp; vi++ ) {

		double u, v;
		double th, ph;
		struct dvec *dir;

		u = (double)ui/samp;
		v = (double)vi/samp;

		th = 2.0 * M_PI * u;
		ph = acos(2.0*v - 1.0);

		dir = &priv->directions[ui + vi*samp];

		dir->x = cos(th) * sin(ph);
		dir->y = sin(th) * sin(th);
		dir->z = cos(ph);

	}
	}

	return (IndexingPrivate *)priv;
}


void reax_index(struct image *image, UnitCell *cell)
{

}
