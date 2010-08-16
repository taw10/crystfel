/*
 * index-priv.h
 *
 * Indexing private data
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef INDEXPRIV_H
#define INDEXPRIV_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "index.h"

struct _indexingprivate
{
	IndexingMethod indm;
};


#endif	/* INDEXPRIV_H */
