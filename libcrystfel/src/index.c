/*
 * index.c
 *
 * Perform indexing (somehow)
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
 *   2010-2011 Richard Kirian <rkirian@asu.edu>
 *   2012      Lorenzo Galli
 *   2013      Cornelius Gati <cornelius.gati@cfel.de>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "dirax.h"
#include "mosflm.h"
#include "xds.h"
#include "detector.h"
#include "index.h"
#include "reax.h"
#include "grainspotter.h"
#include "geometry.h"
#include "cell-utils.h"
#include "grainspotter.h"


static int debug_index(struct image *image)
{
	Crystal *cr = crystal_new();
	UnitCell *cell = cell_new();
	cell_set_reciprocal(cell, +0.0000e9, +0.0000e9, +0.0000e9,
	                          +0.0000e9, +0.0000e9, +0.0000e9,
	                          +0.0000e9, +0.0000e9, +0.0000e9);
	crystal_set_cell(cr, cell);
	image_add_crystal(image, cr);
	return 1;
}


IndexingPrivate **prepare_indexing(IndexingMethod *indm, UnitCell *cell,
                                   struct detector *det, float *ltl)
{
	int n;
	int nm = 0;
	IndexingPrivate **iprivs;

	while ( indm[nm] != INDEXING_NONE ) nm++;
	iprivs = malloc((nm+1) * sizeof(IndexingPrivate *));

	for ( n=0; n<nm; n++ ) {

		int i;
		IndexingMethod in;
		char *str;

		in = indm[n];

		switch ( indm[n] & INDEXING_METHOD_MASK ) {

			case INDEXING_DIRAX :
			iprivs[n] = dirax_prepare(&indm[n], cell, det, ltl);
			break;

			case INDEXING_MOSFLM :
			iprivs[n] = mosflm_prepare(&indm[n], cell, det, ltl);
			break;

			case INDEXING_XDS :
                        iprivs[n] = xds_prepare(&indm[n], cell, det, ltl);
			break;

			case INDEXING_REAX :
			iprivs[n] = reax_prepare(&indm[n], cell, det, ltl);
			break;

			case INDEXING_GRAINSPOTTER :
			iprivs[n] = grainspotter_prepare(&indm[n], cell,
			                                 det, ltl);
			break;

			case INDEXING_DEBUG :
			iprivs[n] = (IndexingPrivate *)strdup("Hello!");
			break;

			default :
			ERROR("Don't know how to prepare indexing method %i\n",
			      indm[n]);
			break;

		}

		if ( iprivs[n] == NULL ) return NULL;

		str = indexer_str(indm[n]);
		STATUS("Prepared indexing method %i: %s\n", n, str);
		free(str);

		if ( in != indm[n] ) {
			ERROR("Note: flags were altered to take into account "
			      "the limitations of the indexing method.\n");
		}

		for ( i=0; i<n; i++ ) {
			if ( indm[i] == indm[n] ) {
				ERROR("Duplicate indexing method.\n");
				ERROR("Have you specified some flags which "
				      "aren't accepted by one of your "
				      "chosen indexing methods?\n");
				return NULL;
			}
		}


	}
	iprivs[n] = NULL;

	return iprivs;
}


void cleanup_indexing(IndexingMethod *indms, IndexingPrivate **privs)
{
	int n = 0;

	if ( indms == NULL ) return;  /* Nothing to do */
	if ( privs == NULL ) return;  /* Nothing to do */

	while ( indms[n] != INDEXING_NONE ) {

		switch ( indms[n] & INDEXING_METHOD_MASK ) {

			case INDEXING_NONE :
			break;

			case INDEXING_DIRAX :
			dirax_cleanup(privs[n]);
			break;

			case INDEXING_MOSFLM :
			mosflm_cleanup(privs[n]);
			break;

                        case INDEXING_XDS :
			xds_cleanup(privs[n]);
			break;

			case INDEXING_REAX :
			reax_cleanup(privs[n]);
			break;

			case INDEXING_GRAINSPOTTER :
			grainspotter_cleanup(privs[n]);
			break;

			case INDEXING_DEBUG :
			free(privs[n]);
			break;

			default :
			ERROR("Don't know how to clean up indexing method %i\n",
			      indms[n]);
			break;

		}

		n++;

	}
}


void map_all_peaks(struct image *image)
{
	int i, n;

	n = image_feature_count(image->features);

	/* Map positions to 3D */
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		struct rvec r;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);
		f->rx = r.u;  f->ry = r.v;  f->rz = r.w;

	}
}


/* Return non-zero for "success" */
static int try_indexer(struct image *image, IndexingMethod indm,
                       IndexingPrivate *ipriv)
{
	switch ( indm & INDEXING_METHOD_MASK ) {

		case INDEXING_NONE :
		return 0;
		break;

		case INDEXING_DIRAX :
		return run_dirax(image, ipriv);
		break;

		case INDEXING_MOSFLM :
		return run_mosflm(image, ipriv);
		break;

		case INDEXING_XDS :
		return run_xds(image, ipriv);
		break;

		case INDEXING_REAX :
		return reax_index(ipriv, image);
		break;

		case INDEXING_GRAINSPOTTER :
		return grainspotter_index(image, ipriv);
		break;

		case INDEXING_DEBUG :
		return debug_index(image);

		default :
		ERROR("Unrecognised indexing method: %i\n", indm);
		break;

	}

	return 0;
}


void index_pattern(struct image *image,
                   IndexingMethod *indms, IndexingPrivate **iprivs)
{
	int n = 0;

	if ( indms == NULL ) return;

	if ( image_feature_count(image->features) > 10000 ) {
		STATUS("WARNING: The number of peaks is very large for '%s'.\n",
		       image->filename);
	}

	map_all_peaks(image);
	image->crystals = NULL;
	image->n_crystals = 0;

	while ( indms[n] != INDEXING_NONE ) {

		if ( try_indexer(image, indms[n], iprivs[n]) ) break;
		n++;

	}

	image->indexed_by = indms[n];
}


/* Set the indexer flags for "raw mode" ("--cell-reduction=none") */
static IndexingMethod set_raw(IndexingMethod a)
{
	/* Disable all unit cell checks */
	a &= ~(INDEXING_CHECK_CELL_COMBINATIONS | INDEXING_CHECK_CELL_AXES);
	return a;
}


/* Set the indexer flags for "bad mode" ("--insane) */
static IndexingMethod set_bad(IndexingMethod a)
{
	/* Disable the peak check */
	return a & ~INDEXING_CHECK_PEAKS;
}


/* Set the indexer flags for "axes mode" ("--cell-reduction=compare") */
static IndexingMethod set_axes(IndexingMethod a)
{
	return (a & ~INDEXING_CHECK_CELL_COMBINATIONS)
	          | INDEXING_CHECK_CELL_AXES;
}


/* Set the indexer flags for "combination mode" ("--cell-reduction=reduce") */
static IndexingMethod set_comb(IndexingMethod a)
{
	return (a & ~INDEXING_CHECK_CELL_AXES)
	          | INDEXING_CHECK_CELL_COMBINATIONS;
}


/* Set the indexer flags for "use no lattice type information" */
static IndexingMethod set_nolattice(IndexingMethod a)
{
	return a & ~INDEXING_USE_LATTICE_TYPE;
}


/* Set the indexer flags for "use lattice type information" */
static IndexingMethod set_lattice(IndexingMethod a)
{
	return a | INDEXING_USE_LATTICE_TYPE;
}


/* Set the indexer flags for "use no unit cell parameters" */
static IndexingMethod set_nocellparams(IndexingMethod a)
{
	return a & ~INDEXING_USE_CELL_PARAMETERS;
}


/* Set the indexer flags for "use unit cell parameters" */
static IndexingMethod set_cellparams(IndexingMethod a)
{
	return a | INDEXING_USE_CELL_PARAMETERS;
}

char *indexer_str(IndexingMethod indm)
{
	char *str;

	str = malloc(32);
	if ( str == NULL ) {
		ERROR("Failed to allocate string.\n");
		return NULL;
	}
	str[0] = '\0';

	switch ( indm & INDEXING_METHOD_MASK ) {

		case INDEXING_NONE :
		strcpy(str, "none");
		return str;

		case INDEXING_DIRAX :
		strcpy(str, "dirax");
		break;

		case INDEXING_MOSFLM :
		strcpy(str, "mosflm");
		break;

		case INDEXING_REAX :
		strcpy(str, "reax");
		break;

		case INDEXING_GRAINSPOTTER :
		strcpy(str, "grainspotter");
		break;

		case INDEXING_XDS :
		strcpy(str, "xds");
		break;

		case INDEXING_SIMULATION :
		strcpy(str, "simulation");
		break;

		case INDEXING_DEBUG :
		strcpy(str, "debug");
		break;

		default :
		ERROR("No test description for indexing method %i\n",
		      indm & INDEXING_METHOD_MASK);
		strcpy(str, "(unknown)");
		break;

	}

	if ( (indm & INDEXING_METHOD_MASK) == INDEXING_SIMULATION ) return str;

	if ( indm & INDEXING_CHECK_CELL_COMBINATIONS ) {
		strcat(str, "-comb");
	} else if ( indm & INDEXING_CHECK_CELL_AXES ) {
		strcat(str, "-axes");
	} else {
		strcat(str, "-raw");
	}

	if ( !(indm & INDEXING_CHECK_PEAKS) ) {
		strcat(str, "-bad");
	}

	if ( indm & INDEXING_USE_LATTICE_TYPE ) {
		strcat(str, "-latt");
	} else {
		strcat(str, "-nolatt");
	}

	if ( indm & INDEXING_USE_CELL_PARAMETERS ) {
		strcat(str, "-cell");
	} else {
		strcat(str, "-nocell");
	}

	return str;
}


IndexingMethod *build_indexer_list(const char *str)
{
	int n, i;
	char **methods;
	IndexingMethod *list;
	int nmeth = 0;

	n = assplode(str, ",-", &methods, ASSPLODE_NONE);
	list = malloc((n+1)*sizeof(IndexingMethod));

	nmeth = -1;  /* So that the first method is #0 */
	for ( i=0; i<n; i++ ) {

		if ( strcmp(methods[i], "dirax") == 0) {
			list[++nmeth] = INDEXING_DEFAULTS_DIRAX;

		} else if ( strcmp(methods[i], "mosflm") == 0) {
			list[++nmeth] = INDEXING_DEFAULTS_MOSFLM;

		} else if ( strcmp(methods[i], "grainspotter") == 0) {
			list[++nmeth] = INDEXING_DEFAULTS_GRAINSPOTTER;

                } else if ( strcmp(methods[i], "xds") == 0) {
			list[++nmeth] = INDEXING_DEFAULTS_XDS;

		} else if ( strcmp(methods[i], "reax") == 0) {
			list[++nmeth] = INDEXING_DEFAULTS_REAX;

		} else if ( strcmp(methods[i], "none") == 0) {
			list[++nmeth] = INDEXING_NONE;

		} else if ( strcmp(methods[i], "simulation") == 0) {
			list[++nmeth] = INDEXING_SIMULATION;
			return list;

		} else if ( strcmp(methods[i], "debug") == 0) {
			list[++nmeth] = INDEXING_DEBUG;
			return list;

		} else if ( strcmp(methods[i], "raw") == 0) {
			list[nmeth] = set_raw(list[nmeth]);

		} else if ( strcmp(methods[i], "bad") == 0) {
			list[nmeth] = set_bad(list[nmeth]);

		} else if ( strcmp(methods[i], "comb") == 0) {
			list[nmeth] = set_comb(list[nmeth]);  /* Default */

		} else if ( strcmp(methods[i], "axes") == 0) {
			list[nmeth] = set_axes(list[nmeth]);

		} else if ( strcmp(methods[i], "latt") == 0) {
			list[nmeth] = set_lattice(list[nmeth]);

		} else if ( strcmp(methods[i], "nolatt") == 0) {
			list[nmeth] = set_nolattice(list[nmeth]);

		} else if ( strcmp(methods[i], "cell") == 0) {
			list[nmeth] = set_cellparams(list[nmeth]);

		} else if ( strcmp(methods[i], "nocell") == 0) {
			list[nmeth] = set_nocellparams(list[nmeth]);

		} else {
			ERROR("Bad list of indexing methods: '%s'\n", str);
			return NULL;
		}

		free(methods[i]);

	}
	free(methods);
	list[++nmeth] = INDEXING_NONE;

	return list;
}
