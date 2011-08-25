/*
 * index.c
 *
 * Perform indexing (somehow)
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 * (c) 2010 Richard Kirian <rkirian@asu.edu>
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
#include <string.h>
#include <assert.h>

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "dirax.h"
#include "mosflm.h"
#include "detector.h"
#include "index.h"
#include "index-priv.h"
#include "reax.h"


/* Base class constructor for unspecialised indexing private data */
IndexingPrivate *indexing_private(IndexingMethod indm)
{
	struct _indexingprivate *priv;
	priv = calloc(1, sizeof(struct _indexingprivate));
	priv->indm = indm;
	return priv;
}


static const char *maybes(int n)
{
	if ( n == 1 ) return "";
	return "s";
}


IndexingPrivate **prepare_indexing(IndexingMethod *indm, UnitCell *cell,
                                   const char *filename, struct detector *det,
                                   double nominal_photon_energy)
{
	int n;
	int nm = 0;
	IndexingPrivate **iprivs;

	while ( indm[nm] != INDEXING_NONE ) nm++;
	STATUS("Preparing %i indexing method%s.\n", nm, maybes(nm));
	iprivs = malloc((nm+1) * sizeof(IndexingPrivate *));

	for ( n=0; n<nm; n++ ) {

		switch ( indm[n] ) {
		case INDEXING_NONE :
			ERROR("Tried to prepare INDEXING_NONE!\n");
			break;
		case INDEXING_DIRAX :
			iprivs[n] = indexing_private(indm[n]);
			break;
		case INDEXING_MOSFLM :
			iprivs[n] = indexing_private(indm[n]);
			break;
		case INDEXING_REAX :
			iprivs[n] = reax_prepare();
			break;
		}

	}
	iprivs[n] = NULL;

	return iprivs;
}


void cleanup_indexing(IndexingPrivate **priv)
{
	int n = 0;

	if ( priv == NULL ) return;  /* Nothing to do */

	while ( priv[n] != NULL ) {

		switch ( priv[n]->indm ) {
		case INDEXING_NONE :
			free(priv[n]);
			break;
		case INDEXING_DIRAX :
			free(priv[n]);
			break;
		case INDEXING_MOSFLM :
			free(priv[n]);
			break;
		case INDEXING_REAX :
			reax_cleanup(priv[n]);
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


void index_pattern(struct image *image, UnitCell *cell, IndexingMethod *indm,
                   int cellr, int verbose, IndexingPrivate **ipriv,
                   int config_insane)
{
	int i;
	int n = 0;

	if ( indm == NULL ) return;

	map_all_peaks(image);
	image->indexed_cell = NULL;

	while ( indm[n] != INDEXING_NONE ) {

		image->ncells = 0;

		/* Index as appropriate */
		switch ( indm[n] ) {
		case INDEXING_NONE :
			return;
		case INDEXING_DIRAX :
			run_dirax(image);
			break;
		case INDEXING_MOSFLM :
			run_mosflm(image, cell);
			break;
		case INDEXING_REAX :
			reax_index(ipriv[n], image, cell);
			break;
		}
		if ( image->ncells == 0 ) {
			n++;
			continue;
		}

		for ( i=0; i<image->ncells; i++ ) {

			UnitCell *new_cell = NULL;
			UnitCell *cand = image->candidate_cells[i];

			if ( verbose ) {
				STATUS("--------------------\n");
				STATUS("Candidate cell %i (before matching):\n",
				       i);
				cell_print(image->candidate_cells[i]);
				STATUS("--------------------\n");
			}

			/* Match or reduce the cell as appropriate */
			switch ( cellr ) {
			case CELLR_NONE :
				new_cell = cell_new_from_cell(cand);
				break;
			case CELLR_REDUCE :
				new_cell = match_cell(cand, cell, verbose, 1);
				break;
			case CELLR_COMPARE :
				new_cell = match_cell(cand, cell, verbose, 0);
				break;
			case CELLR_COMPARE_AB :
				new_cell = match_cell_ab(cand, cell);
				break;
			}

			/* No cell?  Move on to the next candidate */
			if ( new_cell == NULL ) continue;

			/* Sanity check */
			if ( !config_insane &&
			     !peak_sanity_check(image, new_cell, 0, 0.1) ) {
				cell_free(new_cell);
				continue;
			}

			image->indexed_cell = new_cell;
			goto done;  /* Success */

		}

		for ( i=0; i<image->ncells; i++ ) {
			cell_free(image->candidate_cells[i]);
			image->candidate_cells[i] = NULL;
		}

		/* Move on to the next indexing method */
		n++;

	}

done:
	for ( i=0; i<image->ncells; i++ ) {
		/* May free(NULL) if all algorithms were tried and no success */
		cell_free(image->candidate_cells[i]);
	}
}


IndexingMethod *build_indexer_list(const char *str, int *need_cell)
{
	int n, i;
	char **methods;
	IndexingMethod *list;
	int tmp;

	if ( need_cell == NULL ) need_cell = &tmp;
	*need_cell = 0;

	n = assplode(str, ",", &methods, ASSPLODE_NONE);
	list = malloc((n+1)*sizeof(IndexingMethod));

	for ( i=0; i<n; i++ ) {

		if ( strcmp(methods[i], "dirax") == 0) {
			list[i] = INDEXING_DIRAX;
		} else if ( strcmp(methods[i], "mosflm") == 0) {
			list[i] = INDEXING_MOSFLM;
		} else if ( strcmp(methods[i], "reax") == 0) {
			list[i] = INDEXING_REAX;
			*need_cell = 1;
		} else {
			ERROR("Unrecognised indexing method '%s'\n",
			      methods[i]);
			return NULL;
		}

		free(methods[i]);

	}
	free(methods);
	list[i] = INDEXING_NONE;

	return list;
}
