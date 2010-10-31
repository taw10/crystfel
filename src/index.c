/*
 * index.c
 *
 * Perform indexing (somehow)
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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
#include "sfac.h"
#include "detector.h"
#include "index.h"
#include "index-priv.h"
#include "templates.h"


/* Base class constructor for unspecialised indexing private data */
static IndexingPrivate *indexing_private(IndexingMethod indm)
{
	struct _indexingprivate *priv;
	priv = calloc(1, sizeof(struct _indexingprivate));
	priv->indm = indm;
	return priv;
}


IndexingPrivate *prepare_indexing(IndexingMethod indm, UnitCell *cell,
                                  const char *filename, struct detector *det,
                                  double nominal_photon_energy)
{
	switch ( indm ) {
	case INDEXING_NONE :
		return indexing_private(indm);
	case INDEXING_DIRAX :
		return indexing_private(indm);
	case INDEXING_TEMPLATE :
		return generate_templates(cell, filename, det,
		                          nominal_photon_energy);
	}
	return 0;
}


void cleanup_indexing(IndexingPrivate *priv)
{
	switch ( priv->indm ) {
	case INDEXING_NONE :
		free(priv);
		break;
	case INDEXING_DIRAX :
		free(priv);
		break;
	case INDEXING_TEMPLATE :
		free_templates(priv);
	}
}


static void write_drx(struct image *image)
{
	FILE *fh;
	int i;
	char filename[1024];

	snprintf(filename, 1023, "xfel-%i.drx", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file xfel.drx\n");
		return;
	}
	fprintf(fh, "%f\n", 0.5);  /* Lie about the wavelength.  */

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		fprintf(fh, "%10f %10f %10f %8f\n",
		        f->rx/1e10, f->ry/1e10, f->rz/1e10, 1.0);

	}
	fclose(fh);
}


void map_all_peaks(struct image *image)
{
	int i;

	/* Map positions to 3D */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec r;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->x, f->y, 1, NULL, 1.0/image->lambda);
		f->rx = r.u;  f->ry = r.v;  f->rz = r.w;

	}
}


void index_pattern(struct image *image, UnitCell *cell, IndexingMethod indm,
                   int cellr, int verbose, IndexingPrivate *ipriv)
{
	int i;

	map_all_peaks(image);
	write_drx(image);

	image->ncells = 0;

	/* Index (or not) as appropriate */
	switch ( indm ) {
	case INDEXING_NONE :
		return;
	case INDEXING_DIRAX :
		run_dirax(image);
		break;
	case INDEXING_TEMPLATE :
		match_templates(image, ipriv);
		break;
	}

	if ( image->ncells == 0 ) {
		STATUS("No candidate cells found.\n");
		return;
	}

	if ( (cellr == CELLR_NONE) || (indm == INDEXING_TEMPLATE) ) {
		image->indexed_cell = image->candidate_cells[0];
		if ( verbose ) {
			STATUS("--------------------\n");
			STATUS("The indexed cell (matching not performed):\n");
			cell_print(image->indexed_cell);
			STATUS("--------------------\n");
		}
		return;
	}

	for ( i=0; i<image->ncells; i++ ) {

		UnitCell *new_cell = NULL;

		if ( verbose ) {
			STATUS("--------------------\n");
			STATUS("Candidate cell %i (before matching):\n", i);
			cell_print(image->candidate_cells[i]);
			STATUS("--------------------\n");
		}

		/* Match or reduce the cell as appropriate */
		switch ( cellr ) {
		case CELLR_NONE :
			/* Never happens */
			break;
		case CELLR_REDUCE :
			new_cell = match_cell(image->candidate_cells[i],
			                      cell, verbose);
			break;
		case CELLR_COMPARE :
			if ( cells_similar(image->candidate_cells[i], cell) ) {
				new_cell = image->candidate_cells[i];
			}
			break;
		}

		image->indexed_cell = new_cell;
		if ( new_cell != NULL ) {
			STATUS("Matched on attempt %i.\n", i);
			goto done;
		}

	}

done:
	for ( i=0; i<image->ncells; i++ ) {
		cell_free(image->candidate_cells[i]);
	}
}
