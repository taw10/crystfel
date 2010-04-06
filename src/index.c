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


void index_pattern(struct image *image, UnitCell *cell, IndexingMethod indm,
                   int no_match, int verbose)
{
	int i;
	int nc = 0;

	/* Map positions to 3D */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		int c;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		c = map_position(image, f->x, f->y, &f->rx, &f->ry, &f->rz);
		if ( c != 0 ) nc++;

	}
	if ( nc ) {
		ERROR("Failed to map %i reflections\n", nc);
	}

	write_drx(image);

	image->ncells = 0;

	/* Index (or not) as appropriate */
	if ( indm == INDEXING_NONE ) return;
	if ( indm == INDEXING_DIRAX ) run_dirax(image);

	if ( image->ncells == 0 ) {
		STATUS("No candidate cells found.\n");
		return;
	}

	if ( no_match ) {
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

		new_cell = match_cell(image->candidate_cells[i], cell, verbose);
		image->indexed_cell = new_cell;
		if ( new_cell != NULL ) {
			STATUS("Matched on attempt %i.\n", i);
			goto done;
		}

	}

done:
	for ( i=0; i<image->ncells; i++ ) {
		free(image->candidate_cells[i]);
	}
}
