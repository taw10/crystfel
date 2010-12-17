/*
 * index.c
 *
 * Perform indexing (somehow)
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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


IndexingPrivate **prepare_indexing(IndexingMethod *indm, UnitCell *cell,
                                  const char *filename, struct detector *det,
                                  double nominal_photon_energy)
{
	int n;
	int nm = 0;
	IndexingPrivate **iprivs;

	while ( indm[nm] != INDEXING_NONE ) nm++;
	STATUS("Preparing %i indexing methods.\n", nm);
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
		case INDEXING_TEMPLATE :
			iprivs[n] = generate_templates(cell, filename, det,
				                  nominal_photon_energy);
			break;
		}

		n++;

	}
	iprivs[n] = NULL;

	return iprivs;
}


void cleanup_indexing(IndexingPrivate **priv)
{
	int n = 0;

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
		case INDEXING_TEMPLATE :
			free_templates(priv[n]);
		}

		n++;

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


/* write .spt file for mosflm */
/* need to sort mosflm peaks by intensity... */
struct sptline {
	double x; /* x coordinate of peak */
	double y; /* y coordinate of peak */
	double h; /* height of peak */
	double s; /* sigma of peak */
};


static int compare_vals(const void *ap, const void *bp)
{
	const struct sptline a = *(struct sptline *)ap;
	const struct sptline b = *(struct sptline *)bp;

	if ( a.h < b.h ) return 1;
	if ( a.h > b.h ) return -1;
	return 0;
}


static void write_spt(struct image *image)
{
	FILE *fh;
	int i;
	char filename[1024];
	double fclen=67.8;  /* fake camera length in mm */
	double fpix=0.075;  /* fake pixel size in mm */
	double pix;
	double height=100;
	double sigma=1;
	int nPeaks = image_feature_count(image->features);

	snprintf(filename, 1023, "xfel-%i.spt", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file xfel.spt\n");
		return;
	}

	fprintf(fh, "%10d %10d %10.8f %10.6f %10.6f\n", 1, 1, fpix, 1.0, 0.0);
	fprintf(fh, "%10d %10d\n", 1, 1);
	fprintf(fh, "%10.5f %10.5f\n", 0.0, 0.0);

	struct sptline *sptlines;
	sptlines = malloc(sizeof(struct sptline)*nPeaks);

	for ( i=0; i<nPeaks; i++ ) {

		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		struct panel *pan;
		pan = find_panel(image->det,f->x,f->y);
		if ( pan == NULL ) continue;

		pix = 1000/pan->res; /* pixel size in mm */
		height = f->intensity;

		sptlines[i].x = (f->y - pan->cy)*pix*fclen/pan->clen/1000;
		sptlines[i].y = -(f->x - pan->cx)*pix*fclen/pan->clen/1000;
		sptlines[i].h = height;
		sptlines[i].s = sigma;

	}

	qsort(sptlines, nPeaks, sizeof(struct sptline), compare_vals);

	for ( i=0; i<nPeaks; i++ ) {

		fprintf(fh, "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
		        sptlines[i].x, sptlines[i].y,
		        0.0, 0.0,
		        sptlines[i].h, sptlines[i].s);

	}

	fprintf(fh,"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	           -999.0,-999.0,-999.0,-999.0,-999.0,-999.0);
	fclose(fh);
}

/* write a dummy 1x1 pixel image file for mosflm.  Without post refinement,
   mosflm will ignore this, but it must be present.*/
void write_img(struct image *image)
{
	FILE *fh;
	char filename[1024];
	unsigned short int * intimage;

	intimage = malloc(sizeof(unsigned short int));
	intimage[0] = 1;

	snprintf(filename, 1023, "xfel-%i_001.img", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file xfel.spt\n");
		return;
	}

	fprintf(fh,"{\nHEADER_BYTES=512;\n");
	fprintf(fh,"BYTE_ORDER=little_endian;\n");
	fprintf(fh,"TYPE=unsigned_short;\n");
	fprintf(fh,"DIM=2;\n");
	fprintf(fh,"SIZE1=1;\n");
	fprintf(fh,"SIZE2=1;\n");
	fprintf(fh,"}\n");
	while ( ftell(fh) < 512 ) { fprintf(fh," "); };
	fwrite(intimage,sizeof(unsigned short int),1,fh);
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


void index_pattern(struct image *image, UnitCell *cell, IndexingMethod *indm,
                   int cellr, int verbose, IndexingPrivate **ipriv)
{
	int i;
	int n = 0;

	map_all_peaks(image);

	while ( indm[n] != INDEXING_NONE ) {

		image->ncells = 0;

		/* Index as appropriate */
		switch ( indm[n] ) {
		case INDEXING_NONE :
			return;
		case INDEXING_DIRAX :
			STATUS("Running DirAx...\n");
			write_drx(image);
			run_dirax(image);
			break;
		case INDEXING_MOSFLM :
			STATUS("Running MOSFLM...\n");
			write_spt(image);
			write_img(image); /* dummy image */
			run_mosflm(image, cell);
			break;
		case INDEXING_TEMPLATE :
			match_templates(image, ipriv[n]);
			break;
		}

		if ( image->ncells == 0 ) {
			STATUS("No candidate cells found.\n");
			return;
		}

		if ( (cellr == CELLR_NONE) || (indm[n] == INDEXING_TEMPLATE) ) {
			image->indexed_cell = image->candidate_cells[0];
			if ( verbose ) {
				STATUS("--------------------\n");
				STATUS("The indexed cell (matching not"
				       " performed):\n");
				cell_print(image->indexed_cell);
				STATUS("--------------------\n");
			}
			return;
		}

		for ( i=0; i<image->ncells; i++ ) {

			UnitCell *new_cell = NULL;

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
				/* Never happens */
				break;
			case CELLR_REDUCE :
				new_cell = match_cell(image->candidate_cells[i],
					              cell, verbose, 1);
				break;
			case CELLR_COMPARE :
				new_cell = match_cell(image->candidate_cells[i],
					              cell, verbose, 0);
				break;
			}

			image->indexed_cell = new_cell;
			if ( new_cell != NULL ) {
				STATUS("Matched on attempt %i.\n", i);
				goto done;
			}

		}

		/* Move on to the next indexing method */
		n++;

	}

done:
	for ( i=0; i<image->ncells; i++ ) {
		cell_free(image->candidate_cells[i]);
	}
}


IndexingMethod *build_indexer_list(const char *str, int *need_cell)
{
	int n, i;
	char **methods;
	IndexingMethod *list;
	*need_cell = 0;

	n = assplode(str, ",", &methods, ASSPLODE_NONE);
	list = malloc((n+1)*sizeof(IndexingMethod));

	for ( i=0; i<n; i++ ) {

		if ( strcmp(methods[i], "dirax") == 0) {
			list[i] = INDEXING_DIRAX;
		} else if ( strcmp(methods[i], "mosflm") == 0) {
			list[i] = INDEXING_MOSFLM;
		} else if ( strcmp(methods[i], "template") == 0) {
			list[i] = INDEXING_TEMPLATE;
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
