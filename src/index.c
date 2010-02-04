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


/* x,y in pixels relative to central beam */
int map_position(struct image *image, double x, double y,
                 double *rx, double *ry, double *rz)
{
	/* "Input" space */
	double d;

	/* Angular description of reflection */
	double twotheta, psi, k;

	k = 1.0 / image->lambda;

	if ( image->fmode == FORMULATION_CLEN ) {

		/* Convert pixels to metres */
		x /= image->resolution;
		y /= image->resolution;	/* Convert pixels to metres */
		d = sqrt((x*x) + (y*y));
		twotheta = atan2(d, image->camera_len);

	} else if (image->fmode == FORMULATION_PIXELSIZE ) {

		/* Convert pixels to metres^-1 */
		x *= image->pixel_size;
		y *= image->pixel_size;	/* Convert pixels to metres^-1 */
		d = sqrt((x*x) + (y*y));
		twotheta = atan2(d, k);

	} else {
		ERROR("Unrecognised formulation mode in mapping_scale.\n");
		return -1;
	}

	psi = atan2(y, x);

	*rx = k*sin(twotheta)*cos(psi);
	*ry = k*sin(twotheta)*sin(psi);
	*rz = k - k*cos(twotheta);

	return 0;
}


static void write_drx(struct image *image)
{
	FILE *fh;
	int i;

	STATUS("Writing xfel.drx file.  Remember that it uses units of "
	       "reciprocal Angstroms!\n");

	fh = fopen("xfel.drx", "w");
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


void index_pattern(struct image *image, IndexingMethod indm)
{
	int i;
	UnitCell *new_cell = NULL;

	/* Map positions to 3D */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double rx = 0.0;
		double ry = 0.0;
		int p;
		int found = 0;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		for ( p=0; p<image->det.n_panels; p++ ) {
			if ( (f->x >= image->det.panels[p].min_x)
			  && (f->x <= image->det.panels[p].max_x)
			  && (f->y >= image->det.panels[p].min_y)
			  && (f->y <= image->det.panels[p].max_y) ) {
				rx = ((double)f->x - image->det.panels[p].cx);
				ry = ((double)f->y - image->det.panels[p].cy);
				found = 1;
			}
		}
		if ( !found ) {
			ERROR("No mapping found for %f,%f\n", f->x, f->y);
			continue;
		}

		map_position(image, rx, ry, &f->rx, &f->ry, &f->rz);
	}

	write_drx(image);

	/* Index (or not) as appropriate */
	if ( indm == INDEXING_NONE ) return;
	if ( indm == INDEXING_DIRAX ) run_dirax(image);

	new_cell = match_cell(image->indexed_cell,
		              image->molecule->cell);
	free(image->indexed_cell);
	image->indexed_cell = new_cell;
}
