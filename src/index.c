/*
 * index.c
 *
 * Perform indexing (somehow)
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
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


/* x,y in pixels relative to central beam */
static int map_position(struct image *image, double x, double y,
                        double *rx, double *ry, double *rz)
{
	/* "Input" space */
	double d;

	/* Angular description of reflection */
	double theta, psi, k;

	k = 1.0 / image->lambda;

	if ( image->fmode == FORMULATION_CLEN ) {

		/* Convert pixels to metres */
		x /= image->resolution;
		y /= image->resolution;	/* Convert pixels to metres */
		d = sqrt((x*x) + (y*y));
		theta = atan2(d, image->camera_len);

	} else if (image->fmode == FORMULATION_PIXELSIZE ) {

		/* Convert pixels to metres^-1 */
		x *= image->pixel_size;
		y *= image->pixel_size;	/* Convert pixels to metres^-1 */
		d = sqrt((x*x) + (y*y));
		theta = atan2(d, k);

	} else {
		ERROR("Unrecognised formulation mode in mapping_scale.\n");
		return -1;
	}

	psi = atan2(y, x);

	*rx = k*sin(theta)*cos(psi);
	*ry = k*sin(theta)*sin(psi);
	*rz = k - k*cos(theta);

	return 0;
}


void index_pattern(struct image *image, int no_index, int dump_peaks,
                   int use_dirax)
{
	int i;

	/* Perform 'fine' peak search */
	search_peaks(image, dump_peaks);

	/* Map positions to 3D */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;

		f = image_get_feature(image->features, i);

		if ( f->y >=512 ) {
			/* Top half of CCD */
			map_position(image, f->x-UPPER_CX, f->y-UPPER_CY,
			             &f->rx, &f->ry, &f->rz);
		} else {
			/* Lower half of CCD */
			map_position(image, f->x-LOWER_CX, f->y-LOWER_CY,
			             &f->rx, &f->ry, &f->rz);
		}
	}

	if ( use_dirax ) {
		run_dirax(image, no_index);
		return;
	}
}
