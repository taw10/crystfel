/*
 * xds.c
 *
 * Invoke xds for crystal autoindexing
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
 *        2013 Cornelius Gati <cornelius.gati@cfel.de>
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
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/ioctl.h>
#include <errno.h>

#if HAVE_FORKPTY_LINUX
#include <pty.h>
#elif HAVE_FORKPTY_BSD
#include <util.h>
#endif


#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell.h"
#include "cell-utils.h"


#define xds_VERBOSE 0


struct xds_data {

	/* Low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

};

static void write_spot(struct image *image, const char *filename)
{
	FILE *fh;
	int i;
	double fclen = 67.8e-3;  /* fake camera length in m */
	int n;

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}

	n = image_feature_count(image->features);
	for ( i=0; i<n; i++ ) 
        
        {
		struct imagefeature *f;
		struct panel *p;
		double xs, ys, rx, ry, x, y;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		p = find_panel(image->det, f->fs, f->ss);
		if ( p == NULL ) continue;

		xs = (f->fs-p->min_fs)*p->fsx + (f->ss-p->min_ss)*p->ssx;
		ys = (f->fs-p->min_fs)*p->fsy + (f->ss-p->min_ss)*p->ssy;
		rx = (xs + p->cnx) / p->res;
		ry = (ys + p->cny) / p->res;

		x = rx*fclen/p->clen;
		y = ry*fclen/p->clen;  /* Peak positions in m */

		fprintf(fh, "%10.2f %10.2f %10.2f %10.0f",
		        x*1e3, y*1e3, 0.0, 0.0);

	}
	fclose(fh);
}


static char *write_INP(struct image *image)
{
	FILE *fh;
	char *filename;
	//double tt;

	filename = malloc(1024);
	if ( filename == NULL ) return NULL;

	snprintf(filename, 1023, filename, image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		free(filename);
		return NULL;
	}

	//get_q_for_panel(image->det->furthest_out_panel,
	//                    image->det->furthest_out_fs,
	//                    image->det->furthest_out_ss,
	//                    &tt, 1.0/image->lambda);

	fprintf(fh, "JOB= IDXREF\n");
	fprintf(fh, "ORGX= 1540.64\n");
	fprintf(fh, "ORGY= 1543.78\n");
	fprintf(fh, "DETECTOR_DISTANCE= 50.000\n");
	fprintf(fh, "OSCILLATION_RANGE= 0.100\n");
	fprintf(fh, "X-RAY_WAVELENGTH= 1.77000\n");
	fprintf(fh, "NAME_TEMPLATE_OF_DATA_FRAMES=/home/dikay/insu/data_exp1/ins_ssad_1_???.img \n");
	fprintf(fh, "DATA_RANGE=1 1\n");
	fprintf(fh, "SPOT_RANGE=1 1\n");
	fprintf(fh, "SPACE_GROUP_NUMBER=0\n");
	fprintf(fh, "UNIT_CELL_CONSTANTS= 70 80 90 90 90 90\n");
	fprintf(fh, "NX= 3072\n");
	fprintf(fh, "NY= 3072\n");
	fprintf(fh, "QX= .073242\n");
	fprintf(fh, "QY= .073242\n");
	fprintf(fh, "DIRECTION_OF_DETECTOR_X-AXIS=1 0 0\n");
	fprintf(fh, "DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0\n");
	fprintf(fh, "INCIDENT_BEAM_DIRECTION=0 0 1\n");
	fprintf(fh, "ROTATION_AXIS=1 0 0\n");
	fprintf(fh, "DETECTOR= CCDCHESS\n");
	fprintf(fh, "MINIMUM_VALID_PIXEL_VALUE= 1\n");
	fprintf(fh, "OVERLOAD= 65500\n");

	fclose(fh);

	return filename;

void run_xds(struct image *image, UnitCell *cell)
{
	//unsigned int opts;
	//int status;
	//int rval;
	struct xds_data *xds;
	char *INP_filename;
	//char gff_filename[1024];

	write_INP(image);
	INP_filename = write_INP(image);

	if ( INP_filename == NULL ) {
		ERROR("Failed to write XDS.INP file for XDS.\n");
		return;
	}

	xds = malloc(sizeof(struct xds_data));
	if ( xds == NULL ) {
		ERROR("Couldn't allocate memory for xds data.\n");
		return;
	}}
	//free(xds);
//}
