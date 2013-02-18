/*
 * xds.c
 *
 * Invoke xds for crystal autoindexing
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2013 Cornelius Gati
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


/* Global private data, prepared once */
struct xds_private
{
	IndexingMethod indm;
	float *ltl;
	UnitCell *cell;
};


struct xds_data {

	/* Low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

	/* High-level stuff */
	char			newmatfile[128];
	char 			spotfile[128];
	int			step;
	int			finished_ok;
	UnitCell		*target_cell;
};

static void xds_parseline(const char *line, struct image *image,
                             struct xds_data *xds)
{

	char *copy;
	int i;

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='r';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("XDS: %s\n", copy);
	free(copy);
}


static int xds_readable(struct image *image, struct xds_data *xds)
{
	int rval;
	int no_string = 0;

	rval = read(xds->pty, xds->rbuffer+xds->rbufpos,
	            xds->rbuflen-xds->rbufpos);
	if ( (rval == -1) || (rval == 0) ) return 1;

	xds->rbufpos += rval;
	assert(xds->rbufpos <= xds->rbuflen);

	while ( (!no_string) && (xds->rbufpos > 0) ) {

		int i;
		int block_ready = 0;

		/* See if there's a full line in the buffer yet */
		for ( i=0; i<xds->rbufpos-1; i++ ) {
			/* Means the last value looked at is rbufpos-2 */

			if ( (xds->rbuffer[i] == '\r')
			  && (xds->rbuffer[i+1] == '\n') ) {
				block_ready = 1;
				break;
			}

		}

		if ( block_ready ) {

			unsigned int new_rbuflen;
			unsigned int endbit_length;
			char *block_buffer = NULL;

			block_buffer = malloc(i+1);
			memcpy(block_buffer, xds->rbuffer, i);
			block_buffer[i] = '\0';

			if ( block_buffer[0] == '\r' ) {
				memmove(block_buffer, block_buffer+1, i);
			}

			xds_parseline(block_buffer, image, xds);
			free(block_buffer);
			endbit_length = i+2;

			/* Now the block's been parsed, it should be
			 * forgotten about */
			memmove(xds->rbuffer,
			        xds->rbuffer + endbit_length,
			        xds->rbuflen - endbit_length);

			/* Subtract the number of bytes removed */
			xds->rbufpos = xds->rbufpos
			                       - endbit_length;
			new_rbuflen = xds->rbuflen - endbit_length;
			if ( new_rbuflen == 0 ) new_rbuflen = 256;
			xds->rbuffer = realloc(xds->rbuffer,
			                               new_rbuflen);
			xds->rbuflen = new_rbuflen;

		} else {

			if ( xds->rbufpos==xds->rbuflen ) {

				/* More buffer space is needed */
				xds->rbuffer = realloc(
				                    xds->rbuffer,
				                    xds->rbuflen + 256);
				xds->rbuflen = xds->rbuflen + 256;
				/* The new space gets used at the next
				 * read, shortly... */

			}
			no_string = 1;

		}

	}

	return 0;
}


static int check_cell(struct xds_private *xp, struct image *image,
                      UnitCell *cell)
{
	UnitCell *out;
	Crystal *cr;

	if ( xp->indm & INDEXING_CHECK_CELL_COMBINATIONS ) {

		out = match_cell(cell, xp->cell, 0, xp->ltl, 1);
		if ( out == NULL ) return 0;

	} else if ( xp->indm & INDEXING_CHECK_CELL_AXES ) {

		out = match_cell(cell, xp->cell, 0, xp->ltl, 0);
		if ( out == NULL ) return 0;

	} else {
		out = cell_new_from_cell(cell);
	}

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 0;
	}

	crystal_set_cell(cr, out);

	if ( xp->indm & INDEXING_CHECK_PEAKS ) {
		if ( !peak_sanity_check(image, &cr, 1) ) {
			cell_free(out);
			crystal_free(cr);
			return 0;
		}
	}

	image_add_crystal(image, cr);

	return 1;
}


static int read_cell(const char *filename, struct image *image,
                     struct xds_private *xp)
{
	FILE * fh;
	float axstar, aystar, azstar;
	float bxstar, bystar, bzstar;
	float cxstar, cystar, czstar;
	char asx[11], asy[11], asz[11];
	char bsx[11], bsy[11], bsz[11];
	char csx[11], csy[11], csz[11];
	char *rval, line[1024];
	int r;
	UnitCell *cell;

	fh = fopen("IDXREF.LP", "r");
	if ( fh == NULL ) {
		ERROR("Couldn't open 'IDXREF.LP'\n");
		return 0;
	}

	do {
		rval = fgets(line, 1023, fh);

	} while ( strcmp(line, "   #  COORDINATES OF REC. BASIS VECTOR"
	                       "    LENGTH   1/LENGTH\n") != 0 );

	/* Free line after chunk */
	rval = fgets(line, 1023, fh);
	rval = fgets(line, 1023, fh);

	memcpy(asx, line+7, 10);    asx[10] = '\0';
	memcpy(asy, line+17, 10);   asy[10] = '\0';
	memcpy(asz, line+27, 10);   asz[10] = '\0';

	rval = fgets(line, 1023, fh);

	memcpy(bsx, line+7, 10);    bsx[10] = '\0';
	memcpy(bsy, line+17, 10);   bsy[10] = '\0';
	memcpy(bsz, line+27, 10);   bsz[10] = '\0';

	rval = fgets(line, 1023, fh);

	memcpy(csx, line+7, 10);    csx[10] = '\0';
	memcpy(csy, line+17, 10);   csy[10] = '\0';
	memcpy(csz, line+27, 10);   csz[10] = '\0';

	r =  sscanf(asx, "%f", &cxstar);
	r += sscanf(asy, "%f", &cystar);
	r += sscanf(asz, "%f", &czstar);
	r += sscanf(bsx, "%f", &bxstar);
	r += sscanf(bsy, "%f", &bystar);
	r += sscanf(bsz, "%f", &bzstar);
	r += sscanf(csx, "%f", &axstar);
	r += sscanf(csy, "%f", &aystar);
	r += sscanf(csz, "%f", &azstar);

	if ( r != 9 ) {
		STATUS("Fewer than 9 parameters found in NEWMAT file.\n");
		return 0;
	}

	fclose(fh);

	//printf("asx='%s'\n", asx);
	//printf("asy='%s'\n", asy);
	//printf("asz='%s'\n", asz);
	//printf("cxstar='%f'\n", cxstar);
	//printf("cystar='%f'\n", cystar);
	//printf("czstar='%f'\n", czstar);
	//printf("bsx='%s'\n", bsx);
	//printf("bsy='%s'\n", bsy);
	//printf("bsz='%s'\n", bsz);
	//printf("axstar='%f'\n", axstar);
	//printf("aystar='%f'\n", aystar);
	//printf("azstar='%f'\n", azstar);
	//printf("csx='%s'\n", csx);
	//printf("csy='%s'\n", csy);
	//printf("csz='%s'\n", csz);
	//printf("bxstar='%f'\n", bxstar);
	//printf("bystar='%f'\n", bystar);
	//printf("bzstar='%f'\n", bzstar);

	cell = cell_new();
	cell_set_reciprocal(cell,
	                    axstar*10e9, aystar*10e9, azstar*10e9,
	                    bxstar*10e9, bystar*10e9, bzstar*10e9,
	                   -cxstar*10e9, -cystar*10e9, -czstar*10e9);
	r = check_cell(xp, image, cell);
	cell_free(cell);

	return r;
}


static void write_spot(struct image *image, const char *filename)
{
	FILE *fh;
	int i;
	double fclen = 99.0e-3;  /* fake camera length in m */
	int n;

	fh = fopen("SPOT.XDS", "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", "SPOT.XDS");
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
		rx = ((xs + p->cnx) / p->res);
		ry = ((ys + p->cny) / p->res);

		//printf("xs=%f ys=%f  ---->  rx=%f ry=%f\n", xs, ys, rx, ry);

		x = rx*fclen/p->clen;
		y = ry*fclen/p->clen;  /* Peak positions in m */

		//printf("x=%f y=%f\n", x, y);

		x = (rx / 70e-6) + 1500;
		y = (ry / 70e-6) + 1500;

		if (f->intensity <= 0) continue;

		fprintf(fh, "%10.2f %10.2f %10.2f %10.0f.\n",
		        x, y, 0.5, f->intensity);

	}
	fclose(fh);
}


static char *write_inp(struct image *image)
{
	FILE *fh;
	char *filename;

	filename = malloc(1024);
	if ( filename == NULL ) return NULL;

	snprintf(filename, 1023, "XDS.INP");

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		free(filename);
		return NULL;
	}

	fprintf(fh, "JOB= IDXREF\n");
	fprintf(fh, "ORGX= 1500\n");
	fprintf(fh, "ORGY= 1500\n");
	fprintf(fh, "DETECTOR_DISTANCE= 99.00\n"); 	//IMPORTANT
	fprintf(fh, "OSCILLATION_RANGE= 0.300\n");
	fprintf(fh, "X-RAY_WAVELENGTH= 1.32\n");        // IMPORTANT
	fprintf(fh, "NAME_TEMPLATE_OF_DATA_FRAMES=/home/dikay/insu/data_exp1/ins_ssad_1_???.img \n");
	fprintf(fh, "DATA_RANGE=1 1\n");
	fprintf(fh, "SPOT_RANGE=1 1\n");
	fprintf(fh, "SPACE_GROUP_NUMBER= 94\n"); //CatB 94
	fprintf(fh, "UNIT_CELL_CONSTANTS= 126.2 126.2 54.2 90 90 90\n"); //CatB 125.4 125.4 54.56 90 90 90
	//fprintf(fh, "SPACE_GROUP_NUMBER=194\n"); //PS1 194
	//fprintf(fh, "UNIT_CELL_CONSTANTS=281 281 165.2 90 90 120\n"); //PS1 281 281 165.2 90 90 120
	//fprintf(fh, "SPACE_GROUP_NUMBER= 0\n"); //LYS 96
	//fprintf(fh, "UNIT_CELL_CONSTANTS= 0 0 0 0 0 0\n"); //LYS 77.32 77.32 38.16 90 90 90
	fprintf(fh, "NX= 3000\n");
	fprintf(fh, "NY= 3000\n");
	fprintf(fh, "QX= 0.07\n");
	fprintf(fh, "QY= 0.07\n");
	fprintf(fh, "INDEX_ORIGIN=0 0 0\n");
	fprintf(fh, "DIRECTION_OF_DETECTOR_X-AXIS=1 0 0\n");
	fprintf(fh, "DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0\n");
	fprintf(fh, "INCIDENT_BEAM_DIRECTION=0 0 1\n");
	fprintf(fh, "ROTATION_AXIS=0 1 0\n");
	fprintf(fh, "DETECTOR= CSPAD\n");
	fprintf(fh, "MINIMUM_VALID_PIXEL_VALUE= 1\n");
	fprintf(fh, "OVERLOAD= 200000000\n");
	fprintf(fh, "INDEX_ERROR= 0.4\n");
	//fprintf(fh, "INDEX_QUALITY= 0.5\n");
	//fprintf(fh, "REFINE(IDXREF)= ALL\n");
	//fprintf(fh, "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 1\n");
	//fprintf(fh, "MAXIMUM_ERROR_OF_SPOT_POSITION= 20.0\n");

	fclose(fh);

	return filename;
}


int run_xds(struct image *image, IndexingPrivate *priv)
{
	unsigned int opts;
	int status;
	int rval;
	int n;
	struct xds_data *xds;
	char *inp_filename;
	struct xds_private *xp = (struct xds_private *)priv;

	xds = malloc(sizeof(struct xds_data));
	if ( xds == NULL ) {
		ERROR("Couldn't allocate memory for xds data.\n");
		return 0;
	}

	xds->target_cell = xp->cell;

	write_inp(image);
	inp_filename = write_inp(image);

	if ( inp_filename == NULL ) {
		ERROR("Failed to write XDS.INP file for XDS.\n");
		return 0;
	}

	n= image_feature_count(image->features);
	if (n < 25) return 0;

	snprintf(xds->spotfile, 127, "SPOT.XDS");

	write_spot(image, xds->spotfile);

	xds->pid = forkpty(&xds->pty, NULL, NULL, NULL);

	if ( xds->pid == -1 ) {
		ERROR("Failed to fork for XDS\n");
		free(xds);
		return 0;
	}
	if ( xds->pid == 0 ) {

		/* Child process: invoke XDS */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("xds", "", (char *)NULL);
		ERROR("Failed to invoke XDS.\n");
		_exit(0);

	}

	xds->rbuffer = malloc(256);
	xds->rbuflen = 256;
	xds->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(xds->pty, F_GETFL);
	fcntl(xds->pty, F_SETFL, opts | O_NONBLOCK);

	//xds->step = 1;	/* This starts the "initialisation" procedure */
	xds->finished_ok = 0;

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(xds->pty, &fds);

		tv.tv_sec = 30;
		tv.tv_usec = 0;

		sval = select(xds->pty+1, &fds, NULL, NULL, &tv);

		if ( sval == -1 ) {

			const int err = errno;

			switch ( err ) {

				case EINTR:
				STATUS("Restarting select()\n");
				break;

				default:
				ERROR("select() failed: %s\n", strerror(err));
				rval = 1;
				break;

			}

		} else if ( sval != 0 ) {
			rval = xds_readable(image, xds);
		} else {
			ERROR("No response from XDS..\n");
			rval = 1;
		}

	} while ( !rval );

	close(xds->pty);
	free(xds->rbuffer);
	waitpid(xds->pid, &status, 0);

				//FIXME//
	//if ( xds->finished_ok == 1 ) {
	//	ERROR("XDS doesn't seem to be working properly.\n");
	//} else {
		/* Read the XDS NEWMAT file and get cell if found */
	rval = read_cell(xds->newmatfile, image, xp);

	//}

	free(xds);

	return rval;
}


IndexingPrivate *xds_prepare(IndexingMethod *indm, UnitCell *cell,
                             const char *filename,
                             struct detector *det,
                             struct beam_params *beam, float *ltl)
{
	struct xds_private *xp;

	if ( cell == NULL ) {
		ERROR("XDS needs a unit cell.\n");
		return NULL;
	}

	xp = calloc(1, sizeof(*xp));
	if ( xp == NULL ) return NULL;

	/* Flags that XDS knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_CELL_COMBINATIONS
	          | INDEXING_CHECK_CELL_AXES | INDEXING_USE_LATTICE_TYPE
	          | INDEXING_CHECK_PEAKS;

	xp->ltl = ltl;
	xp->cell = cell;
	xp->indm = *indm;

	return (IndexingPrivate *)xp;
}


void xds_cleanup(IndexingPrivate *pp)
{
	struct xds_private *xp;

	xp = (struct xds_private *)pp;
	free(xp);
}
