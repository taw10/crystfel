/*
 * xds.c
 *
 * Invoke xds for crystal autoindexing
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2013 Cornelius Gati
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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
#include <pty.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/ioctl.h>
#include <errno.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "cell.h"
#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "detector.h"
#include "cell-utils.h"


#define XDS_VERBOSE 0


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
	int			step;
	int			finished_ok;
	UnitCell		*target_cell;
};

static void xds_parseline(const char *line, struct image *image,
                             struct xds_data *xds)
{
#if XDS_VERBOSE
	char *copy;
	int i;

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='r';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("XDS: %s\n", copy);
	free(copy);
#endif
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


static int read_cell(struct image *image, struct xds_private *xp)
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
		if ( rval == NULL ) return 0;

	} while ( strcmp(line, "   #  COORDINATES OF REC. BASIS VECTOR"
	                       "    LENGTH   1/LENGTH\n") != 0 );

	/* Free line after chunk */
	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) return 0;
	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) return 0;

	memcpy(asx, line+7, 10);    asx[10] = '\0';
	memcpy(asy, line+17, 10);   asy[10] = '\0';
	memcpy(asz, line+27, 10);   asz[10] = '\0';

	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) return 0;

	memcpy(bsx, line+7, 10);    bsx[10] = '\0';
	memcpy(bsy, line+17, 10);   bsy[10] = '\0';
	memcpy(bsz, line+27, 10);   bsz[10] = '\0';

	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) return 0;

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

	cell = cell_new();
	cell_set_reciprocal(cell,
	                    axstar*10e9, aystar*10e9, azstar*10e9,
	                    bxstar*10e9, bystar*10e9, bzstar*10e9,
	                   -cxstar*10e9, -cystar*10e9, -czstar*10e9);
	r = check_cell(xp, image, cell);
	cell_free(cell);

	return r;
}


static void write_spot(struct image *image)
{
	FILE *fh;
	int i;
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

		x = rx;
		y = ry; /* Peak positions in m */

		//printf("x=%f y=%f\n", x, y);

		x = (rx / 70e-6) + 1500;
		y = (ry / 70e-6) + 1500;

		if (f->intensity <= 0) continue;

		fprintf(fh, "%10.2f %10.2f %10.2f %10.0f.\n",
		        x, y, 0.5, f->intensity);

	}
	fclose(fh);
}


/* Turn what we know about the unit cell into something which we can give to
 * XDS to make it give us only indexing results compatible with the cell. */
static const char *spacegroup_for_lattice(UnitCell *cell)
{
	LatticeType latt;
	char centering;
	char *g = NULL;

	latt = cell_get_lattice_type(cell);
	centering = cell_get_centering(cell);

	switch ( latt )
	{
		case L_TRICLINIC :
		g = "1";
		break;

		case L_MONOCLINIC :
		if ( centering == 'P' )	{
			g = "3";
		} else {
			g = "5";
		}
		break;

		case L_ORTHORHOMBIC :
		if ( centering == 'P' ) {
			g = "16";
		} else if ( centering == 'C' ) {
			g = "20";
		} else if ( centering == 'F' ) {
			g = "22";
		} else {
			g = "23";
		}
		break;

		case L_TETRAGONAL :
		if ( centering == 'P' ) {
			g = "75";
		} else {
			g = "79";
		}
		break;

		case L_RHOMBOHEDRAL :
		if ( centering == 'P' ) {
			g = "143";
		} else {
			g = "146";
		}
		break;

		case L_HEXAGONAL :
		g = "168";
		break;

		case L_CUBIC :
		if ( centering == 'P' ) {
			g = "195";
		} else if ( centering == 'F' ) {
			g = "196";
		} else {
			g = "197";
		}
		break;
	}
	assert(g != NULL);

	return g;
}


static int write_inp(struct image *image, struct xds_private *xp)
{
	FILE *fh;

	fh = fopen("XDS.INP", "w");
	if ( !fh ) {
		ERROR("Couldn't open XDS.INP\n");
		return 1;
	}

	fprintf(fh, "JOB= IDXREF\n");
	fprintf(fh, "ORGX= 1500\n");
	fprintf(fh, "ORGY= 1500\n");
	fprintf(fh, "DETECTOR_DISTANCE= %f\n", image->det->panels[0].clen*1000);
	fprintf(fh, "OSCILLATION_RANGE= 0.300\n");
	fprintf(fh, "X-RAY_WAVELENGTH= %.6f\n", image->lambda*1e10);
	fprintf(fh, "NAME_TEMPLATE_OF_DATA_FRAMES=???.img \n");
	fprintf(fh, "DATA_RANGE=1 1\n");
	fprintf(fh, "SPOT_RANGE=1 1\n");

	if ( xp->indm & INDEXING_USE_LATTICE_TYPE ) {
		fprintf(fh, "SPACE_GROUP_NUMBER= %s\n",
			    spacegroup_for_lattice(xp->cell));
	} else {
		fprintf(fh, "SPACE_GROUP_NUMBER= 0\n");
	}

	if ( xp->indm & INDEXING_USE_CELL_PARAMETERS ) {

		double a, b, c, al, be, ga;

		cell_get_parameters(xp->cell, &a, &b, &c, &al, &be, &ga);

		fprintf(fh, "UNIT_CELL_CONSTANTS= "
		            "%.2f %.2f %.2f %.2f %.2f %.2f\n",
		            a*1e10, b*1e10, c*1e10,
		            rad2deg(al), rad2deg(be), rad2deg(ga));

        } else {
		fprintf(fh, "UNIT_CELL_CONSTANTS= 0 0 0 0 0 0\n");
	}

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
	fprintf(fh, "INDEX_ERROR= 0.05\n");
	//fprintf(fh, "INDEX_QUALITY= 0.5\n");
	fprintf(fh, "REFINE(IDXREF)= CELL ORIENTATION\n");
	//fprintf(fh, "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= 1\n");
	//fprintf(fh, "MAXIMUM_ERROR_OF_SPOT_POSITION= 20.0\n");

	fclose(fh);

	return 0;
}


int run_xds(struct image *image, IndexingPrivate *priv)
{
	unsigned int opts;
	int status;
	int rval;
	int n;
	struct xds_data *xds;
	struct xds_private *xp = (struct xds_private *)priv;

	xds = malloc(sizeof(struct xds_data));
	if ( xds == NULL ) {
		ERROR("Couldn't allocate memory for xds data.\n");
		return 0;
	}

	xds->target_cell = xp->cell;

	if ( write_inp(image, xp) ) {
		ERROR("Failed to write XDS.INP file for XDS.\n");
		return 0;
	}

	n = image_feature_count(image->features);
	if (n < 25) return 0;

	write_spot(image);

	/* Delete any old indexing result which may exist */
	remove("IDXREF.LP");

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
				rval = 0;
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

	rval = read_cell(image, xp);

	free(xds);

	return rval;
}


IndexingPrivate *xds_prepare(IndexingMethod *indm, UnitCell *cell,
                             struct detector *det, float *ltl)
{
	struct xds_private *xp;

	/* Either cell,latt and cell provided, or nocell-nolatt and no cell
	 * - complain about anything else.  Could figure this out automatically,
	 * but we'd have to decide whether the user just forgot the cell, or
	 * forgot "-nolatt", or whatever. */
	if ( ((*indm & INDEXING_USE_LATTICE_TYPE)
	  || (*indm & INDEXING_USE_CELL_PARAMETERS))
	  && !cell_has_parameters(cell) )
	{
		ERROR("No cell parameters provided.  If you wanted to use XDS "
		      "without prior cell information, use "
		      "xds-nolatt-nocell.\n");
		return NULL;
	}

	if ( (*indm & INDEXING_USE_LATTICE_TYPE)
	  && !(*indm & INDEXING_USE_CELL_PARAMETERS) ) {
		ERROR("Invalid XDS options (-latt-nocell): "
		      "try xds-nolatt-nocell.\n");
		return NULL;
	}

	if ( (*indm & INDEXING_USE_CELL_PARAMETERS)
	  && !(*indm & INDEXING_USE_LATTICE_TYPE) ) {
		ERROR("Invalid XDS options (-cell-nolatt): "
		      "try xds-nolatt-nocell.\n");
		return NULL;
	}

	if ( ((*indm & INDEXING_USE_CELL_PARAMETERS)
	  || (*indm & INDEXING_USE_LATTICE_TYPE))
	  && !(*indm & INDEXING_CHECK_CELL_AXES)
	  && !(*indm & INDEXING_CHECK_CELL_COMBINATIONS) ) {
		ERROR("The cell from xds-raw-cell or xds-raw-latt may have had"
		      " its axes permuted from the cell you provided.  If this"
		      " is a problem, consider using xds-axes-cell.\n");
	}

	xp = calloc(1, sizeof(*xp));
	if ( xp == NULL ) return NULL;

	/* Flags that XDS knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_CELL_COMBINATIONS
	          | INDEXING_CHECK_CELL_AXES | INDEXING_USE_LATTICE_TYPE
	          | INDEXING_CHECK_PEAKS | INDEXING_USE_CELL_PARAMETERS;

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
