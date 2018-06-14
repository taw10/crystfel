/*
 * xds.c
 *
 * Invoke xds for crystal autoindexing
 *
 * Copyright © 2013-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2013 Cornelius Gati
 *
 * Authors:
 *   2010-2018 Thomas White <taw@physics.org>
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

#ifdef HAVE_FORKPTY_PTY_H
#include <pty.h>
#endif
#ifdef HAVE_FORKPTY_UTIL_H
#include <util.h>
#endif

#include "xds.h"
#include "cell.h"
#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "detector.h"
#include "cell-utils.h"


/* Fake pixel size and camera length, both in metres */
#define FAKE_PIXEL_SIZE (70e-6)
#define FAKE_CLEN (0.1)


/* Global private data, prepared once */
struct xds_private
{
	IndexingMethod indm;
	UnitCell *cell;
};


static int read_cell(struct image *image)
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
	Crystal *cr;

	fh = fopen("IDXREF.LP", "r");
	if ( fh == NULL ) return 0; /* Not indexable */

	do {
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) {
			fclose(fh);
			return 0;
		}

	} while ( strcmp(line, "   #  COORDINATES OF REC. BASIS VECTOR"
	                       "    LENGTH   1/LENGTH\n") != 0 );

	/* Free line after chunk */
	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) {
		fclose(fh);
		return 0;
	}

	/* Get first vector */
	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) {
		fclose(fh);
		return 0;
	}
	if ( line[4] != '1' ) {
		ERROR("No first vector from XDS.\n");
		return 0;
	}
	memcpy(asx, line+7, 10);    asx[10] = '\0';
	memcpy(asy, line+17, 10);   asy[10] = '\0';
	memcpy(asz, line+27, 10);   asz[10] = '\0';

	/* Get second vector */
	rval = fgets(line, 1023, fh);
	if ( rval == NULL ) {
		fclose(fh);
		return 0;
	}
	if ( line[4] != '2' ) {
		ERROR("No second vector from XDS.\n");
		return 0;
	}
	memcpy(bsx, line+7, 10);    bsx[10] = '\0';
	memcpy(bsy, line+17, 10);   bsy[10] = '\0';
	memcpy(bsz, line+27, 10);   bsz[10] = '\0';

	/* Get third vector */
	rval = fgets(line, 1023, fh);
	fclose(fh);
	if ( rval == NULL ) return 0;
	if ( line[4] != '3' ) return 0; /* No error message this time
	                                 * - happens a lot */
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

	cell = cell_new();
	cell_set_reciprocal(cell,
	                    axstar*10e9, aystar*10e9, azstar*10e9,
	                    bxstar*10e9, bystar*10e9, bzstar*10e9,
	                   -cxstar*10e9, -cystar*10e9, -czstar*10e9);

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 0;
	}
	crystal_set_cell(cr, cell);
	image_add_crystal(image, cr);

	return 1;
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
		double ttx, tty, x, y;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		if ( f->intensity <= 0 ) continue;

		ttx = angle_between_2d(0.0, 1.0,
		                       f->rx, 1.0/image->lambda + f->rz);
		tty = angle_between_2d(0.0, 1.0,
		                       f->ry, 1.0/image->lambda + f->rz);
		if ( f->rx < 0.0 ) ttx *= -1.0;
		if ( f->ry < 0.0 ) tty *= -1.0;
		x = tan(ttx)*FAKE_CLEN;
		y = tan(tty)*FAKE_CLEN;

		x = (x / FAKE_PIXEL_SIZE) + 1500;
		y = (y / FAKE_PIXEL_SIZE) + 1500;

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
	fprintf(fh, "DETECTOR_DISTANCE= %f\n", FAKE_CLEN*1e3);
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
	fprintf(fh, "QX= %f\n", FAKE_PIXEL_SIZE*1e3);  /* Pixel size in mm */
	fprintf(fh, "QY= %f\n", FAKE_PIXEL_SIZE*1e3);  /* Pixel size in mm */
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


int run_xds(struct image *image, void *priv)
{
	int status;
	int rval;
	int n;
	pid_t pid;
	int pty;
	struct xds_private *xp = (struct xds_private *)priv;

	if ( write_inp(image, xp) ) {
		ERROR("Failed to write XDS.INP file for XDS.\n");
		return 0;
	}

	n = image_feature_count(image->features);
	if ( n < 25 ) return 0;

	write_spot(image);

	/* Delete any old indexing result which may exist */
	remove("IDXREF.LP");

	pid = forkpty(&pty, NULL, NULL, NULL);

	if ( pid == -1 ) {
		ERROR("Failed to fork for XDS\n");
		return 0;
	}
	if ( pid == 0 ) {

		/* Child process: invoke XDS */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("xds", "xds", (char *)NULL);
		ERROR("Failed to invoke XDS.\n");
		_exit(0);

	}
	waitpid(pid, &status, 0);

	close(pty);
	rval = read_cell(image);

	return rval;
}


void *xds_prepare(IndexingMethod *indm, UnitCell *cell)
{
	struct xds_private *xp;

	if ( xds_probe(cell) == NULL ) {
		ERROR("XDS does not appear to run properly.\n");
		ERROR("Please check your XDS installation.\n");
		return NULL;
	}

	/* Either cell,latt and cell provided, or nocell-nolatt and no cell
	 * - complain about anything else.  Could figure this out automatically,
	 * but we'd have to decide whether the user just forgot the cell, or
	 * forgot "-nolatt", or whatever. */
	if ( (*indm & INDEXING_USE_LATTICE_TYPE)
	  && !(*indm & INDEXING_USE_CELL_PARAMETERS) )
	{
		ERROR("Invalid XDS options (-latt-nocell): "
		      "try xds-nolatt-nocell.\n");
		return NULL;
	}

	if ( (*indm & INDEXING_USE_CELL_PARAMETERS)
	  && !(*indm & INDEXING_USE_LATTICE_TYPE) )
	{
		ERROR("Invalid XDS options (-cell-nolatt): "
		      "try xds-nolatt-nocell.\n");
		return NULL;
	}

	xp = calloc(1, sizeof(*xp));
	if ( xp == NULL ) return NULL;

	/* Flags that XDS knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_USE_LATTICE_TYPE
	          | INDEXING_USE_CELL_PARAMETERS;

	xp->cell = cell;
	xp->indm = *indm;

	return xp;
}


void xds_cleanup(void *pp)
{
	struct xds_private *xp;

	xp = (struct xds_private *)pp;
	free(xp);
}


const char *xds_probe(UnitCell *cell)
{
	pid_t pid;
	int pty;
	int status;
	FILE *fh;
	char line[1024];
	int ok = 0;
	int l;

	pid = forkpty(&pty, NULL, NULL, NULL);
	if ( pid == -1 ) {
		return NULL;
	}
	if ( pid == 0 ) {

		/* Child process */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("xds", "xds", (char *)NULL);
		_exit(1);

	}

	fh = fdopen(pty, "r");

	for ( l=0; l<10; l++ ) {
		char *pos;
		if ( fgets(line, 1024, fh) != NULL ) {
			pos = strstr(line, "** XDS **");
			if ( pos != NULL ) {
				ok = 1;
			}
		}
	}

	fclose(fh);
	close(pty);
	waitpid(pid, &status, 0);

	if ( !ok ) return NULL;

	if ( cell_has_parameters(cell) ) return "xds-cell-latt";
	return "xds-nocell-nolatt";
}
