/*
 * mosflm.c
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010      Richard Kirian <rkirian@asu.edu>
 *   2010-2013 Thomas White <taw@physics.org>
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

/* TODO:
 *
 * Properly read the newmat file (don't use fscanf-- spaces between numers
 *  are not guaranteed)
 *
 * "Success" is indicated by existence of NEWMAT file written by mosflm.
 *  Better to interact with mosflm directly in order to somehow verify success.
 *
 * Investigate how these keywords affect mosflms behavior:
 *
 *  MOSAICITY
 *  DISPERSION
 *  DIVERGENCE
 *  POLARISATION
 *  POSTREF BEAM
 *  POSTREF USEBEAM OFF
 *  PREREFINE ON
 *  EXTRA ON
 *  POSTREF ON
 *
 *  These did not seem to affect the results by my (Rick's) experience, probably
 *  because they are only used conjunction with image intensity data, but it's
 *  worth another look at the documentation.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pty.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "image.h"
#include "mosflm.h"
#include "utils.h"
#include "peaks.h"
#include "cell-utils.h"


#define MOSFLM_VERBOSE 0


typedef enum {
	MOSFLM_INPUT_NONE,
	MOSFLM_INPUT_LINE,
	MOSFLM_INPUT_PROMPT
} MOSFLMInputType;



struct mosflm_private {
	IndexingMethod          indm;
	float                   *ltl;
	UnitCell                *template;
};


struct mosflm_data {

	/* MOSFLM auto-indexing low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

	/* MOSFLM high-level stuff */
	char                    newmatfile[128];
	char                    imagefile[128];
	char                    sptfile[128];
	int                     step;
	int                     finished_ok;
	int                     done;
	int                     success;

	struct mosflm_private  *mp;

};

static int check_cell(struct mosflm_private *mp, struct image *image,
                      UnitCell *cell)
{
	UnitCell *out;
	Crystal *cr;

	/* If we sent lattice information, make sure that we got back what we
	 * asked for, not (e.g.) some "H" version of a rhombohedral R cell */
	if ( mp->indm & INDEXING_USE_LATTICE_TYPE ) {

		LatticeType latt_m, latt_r;
		char cen_m, cen_r;

		/* What we asked for */
		latt_r = cell_get_lattice_type(mp->template);
		cen_r = cell_get_centering(mp->template);

		/* What we got back */
		latt_m = cell_get_lattice_type(cell);
		cen_m = cell_get_centering(cell);

		if ( latt_r != latt_m ) {
			ERROR("Lattice type produced by MOSFLM (%i) does not "
			      "match what was requested (%i).  "
			      "Please report this.\n", latt_m, latt_r);
			return 0;
		}

		if ( (latt_m != L_MONOCLINIC) && (cen_r != cen_m) ) {
			ERROR("Centering produced by MOSFLM (%c) does not "
			      "match what was requested (%c).  "
			      "Please report this.\n", cen_m, cen_r);
			return 0;
		}
		/* If it's monoclinic, see the warning in mosflm_prepare() */

	}

	if ( mp->indm & INDEXING_CHECK_CELL_COMBINATIONS ) {

		out = match_cell(cell, mp->template, 0, mp->ltl, 1);
		if ( out == NULL ) return 0;

	} else if ( mp->indm & INDEXING_CHECK_CELL_AXES ) {

		out = match_cell(cell, mp->template, 0, mp->ltl, 0);
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

	if ( mp->indm & INDEXING_CHECK_PEAKS ) {
		if ( !peak_sanity_check(image, &cr, 1) ) {
			cell_free(out);
			crystal_free(cr);
			return 0;
		}
	}

	image_add_crystal(image, cr);

	return 1;
}


static void mosflm_parseline(const char *line, struct image *image,
                             struct mosflm_data *dirax)
{
	#if MOSFLM_VERBOSE
	char *copy;
	int i;

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='r';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("MOSFLM: %s\n", copy);
	free(copy);
	#endif
}


/* This is the opposite of spacegroup_for_lattice() below.
 * Note that this is not general, just a set of rules for interpreting MOSFLM's
 * output. */
static LatticeType spacegroup_to_lattice(const char *sg, char *ua, char *cen)
{
	LatticeType latt;

	*cen = sg[0];

	if ( sg[1] == '1' ) {
		latt = L_TRICLINIC;
		*ua = '*';
	} else if ( strncmp(sg+1, "23", 2) == 0 ) {
		latt = L_CUBIC;
		*ua = '*';
	} else if ( strncmp(sg+1, "222", 3) == 0 ) {
		latt = L_ORTHORHOMBIC;
		*ua = '*';
	} else if ( sg[1] == '2' ) {
		latt = L_MONOCLINIC;
		*ua = 'b';
	} else if ( sg[1] == '4' ) {
		latt = L_TETRAGONAL;
		*ua = 'c';
	} else if ( sg[1] == '6' ) {
		latt = L_HEXAGONAL;
		*ua = 'c';
	} else if ( sg[1] == '3' ) {
		if ( sg[0] == 'H' ) {
			latt = L_HEXAGONAL;
			*ua = 'c';
		} else {
			latt = L_RHOMBOHEDRAL;
			*ua = '*';
		}
	} else {
		ERROR("Couldn't understand '%s'\n", sg);
		latt = L_TRICLINIC;
	}

	return latt;
}


static int read_newmat(struct mosflm_data *mosflm, const char *filename,
                       struct image *image)
{
	FILE *fh;
	float asx, asy, asz;
	float bsx, bsy, bsz;
	float csx, csy, csz;
	int n;
	double c;
	UnitCell *cell;
	char symm[32];
	char *rval;
	int i;
	char cen;
	LatticeType latt;
	char ua = '?';

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		return 1;
	}
	n  = fscanf(fh, "%f %f %f\n", &asx, &bsx, &csx);
	n += fscanf(fh, "%f %f %f\n", &asy, &bsy, &csy);
	n += fscanf(fh, "%f %f %f\n", &asz, &bsz, &csz);
	if ( n != 9 ) {
		STATUS("Fewer than 9 parameters found in NEWMAT file.\n");
		return 1;
	}

	/* Skip the next six lines */
	for ( i=0; i<6; i++ ) {
		char tmp[1024];
		rval = fgets(tmp, 1024, fh);
		if ( rval == NULL ) {
			ERROR("Failed to read newmat file.\n");
			return 1;
		}
	}

	rval = fgets(symm, 32, fh);
	if ( rval == NULL ) {
		ERROR("Failed to read newmat file.\n");
		return 1;
	}

	fclose(fh);

	chomp(symm);
	if ( strncmp(symm, "SYMM ", 5) != 0 ) {
		ERROR("Bad 'SYMM' line from MOSFLM.\n");
		return 1;
	}
	//STATUS("MOSFLM says '%s'\n", symm);
	latt = spacegroup_to_lattice(symm+5, &ua, &cen);

	/* MOSFLM "A" matrix is multiplied by lambda, so fix this */
	c = 1.0/image->lambda;

	cell = cell_new();

	/* The relationship between the coordinates in the spot file and the
	 * resulting matrix is diabolically complicated.  This transformation
	 * seems to work, but is not derived by working through all the
	 * transformations. */
	cell_set_reciprocal(cell,
	                    -asy*c, -asz*c, asx*c,
	                    -bsy*c, -bsz*c, bsx*c,
	                    -csy*c, -csz*c, csx*c);
	cell_set_centering(cell, cen);
	cell_set_lattice_type(cell, latt);
	cell_set_unique_axis(cell, ua);
	//STATUS("My cell:\n");
	//cell_print(cell);

	if ( check_cell(mosflm->mp, image, cell) ) {
		mosflm->success = 1;
		mosflm->done = 1;
	}

	cell_free(cell);

        return 0;
}


/* Write .spt file for mosflm */
static void write_spt(struct image *image, const char *filename)
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

	/* Number of pixels in x, number of pixels in y, pixel size (mm),
	 * YSCALE, OMEGA */
	fprintf(fh, "%10d %10d %10.8f %10.6f %10.6f\n", 1, 1, 0.0, 1.0, 0.0);

	/* INVERTX, ISWUNG */
	fprintf(fh, "%10d %10d\n", 0, 1);

	/* XBEAM, YBEAM */
	fprintf(fh, "%10.5f %10.5f\n", 0.0, 0.0);

	n = image_feature_count(image->features);
	for ( i=0; i<n; i++ ) {

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

		x = -rx*fclen/p->clen;
		y = ry*fclen/p->clen;  /* Peak positions in m */

		fprintf(fh, "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
		        x*1e3, y*1e3, 0.0, 0.0, 1000.0, 10.0);

	}

	fprintf(fh,"%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	           -999.0,-999.0,-999.0,-999.0,-999.0,-999.0);
	fclose(fh);
}


/* Write a dummy 1x1 pixel image file for MOSFLM.  Without post refinement,
 * MOSFLM will ignore this, but it must be present. */
static void write_img(struct image *image, const char *filename)
{
	FILE *fh;
	unsigned short int *intimage;

	intimage = malloc(sizeof(unsigned short int));
	intimage[0] = 1;

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}

	/* Write header */
	fprintf(fh, "{\nHEADER_BYTES=512;\n");
	fprintf(fh, "BYTE_ORDER=little_endian;\n");
	fprintf(fh, "TYPE=unsigned_short;\n");
	fprintf(fh, "DIM=2;\n");
	fprintf(fh, "SIZE1=1;\n");
	fprintf(fh, "SIZE2=1;\n");
	fprintf(fh, "}\n");

	/* Header padding */
	while ( ftell(fh) < 512 ) fprintf(fh," ");

	fwrite(intimage, sizeof(unsigned short int), 1, fh);
	free(intimage);
	fclose(fh);
}


static void mosflm_sendline(const char *line, struct mosflm_data *mosflm)
{
	#if MOSFLM_VERBOSE
	char *copy;
	int i;

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='\0';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("To MOSFLM: '%s'\n", copy);
	free(copy);
	#endif

	write(mosflm->pty, line, strlen(line));
}


/* Turn what we know about the unit cell into something which we can give to
 * MOSFLM to make it give us only indexing results compatible with the cell. */
static char *spacegroup_for_lattice(UnitCell *cell)
{
	LatticeType latt;
	char centering;
	char *g = NULL;
	char *result;
	char ua;

	latt = cell_get_lattice_type(cell);
	centering = cell_get_centering(cell);
	ua = cell_get_unique_axis(cell);

	switch ( latt )
	{
		case L_TRICLINIC :
		g = "1";
		break;

		case L_MONOCLINIC :
		switch ( ua ) {
			case 'a' : g = "211"; break;
			case 'b' : g = "121"; break;
			case 'c' : g = "112"; break;
		}
		break;

		case L_ORTHORHOMBIC :
		g = "222";
		break;

		case L_TETRAGONAL :
		g = "4";
		break;

		case L_RHOMBOHEDRAL :
		g = "3";
		break;

		case L_HEXAGONAL :
		if ( centering != 'H' ) {
			g = "6";
		} else {
			g = "3";
		}
		break;

		case L_CUBIC :
		g = "23";
		break;
	}
	assert(g != NULL);

	result = malloc(32);
	if ( result == NULL ) return NULL;

	snprintf(result, 31, "%c%s", centering, g);

	return result;
}


static void mosflm_send_next(struct image *image, struct mosflm_data *mosflm)
{
	char tmp[256];
	double wavelength;

	switch ( mosflm->step )
	{
		case 1 :
		mosflm_sendline("DETECTOR ROTATION HORIZONTAL"
		                " ANTICLOCKWISE ORIGIN LL FAST HORIZONTAL"
		                " RECTANGULAR\n", mosflm);
		break;

		case 2 :
		if ( (mosflm->mp->indm & INDEXING_USE_LATTICE_TYPE)
		  && (mosflm->mp->template != NULL) )
		{
			char *symm;

			if ( cell_get_lattice_type(mosflm->mp->template)
			     == L_RHOMBOHEDRAL ) {
				mosflm_sendline("CRYSTAL R\n", mosflm);
			}

			symm = spacegroup_for_lattice(mosflm->mp->template);
			snprintf(tmp, 255, "SYMM %s\n", symm);
			//STATUS("Asking MOSFLM for '%s'\n", symm);
			free(symm);
			mosflm_sendline(tmp, mosflm);

		} else {
			mosflm_sendline("\n", mosflm);
		}
		break;

		case 3 :
		mosflm_sendline("DISTANCE 67.8\n", mosflm);
		break;

		case 4 :
		mosflm_sendline("BEAM 0.0 0.0\n", mosflm);
		break;

		case 5 :
		wavelength = image->lambda*1e10;
		snprintf(tmp, 255, "WAVELENGTH %10.5f\n", wavelength);
		mosflm_sendline(tmp, mosflm);
		break;

		case 6 :
		snprintf(tmp, 255, "NEWMAT %s\n", mosflm->newmatfile);
		mosflm_sendline(tmp, mosflm);
		break;

		case 7 :
		snprintf(tmp, 255, "IMAGE %s phi 0 0\n", mosflm->imagefile);
		mosflm_sendline(tmp, mosflm);
		break;

		case 8 :
		snprintf(tmp, 255, "AUTOINDEX DPS FILE %s"
		                   " IMAGE 1 MAXCELL 1000 REFINE\n",
		         mosflm->sptfile);

		/* "This option is only available if you e-mail Andrew Leslie
		 * and ask for it." - MOSFLM
		 * snprintf(tmp, 255, "AUTOINDEX NODISPLAY IMAGE 1 FILE %s\n",
		 *         mosflm->sptfile); */
		mosflm_sendline(tmp, mosflm);
		break;

		case 9 :
		mosflm_sendline("GO\n", mosflm);
		mosflm->finished_ok = 1;
		break;

		default:
		mosflm_sendline("exit\n", mosflm);
		return;
	}

	mosflm->step++;
}


static int mosflm_readable(struct image *image, struct mosflm_data *mosflm)
{
	int rval;
	int no_string = 0;

	rval = read(mosflm->pty, mosflm->rbuffer+mosflm->rbufpos,
	            mosflm->rbuflen-mosflm->rbufpos);
	if ( (rval == -1) || (rval == 0) ) return 1;

	mosflm->rbufpos += rval;
	assert(mosflm->rbufpos <= mosflm->rbuflen);

	while ( (!no_string) && (mosflm->rbufpos > 0) ) {

		int i;
		int block_ready = 0;
		MOSFLMInputType type = MOSFLM_INPUT_NONE;

		/* See if there's a full line in the buffer yet */
		for ( i=0; i<mosflm->rbufpos-1; i++ ) {
			/* Means the last value looked at is rbufpos-2 */

			/* Is there a prompt in the buffer? */
			if ( (i+10 <= mosflm->rbufpos)
			  && (!strncmp(mosflm->rbuffer+i, "MOSFLM => ", 10)) ) {
				block_ready = 1;
				type = MOSFLM_INPUT_PROMPT;
				break;
			}

			if ( (mosflm->rbuffer[i] == '\r')
			  && (mosflm->rbuffer[i+1] == '\n') ) {
				block_ready = 1;
				type = MOSFLM_INPUT_LINE;
				break;
			}

		}

		if ( block_ready ) {

			unsigned int new_rbuflen;
			unsigned int endbit_length;
			char *block_buffer = NULL;

			switch ( type ) {

				case MOSFLM_INPUT_LINE :
				block_buffer = malloc(i+1);
				memcpy(block_buffer, mosflm->rbuffer, i);
				block_buffer[i] = '\0';

				if ( block_buffer[0] == '\r' ) {
					memmove(block_buffer, block_buffer+1, i);
				}

				mosflm_parseline(block_buffer, image, mosflm);
				free(block_buffer);
				endbit_length = i+2;
				break;

				case MOSFLM_INPUT_PROMPT :
				mosflm_send_next(image, mosflm);
				endbit_length = i+7;
				break;

				default :

				/* Obviously, this never happens :) */
				ERROR("Unrecognised MOSFLM input mode! "
				      "I don't know how to understand MOSFLM\n");
				return 1;

			}

			/* Now the block's been parsed, it should be
			 * forgotten about */
			memmove(mosflm->rbuffer,
			        mosflm->rbuffer + endbit_length,
			        mosflm->rbuflen - endbit_length);

			/* Subtract the number of bytes removed */
			mosflm->rbufpos = mosflm->rbufpos
			                       - endbit_length;
			new_rbuflen = mosflm->rbuflen - endbit_length;
			if ( new_rbuflen == 0 ) new_rbuflen = 256;
			mosflm->rbuffer = realloc(mosflm->rbuffer,
			                               new_rbuflen);
			mosflm->rbuflen = new_rbuflen;

		} else {

			if ( mosflm->rbufpos==mosflm->rbuflen ) {

				/* More buffer space is needed */
				mosflm->rbuffer = realloc(
				                    mosflm->rbuffer,
				                    mosflm->rbuflen + 256);
				mosflm->rbuflen = mosflm->rbuflen + 256;
				/* The new space gets used at the next
				 * read, shortly... */

			}
			no_string = 1;

		}

	}

	return 0;
}


int run_mosflm(struct image *image, IndexingPrivate *ipriv)
{
	struct mosflm_data *mosflm;
	unsigned int opts;
	int status;
	int rval;

	mosflm = malloc(sizeof(struct mosflm_data));
	if ( mosflm == NULL ) {
		ERROR("Couldn't allocate memory for MOSFLM data.\n");
		return 0;
	}

	snprintf(mosflm->imagefile, 127, "xfel-%i_001.img", image->id);
	write_img(image, mosflm->imagefile); /* Dummy image */

	snprintf(mosflm->sptfile, 127, "xfel-%i_001.spt", image->id);
	write_spt(image, mosflm->sptfile);

	snprintf(mosflm->newmatfile, 127, "xfel-%i.newmat", image->id);
	remove(mosflm->newmatfile);

	mosflm->pid = forkpty(&mosflm->pty, NULL, NULL, NULL);

	if ( mosflm->pid == -1 ) {
		ERROR("Failed to fork for MOSFLM: %s\n", strerror(errno));
		free(mosflm);
		return 0;
	}
	if ( mosflm->pid == 0 ) {

		/* Child process: invoke MOSFLM */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("ipmosflm", "", (char *)NULL);
		ERROR("Failed to invoke MOSFLM.\n");
		_exit(0);

	}

	mosflm->rbuffer = malloc(256);
	mosflm->rbuflen = 256;
	mosflm->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(mosflm->pty, F_GETFL);
	fcntl(mosflm->pty, F_SETFL, opts | O_NONBLOCK);

	mosflm->step = 1;	/* This starts the "initialisation" procedure */
	mosflm->finished_ok = 0;
	mosflm->mp = (struct mosflm_private *)ipriv;
	mosflm->done = 0;
	mosflm->success = 0;

	rval = 0;
	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(mosflm->pty, &fds);

		tv.tv_sec = 30;
		tv.tv_usec = 0;

		sval = select(mosflm->pty+1, &fds, NULL, NULL, &tv);

		if ( sval == -1 ) {

			const int err = errno;

			switch ( err ) {

				case EINTR:
				STATUS("Restarting select()\n");
				break;

				default:
				ERROR("select() failed: %s\n", strerror(err));
				rval = 1;

			}

		} else if ( sval != 0 ) {
			rval = mosflm_readable(image, mosflm);
		} else {
			ERROR("No response from MOSFLM..\n");
			rval = 1;
		}

	} while ( !rval );

	close(mosflm->pty);
	free(mosflm->rbuffer);
	waitpid(mosflm->pid, &status, 0);

	if ( mosflm->finished_ok == 0 ) {
		ERROR("MOSFLM doesn't seem to be working properly.\n");
	} else {
		/* Read the mosflm NEWMAT file and get cell if found */
		read_newmat(mosflm, mosflm->newmatfile, image);
	}

	rval = mosflm->success;
	free(mosflm);
	return rval;
}


IndexingPrivate *mosflm_prepare(IndexingMethod *indm, UnitCell *cell,
                                struct detector *det, struct beam_params *beam,
                                float *ltl)
{
	struct mosflm_private *mp;
	int need_cell = 0;

	if ( *indm & INDEXING_CHECK_CELL_COMBINATIONS ) need_cell = 1;
	if ( *indm & INDEXING_CHECK_CELL_AXES ) need_cell = 1;
	if ( *indm & INDEXING_USE_LATTICE_TYPE ) need_cell = 1;

	if ( need_cell && (cell == NULL) ) {
		ERROR("Altering your MOSFLM flags because no PDB file was"
		      " provided.\n");
		*indm &= ~INDEXING_CHECK_CELL_COMBINATIONS;
		*indm &= ~INDEXING_CHECK_CELL_AXES;
		*indm &= ~INDEXING_USE_LATTICE_TYPE;
	}

	/* Flags that MOSFLM knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_CELL_COMBINATIONS
	       | INDEXING_CHECK_CELL_AXES | INDEXING_CHECK_PEAKS
	       | INDEXING_USE_LATTICE_TYPE;

	if ( *indm & INDEXING_USE_LATTICE_TYPE ) {
		if ( !((*indm & INDEXING_CHECK_CELL_COMBINATIONS)
		    || (*indm & INDEXING_CHECK_CELL_AXES)) ) {
			ERROR("WARNING: The unit cell from %s might have had "
			      "its axes permuted from the unit cell you gave.\n"
			      "If this is a problem, consider using "
			      "mosflm-axes-latt or mosflm-comb-latt instead of "
			      "mosflm-raw-latt.\n", indexer_str(*indm));
		}

	}

	mp = malloc(sizeof(struct mosflm_private));
	if ( mp == NULL ) return NULL;

	mp->ltl = ltl;
	mp->template = cell;
	mp->indm = *indm;

	return (IndexingPrivate *)mp;
}


void mosflm_cleanup(IndexingPrivate *pp)
{
	struct mosflm_private *p;
	p = (struct mosflm_private *)pp;
	free(p);
}
