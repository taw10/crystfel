/*
 * mosflm.c
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010      Richard Kirian <rkirian@asu.edu>
 *   2010-2012 Thomas White <taw@physics.org>
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
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>

#if HAVE_FORKPTY_LINUX
#include <pty.h>
#elif HAVE_FORKPTY_BSD
#include <util.h>
#endif


#include "image.h"
#include "mosflm.h"
#include "utils.h"
#include "peaks.h"


#define MOSFLM_VERBOSE 0


typedef enum {
	MOSFLM_INPUT_NONE,
	MOSFLM_INPUT_LINE,
	MOSFLM_INPUT_PROMPT
} MOSFLMInputType;


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
	UnitCell                *target_cell;

};


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


static int read_newmat(const char *filename, struct image *image)
{
	FILE * fh;
	float asx, asy, asz;
	float bsx, bsy, bsz;
	float csx, csy, csz;
	int n;
	double c;

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
	fclose(fh);

	/* MOSFLM "A" matrix is multiplied by lambda, so fix this */
	c = 1.0/image->lambda;

	image->candidate_cells[0] = cell_new();

	/* The relationship between the coordinates in the spot file and the
	 * resulting matrix is diabolically complicated.  This transformation
	 * seems to work, but is not derived by working through all the
	 * transformations. */
	cell_set_reciprocal(image->candidate_cells[0],
	                    -asy*c, -asz*c, asx*c,
	                    -bsy*c, -bsz*c, bsx*c,
	                    -csy*c, -csz*c, csx*c);

        image->ncells = 1;

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
static const char *spacegroup_for_lattice(UnitCell *cell)
{
	LatticeType latt;
	char centering;
	char *g = NULL;
	char *result;

	latt = cell_get_lattice_type(cell);
	centering = cell_get_centering(cell);

	switch ( latt )
	{
		case L_TRICLINIC :
		g = "1";
		break;

		case L_MONOCLINIC :
		g = "2";
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
			g = "32";
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

	switch ( mosflm->step ) {

		case 1 :
		mosflm_sendline("DETECTOR ROTATION HORIZONTAL"
		                " ANTICLOCKWISE ORIGIN LL FAST HORIZONTAL"
		                " RECTANGULAR\n", mosflm);
		break;

		case 2 :
		if ( mosflm->target_cell != NULL ) {
			const char *symm;
			symm = spacegroup_for_lattice(mosflm->target_cell);
			snprintf(tmp, 255, "SYMM %s\n", symm);
			mosflm_sendline(tmp, mosflm);
		} else {
			mosflm_sendline("SYMM P1\n", mosflm);
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


void run_mosflm(struct image *image, UnitCell *cell)
{
	struct mosflm_data *mosflm;
	unsigned int opts;
	int status;
	int rval;

	mosflm = malloc(sizeof(struct mosflm_data));
	if ( mosflm == NULL ) {
		ERROR("Couldn't allocate memory for MOSFLM data.\n");
		return;
	}

	mosflm->target_cell = cell;

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
		return;
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
			int err = errno;
			ERROR("select() failed: %s\n", strerror(err));
			rval = 1;
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
		read_newmat(mosflm->newmatfile, image);
	}

	free(mosflm);
}
