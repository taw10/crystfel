/*
 * mosflm.c
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010      Richard Kirian <rkirian@asu.edu>
 *   2010-2018 Thomas White <taw@physics.org>
 *   2014      Takanori Nakane <nakane.t@gmail.com>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>

#ifdef HAVE_FORKPTY_PTY_H
#include <pty.h>
#endif
#ifdef HAVE_FORKPTY_UTIL_H
#include <util.h>
#endif

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

/** \file mosflm.h */

#define MOSFLM_VERBOSE 0
#define FAKE_CLEN (0.1)


typedef enum {
	MOSFLM_INPUT_NONE,
	MOSFLM_INPUT_LINE,
	MOSFLM_INPUT_PROMPT
} MOSFLMInputType;



struct mosflm_private {
	IndexingMethod          indm;
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

static int check_mosflm_cell(struct mosflm_private *mp, struct image *image,
                             UnitCell *cell)
{
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

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 0;
	}

	crystal_set_cell(cr, cell);

	image_add_crystal(image, cr);

	return 1;
}


static void mosflm_parseline(const char *line, struct image *image,
                             struct mosflm_data *dirax)
{
	if ( MOSFLM_VERBOSE || (strncmp(line, "Invocation:", 11) == 0) ) {
		char *copy;
		int i;

		copy = cfstrdup(line);
		for ( i=0; i<strlen(copy); i++ ) {
			if ( copy[i] == '\r' ) copy[i]='r';
			if ( copy[i] == '\n' ) copy[i]='\0';
		}
		STATUS("MOSFLM: %s\n", copy);
		cffree(copy);
	}
}


/* This is the opposite of mosflm_spacegroup_for_lattice() below.
 * Note that this is not general, just a set of rules for interpreting MOSFLM's
 * output. */
static LatticeType mosflm_spacegroup_to_lattice(const char *sg,
                                                char *ua, char *cen)
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
		if ( (sg[0] == 'H') || (sg[0] == 'P') ) {
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
		fclose(fh);
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
	latt = mosflm_spacegroup_to_lattice(symm+5, &ua, &cen);

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

	if ( check_mosflm_cell(mosflm->mp, image, cell) ) {
		mosflm->success = 1;
		mosflm->done = 1;
	}

        return 0;
}


/* Write .spt file for mosflm */
static void write_spt(struct image *image, const char *filename)
{
	FILE *fh;
	int i;
	int n;

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}

	/* Number of pixels in x, number of pixels in y, pixel size (mm),
	 * YSCALE, OMEGA */
	fputs("1 1 0.0 1.0 0.0\n", fh);

	/* INVERTX, ISWUNG */
	fputs("0 1\n", fh);

	/* XBEAM, YBEAM */
	fputs("0.0 0.0\n", fh);

	n = image_feature_count(image->features);
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		double ttx, tty, x, y;
		double r[3];

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		detgeom_transform_coords(&image->detgeom->panels[f->pn],
		                         f->fs, f->ss, image->lambda,
		                         0.0, 0.0, r);

		ttx = angle_between_2d(0.0, 1.0,
		                       r[0], 1.0/image->lambda + r[2]);
		tty = angle_between_2d(0.0, 1.0,
		                       r[1], 1.0/image->lambda + r[2]);
		if ( r[0] < 0.0 ) ttx *= -1.0;
		if ( r[1] < 0.0 ) tty *= -1.0;
		x = -tan(ttx)*FAKE_CLEN;
		y = tan(tty)*FAKE_CLEN;

		fprintf(fh, "%10.2f %10.2f 0.0 0.0 1000.0 10.0\n",
		        x*1e3, y*1e3);

	}

	fputs("-999.0 -999.0 -999.0 -999.0 -999.0 -999.0\n", fh);

	fclose(fh);
}


/* Write a dummy 1x1 pixel image file for MOSFLM.  Without post refinement,
 * MOSFLM will ignore this, but it must be present. */
static void write_img(struct image *image, const char *filename)
{
	FILE *fh;
	unsigned short int *intimage;

	intimage = cfmalloc(sizeof(unsigned short int));
	intimage[0] = 1;

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		cffree(intimage);
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
	cffree(intimage);
	fclose(fh);
}


static void mosflm_sendline(const char *line, struct mosflm_data *mosflm)
{
	#if MOSFLM_VERBOSE
	char *copy;
	int i;

	copy = cfstrdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='\0';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("To MOSFLM: '%s'\n", copy);
	cffree(copy);
	#endif

	if ( write(mosflm->pty, line, strlen(line)) == -1 ) {
		ERROR("write() to MOSFLM failed: %s\n", strerror(errno));
	}
}


/* Turn what we know about the unit cell into something which we can give to
 * MOSFLM to make it give us only indexing results compatible with the cell. */
static char *mosflm_spacegroup_for_lattice(UnitCell *cell)
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

	result = cfmalloc(32);
	if ( result == NULL ) return NULL;

	snprintf(result, 31, "%c%s", centering, g);

	return result;
}


static void mosflm_send_next(struct image *image, struct mosflm_data *mosflm)
{
	char tmp[256];
	double wavelength;
	double a = 0, b = 0, c = 0, alpha = 0, beta = 0, gamma = 0;

	switch ( mosflm->step )
	{
		case 1 :
		/* Backwards-compatible workaround for different Mosflm behaviour
		 * in version 7.2.2 */
		mosflm_sendline("DETECTOR ADSC\n", mosflm);
		break;

		case 2 :
		mosflm_sendline("DETECTOR ROTATION HORIZONTAL"
		                " ANTICLOCKWISE ORIGIN LL FAST HORIZONTAL"
		                " RECTANGULAR\n", mosflm);
		break;

		case 3 :
		if ( (mosflm->mp->indm & INDEXING_USE_LATTICE_TYPE)
		  && (mosflm->mp->template != NULL) )
		{
			char *symm;

			if ( cell_get_lattice_type(mosflm->mp->template)
			     == L_RHOMBOHEDRAL ) {
				mosflm_sendline("CRYSTAL R\n", mosflm);
			}

			symm = mosflm_spacegroup_for_lattice(mosflm->mp->template);
			snprintf(tmp, 255, "SYMM %s\n", symm);
			//STATUS("Asking MOSFLM for '%s'\n", symm);
			cffree(symm);
			mosflm_sendline(tmp, mosflm);

		} else {
			mosflm_sendline("\n", mosflm);
		}
		break;

		case 4 :
		snprintf(tmp, 255, "DISTANCE %f\n", FAKE_CLEN*1e3);
		mosflm_sendline(tmp, mosflm);
		break;

		case 5 :
		mosflm_sendline("BEAM 0.0 0.0\n", mosflm);
		break;

		case 6 :
		wavelength = image->lambda*1e10;
		snprintf(tmp, 255, "WAVELENGTH %10.5f\n", wavelength);
		mosflm_sendline(tmp, mosflm);
		break;

		case 7 :
		snprintf(tmp, 255, "NEWMAT %s\n", mosflm->newmatfile);
		mosflm_sendline(tmp, mosflm);
		break;

		case 8 :
		snprintf(tmp, 255, "IMAGE %s phi 0 0\n", mosflm->imagefile);
		mosflm_sendline(tmp, mosflm);
		break;

		case 9 :
		if ( mosflm->mp->indm & INDEXING_USE_CELL_PARAMETERS ) {
			char cen;
			cell_get_parameters(mosflm->mp->template,
			                    &a, &b, &c, &alpha, &beta, &gamma);
			cen = cell_get_centering(mosflm->mp->template);
			/* Old MOSFLM simply ignores CELL and CENTERING subkeywords.
			 * So this doesn't do any harm.  */
			snprintf(tmp, 255, "AUTOINDEX DPS FILE %s IMAGE 1 "
			         "MAXCELL 1000 REFINE "
			         "CELL %.1f %.1f %.1f %.1f %.1f %.1f "
			         "CENTERING %c\n",
			         mosflm->sptfile, a*1e10, b*1e10, c*1e10,
			         rad2deg(alpha), rad2deg(beta), rad2deg(gamma),
			         cen);
                } else {
			snprintf(tmp, 255, "AUTOINDEX DPS FILE %s IMAGE 1 "
			         "MAXCELL 1000 REFINE\n", mosflm->sptfile);
		}
		mosflm_sendline(tmp, mosflm);
		break;

		case 10 :
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
				block_buffer = cfmalloc(i+1);
				memcpy(block_buffer, mosflm->rbuffer, i);
				block_buffer[i] = '\0';

				if ( block_buffer[0] == '\r' ) {
					memmove(block_buffer, block_buffer+1, i);
				}

				mosflm_parseline(block_buffer, image, mosflm);
				cffree(block_buffer);
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
			mosflm->rbuffer = cfrealloc(mosflm->rbuffer,
			                            new_rbuflen);
			mosflm->rbuflen = new_rbuflen;

		} else {

			if ( mosflm->rbufpos==mosflm->rbuflen ) {

				/* More buffer space is needed */
				mosflm->rbuffer = cfrealloc(mosflm->rbuffer,
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


int run_mosflm(struct image *image, void *ipriv)
{
	struct mosflm_data *mosflm;
	unsigned int opts;
	int status;
	int rval;

	mosflm = cfmalloc(sizeof(struct mosflm_data));
	if ( mosflm == NULL ) {
		ERROR("Couldn't allocate memory for MOSFLM data.\n");
		return 0;
	}

	snprintf(mosflm->imagefile, 127, "xfel_001.img");
	write_img(image, mosflm->imagefile); /* Dummy image */

	snprintf(mosflm->sptfile, 127, "xfel_001.spt");
	write_spt(image, mosflm->sptfile);

	snprintf(mosflm->newmatfile, 127, "xfel.newmat");
	remove(mosflm->newmatfile);

	mosflm->pid = forkpty(&mosflm->pty, NULL, NULL, NULL);

	if ( mosflm->pid == -1 ) {
		ERROR("Failed to fork for MOSFLM: %s\n", strerror(errno));
		cffree(mosflm);
		return 0;
	}
	if ( mosflm->pid == 0 ) {

		/* Child process: invoke MOSFLM */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("mosflm", "mosflm", "-n", (char *)NULL);
		execlp("ipmosflm", "ipmosflm", "-n", (char *)NULL);
		ERROR("Invocation: Failed to invoke MOSFLM: %s\n",
		      strerror(errno));
		_exit(0);

	}

	mosflm->rbuffer = cfmalloc(256);
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
	cffree(mosflm->rbuffer);
	waitpid(mosflm->pid, &status, 0);

	if ( mosflm->finished_ok == 0 ) {
		ERROR("MOSFLM doesn't seem to be working properly.\n");
	} else {
		/* Read the mosflm NEWMAT file and get cell if found */
		read_newmat(mosflm, mosflm->newmatfile, image);
	}

	rval = mosflm->success;
	cffree(mosflm);
	return rval;
}


void *mosflm_prepare(IndexingMethod *indm, UnitCell *cell)
{
	struct mosflm_private *mp;

	if ( mosflm_probe(cell) == NULL ) {
		ERROR("Mosflm does not appear to run properly.\n");
		ERROR("Please check your Mosflm installation.\n");
		return NULL;
	}

	/* Flags that MOSFLM knows about */
	*indm &= INDEXING_METHOD_MASK
	       | INDEXING_USE_LATTICE_TYPE | INDEXING_USE_CELL_PARAMETERS;

	if ( (cell != NULL) && (cell_get_centering(cell) == 'I')
	  && (cell_get_lattice_type(cell) == L_MONOCLINIC) )
	{
		ERROR("WARNING: Mosflm always gives the monoclinic C cell in "
		      "preference to the monoclinic I cell choice.\n");
		ERROR("To get a higher indexing rate, convert your cell to the "
		      "monoclinic C cell choice.\n");
	}

	mp = cfmalloc(sizeof(struct mosflm_private));
	if ( mp == NULL ) return NULL;

	mp->template = cell;
	mp->indm = *indm;

	return (IndexingPrivate *)mp;
}


void mosflm_cleanup(void *pp)
{
	struct mosflm_private *p;
	p = (struct mosflm_private *)pp;
	cffree(p);
}


const char *mosflm_probe(UnitCell *cell)
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

		execlp("mosflm", "mosflm", "-n", (char *)NULL);
		execlp("ipmosflm", "ipmosflm", "-n", (char *)NULL);
		_exit(1);

	}

	fh = fdopen(pty, "r");

	for ( l=0; l<10; l++ ) {
		if ( fgets(line, 1024, fh) != NULL ) {
			char *pos = strstr(line, "Mosflm version ");
			if ( pos != NULL ) ok = 1;
		}
	}

	fclose(fh);
	close(pty);
	waitpid(pid, &status, 0);

	if ( !ok ) return NULL;

	if ( cell_has_parameters(cell) ) return "mosflm-cell-nolatt,mosflm-latt-nocell";
	if ( cell != NULL ) return "mosflm-latt-nocell";
	return "mosflm-nolatt-nocell";
}
