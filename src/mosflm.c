/*
 * mosflm.c
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * (c) 2010 Richard Kirian <rkirian@asu.edu>
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
#include "sfac.h"
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
		STATUS("No NEWMAT file (autoindexing was unsuccessful).\n");
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
	c = 1/image->lambda;

	image->candidate_cells[0] = cell_new();

	cell_set_reciprocal(image->candidate_cells[0],
                            asz*c, asy*c, asx*c,
                            bsz*c, bsy*c, bsx*c,
                            csz*c, csy*c, csx*c);

        image->ncells = 1;

        return 0;
}


/* Need to sort mosflm peaks by intensity... */
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


/* Write .spt file for mosflm */
static void write_spt(struct image *image, const char *filename)
{
	FILE *fh;
	int i;
	double fclen=67.8;  /* fake camera length in mm */
	double fpix=0.075;  /* fake pixel size in mm */
	double pix;
	double height=100;
	double sigma=1;
	int nPeaks = image_feature_count(image->features);
	struct sptline *sptlines;

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}

	fprintf(fh, "%10d %10d %10.8f %10.6f %10.6f\n", 1, 1, fpix, 1.0, 0.0);
	fprintf(fh, "%10d %10d\n", 1, 1);
	fprintf(fh, "%10.5f %10.5f\n", 0.0, 0.0);

	sptlines = malloc(sizeof(struct sptline)*nPeaks);

	for ( i=0; i<nPeaks; i++ ) {

		struct imagefeature *f;
		struct panel *p;
		double xs, ys, rx, ry;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		p = find_panel(image->det, f->fs, f->ss);
		if ( p == NULL ) continue;

		pix = 1000.0/p->res; /* pixel size in mm */
		height = f->intensity;

		xs = (f->fs-p->min_fs)*p->fsx + (f->ss-p->min_ss)*p->ssx;
		ys = (f->ss-p->min_fs)*p->fsy + (f->ss-p->min_ss)*p->ssy;
		rx = xs + p->cnx;
		ry = ys + p->cny;

		sptlines[i].x = ry*pix*fclen/p->clen/1000.0;
		sptlines[i].y = -rx*pix*fclen/p->clen/1000.0;
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

	free(sptlines);

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


static void mosflm_send_next(struct image *image, struct mosflm_data *mosflm)
{
	char tmp[256];
	char symm[32];
	const char *sg;
	double wavelength;
	double a, b, c, alpha, beta, gamma;
	int i, j;

	switch ( mosflm->step ) {

	case 1 :
		mosflm_sendline("DETECTOR ROTATION HORIZONTAL"
		                " ANTICLOCKWISE ORIGIN LL FAST HORIZONTAL"
		                " RECTANGULAR\n", mosflm);
		break;

	case 2 :
		if ( mosflm->target_cell != NULL ) {
			cell_get_parameters(mosflm->target_cell, &a, &b, &c,
				            &alpha, &beta, &gamma);
			snprintf(tmp, 255,
			         "CELL %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
		                 a*1e10, b*1e10, c*1e10,
		                 rad2deg(alpha), rad2deg(beta), rad2deg(gamma));
			mosflm_sendline(tmp, mosflm);
		} else {
			mosflm_sendline("# Do nothing\n", mosflm);
		}
		break;

	case 3 :
		if ( mosflm->target_cell != NULL ) {
			sg = cell_get_spacegroup(mosflm->target_cell);
			/* Remove white space from space group */
			j = 0;
			for ( i=0; i<strlen(sg); i++ ) {
				if (sg[i] != ' ') {
					symm[j++] = sg[i];
				}
			}
			symm[j] = '\0';
			snprintf(tmp, 255, "SYMM %s\n", symm);
			mosflm_sendline(tmp, mosflm);
		} else {
			mosflm_sendline("SYMM P1\n", mosflm);
		}
		break;

	case 4 :
		mosflm_sendline("DISTANCE 67.8\n", mosflm);
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
		snprintf(tmp, 255, "AUTOINDEX DPS FILE %s IMAGE 1\n",
		         mosflm->sptfile);
		mosflm_sendline(tmp, mosflm);
		break;

	case 10 :
		mosflm_sendline("GO\n", mosflm);
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
				mosflm->rbuflen = mosflm->rbuflen
				                                          + 256;
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
	int fail;
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
		ERROR("Failed to fork for MOSFLM\n");
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

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(mosflm->pty, &fds);

		tv.tv_sec = 10;
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

	/* Read the mosflm NEWMAT file and set cell candidate *
	 * Existence of this file means possible success. Pretty shady. */
	fail = read_newmat(mosflm->newmatfile, image);
	if ( fail ) {
		ERROR("Failed to read MOSFLM's NEWMAT file.\n");
		return;
	}

	free(mosflm);
}
