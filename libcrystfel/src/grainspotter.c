/*
 * grainspotter.c
 *
 * Invoke GrainSpotter for multi-crystal autoindexing
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
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

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell.h"
#include "cell-utils.h"
#include "grainspotter.h"


#define GRAINSPOTTER_VERBOSE 0


/* Global private data, prepared once */
struct grainspotter_private
{
	IndexingMethod indm;
	UnitCell *cell;
};


/* Data needed during the run of Grainspotter */
struct grainspotter_data {

	struct grainspotter_private *gp;

	/* Low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

};


static int read_matrix(struct grainspotter_private *gp, struct image *image,
                       char *filename)
{
	FILE *fh;
	int d1;
	float d2;
	float ubi11, ubi12, ubi13;
	float ubi21, ubi22, ubi23;
	float ubi31, ubi32, ubi33;
	char line[1024];
	int r;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Can't open '%s'\n", filename);
		return 1;
	}

	/* Read and discard first line */
	if ( fgets(line, 1024, fh) == NULL ) {
		ERROR("Failed to read GFF file.\n");
		return 1;
	}

	do {

		Crystal *cr;
		UnitCell *cell;

		/* One line per grain */
		if ( fgets(line, 1024, fh) == NULL ) {
			ERROR("Failed to read GFF file.\n");
			return 1;
		}

		STATUS("'%s'\n", line);

		r = sscanf(line, "%i %f %f %f %f %f %f %f %f %f %f %f %f"
			         "%f %f %f %f %f %f %f %f %f %f %f",
			         &d1, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2,
			         &d2, &d2, &d2, &d2, &d2, &d2,
			         &ubi11, &ubi12, &ubi13,
			         &ubi21, &ubi22, &ubi23,
			         &ubi31, &ubi32, &ubi33);

		if ( r != 24 ) {
			ERROR("Only %i parameters in GFF file\n", r);
			return 1;
		}

		cell = cell_new();

		cell_set_cartesian(cell, ubi12/1e10, ubi13/1e10, ubi11/1e10,
			                 ubi22/1e10, ubi23/1e10, ubi21/1e10,
			                 ubi32/1e10, ubi33/1e10, ubi31/1e10);

		cr = crystal_new();
		if ( cr == NULL ) {
			ERROR("Failed to allocate crystal.\n");
			return 0;
		}

		crystal_set_cell(cr, cell);
		image_add_crystal(image, cr);

	} while ( !feof(fh) );

	fclose(fh);

	if ( gp->indm & INDEXING_CHECK_PEAKS ) {
		if ( !peak_sanity_check(image, image->crystals,
		                        image->n_crystals) )
		{
			free_all_crystals(image);
			return 0;
		}
	}

        return 0;
}


static void gs_parseline(char *line, struct image *image,
                         struct grainspotter_data *gs)
{
	#if GRAINSPOTTER_VERBOSE
	STATUS("%s\n", line);
	#endif
}


static int grainspotter_readable(struct image *image,
                                 struct grainspotter_data *gs)
{
	int rval;
	int no_string = 0;

	rval = read(gs->pty, gs->rbuffer+gs->rbufpos, gs->rbuflen-gs->rbufpos);

	if ( (rval == -1) || (rval == 0) ) return 1;

	gs->rbufpos += rval;
	assert(gs->rbufpos <= gs->rbuflen);

	while ( (!no_string) && (gs->rbufpos > 0) ) {

		int i;
		int block_ready = 0;

		/* See if there's a full line in the buffer yet */
		for ( i=0; i<gs->rbufpos-1; i++ ) {
			/* Means the last value looked at is rbufpos-2 */

			if ( (gs->rbuffer[i] == '\r')
			  && (gs->rbuffer[i+1] == '\n') ) {
				block_ready = 1;
				break;
			}

		}

		if ( block_ready ) {

			unsigned int new_rbuflen;
			unsigned int endbit_length;
			char *block_buffer = NULL;

			block_buffer = malloc(i+1);
			memcpy(block_buffer, gs->rbuffer, i);
			block_buffer[i] = '\0';

			if ( block_buffer[0] == '\r' ) {
				memmove(block_buffer, block_buffer+1, i);
			}

			gs_parseline(block_buffer, image, gs);
			free(block_buffer);
			endbit_length = i+2;

			/* Now the block's been parsed, it should be
			 * forgotten about */
			memmove(gs->rbuffer,
			        gs->rbuffer + endbit_length,
			        gs->rbuflen - endbit_length);

			/* Subtract the number of bytes removed */
			gs->rbufpos = gs->rbufpos  - endbit_length;
			new_rbuflen = gs->rbuflen - endbit_length;
			if ( new_rbuflen == 0 ) new_rbuflen = 256;
			gs->rbuffer = realloc(gs->rbuffer, new_rbuflen);
			gs->rbuflen = new_rbuflen;

		} else {

			if ( gs->rbufpos==gs->rbuflen ) {

				/* More buffer space is needed */
				gs->rbuffer = realloc(gs->rbuffer,
				                      gs->rbuflen + 256);
				gs->rbuflen = gs->rbuflen + 256;
				/* The new space gets used at the next
				 * read, shortly... */

			}
			no_string = 1;

		}

	}

	return 0;
}


static void write_gve(struct image *image, struct grainspotter_private *gp)
{
	FILE *fh;
	int i;
	char filename[1024];
	double a, b, c, al, be, ga;

	snprintf(filename, 1023, "xfel-%i.gve", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}

	cell_get_parameters(gp->cell, &a, &b, &c, &al, &be, &ga);
	fprintf(fh, "%.6f %.6f %.6f %.6f %.6f %.6f P\n", a*1e10, b*1e10, c*1e10,
	                     rad2deg(al), rad2deg(be), rad2deg(ga));
	fprintf(fh, "# wavelength = %.6f\n", image->lambda*1e10);
	fprintf(fh, "# wedge = 0.000000\n");
	fprintf(fh, "# ds h k l\n");
	fprintf(fh, "# xr yr zr xc yc ds eta omega\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		fprintf(fh, "%.6f %.6f %.6f 0 0 %.6f %.6f %.6f 0\n",
		        f->rz/1e10, f->rx/1e10, f->ry/1e10,
		        modulus(f->rx, f->ry, f->rz)/1e10, /* dstar */
		        rad2deg(atan2(f->ry, f->rx)), 0.0);   /* eta, omega */

	}
	fclose(fh);
}


static char *write_ini(struct image *image)
{
	FILE *fh;
	char *filename;
	double tt;

	filename = malloc(1024);
	if ( filename == NULL ) return NULL;

	snprintf(filename, 1023, "xfel-%i.ini", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		free(filename);
		return NULL;
	}

	get_q_for_panel(image->det->furthest_out_panel,
	                image->det->furthest_out_fs,
	                image->det->furthest_out_ss,
	                &tt, 1.0/image->lambda);

	fprintf(fh, "spacegroup 96\n");
	fprintf(fh, "!dsrange 0 1.3\n");
	fprintf(fh, "tthrange 0 %.2f\n", 20.0);
	fprintf(fh, "etarange 0 360\n");
	fprintf(fh, "domega 1\n");
	fprintf(fh, "omegarange -0.5 0.5\n");
	fprintf(fh, "filespecs xfel-%i.gve xfel-%i.log\n",
	            image->id, image->id);
	fprintf(fh, "cuts 3 0.1 0.5\n");
	fprintf(fh, "eulerstep 8\n");
	fprintf(fh, "uncertainties 0.1 0.2 .5\n");
	fprintf(fh, "nsigmas 2\n");
	fprintf(fh, "minfracg 0.95\n");
	fprintf(fh, "Nhkls_in_indexing 500\n");
	fprintf(fh, "random 10000\n");
	fprintf(fh, "!positionfit\n");
	fprintf(fh, "genhkl\n");
	fprintf(fh, "xfelmode\n");

	fclose(fh);

	return filename;
}


int grainspotter_index(struct image *image, IndexingPrivate *ipriv)
{
	unsigned int opts;
	int status;
	int rval;
	struct grainspotter_data *grainspotter;
	struct grainspotter_private *gp = (struct grainspotter_private *)ipriv;
	char *ini_filename;
	char gff_filename[1024];

	write_gve(image, gp);
	ini_filename = write_ini(image);

	if ( ini_filename == NULL ) {
		ERROR("Failed to write ini file for GrainSpotter.\n");
		return 0;
	}

	grainspotter = malloc(sizeof(struct grainspotter_data));
	if ( grainspotter == NULL ) {
		ERROR("Couldn't allocate memory for GrainSpotter data.\n");
		return 0;
	}

	grainspotter->gp = gp;

	snprintf(gff_filename, 1023, "xfel-%i.gff", image->id);
	remove(gff_filename);

	grainspotter->pid = forkpty(&grainspotter->pty, NULL, NULL, NULL);
	if ( grainspotter->pid == -1 ) {
		ERROR("Failed to fork for GrainSpotter: %s\n", strerror(errno));
		return 0;
	}
	if ( grainspotter->pid == 0 ) {

		/* Child process: invoke GrainSpotter */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		STATUS("Running GrainSpotter.0.93 '%s'\n", ini_filename);
		execlp("GrainSpotter.0.93", "", ini_filename, (char *)NULL);
		ERROR("Failed to invoke GrainSpotter.\n");
		_exit(0);

	}
	free(ini_filename);

	grainspotter->rbuffer = malloc(256);
	grainspotter->rbuflen = 256;
	grainspotter->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(grainspotter->pty, F_GETFL);
	fcntl(grainspotter->pty, F_SETFL, opts | O_NONBLOCK);

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(grainspotter->pty, &fds);

		tv.tv_sec = 30;
		tv.tv_usec = 0;

		sval = select(grainspotter->pty+1, &fds, NULL, NULL, &tv);

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

			}

		} else if ( sval != 0 ) {
			rval = grainspotter_readable(image, grainspotter);
		} else {
			ERROR("No response from GrainSpotter..\n");
			rval = 1;
		}

	} while ( !rval );

	close(grainspotter->pty);
	free(grainspotter->rbuffer);
	waitpid(grainspotter->pid, &status, 0);

	if ( status != 0 ) {
		ERROR("GrainSpotter doesn't seem to be working properly.\n");
		free(grainspotter);
		return 0;
	}

	if ( read_matrix(gp, image, gff_filename) != 0 ) {
		free(grainspotter);
		return 0;
	}

	/* Success! */
	free(grainspotter);
	return 1;
}



IndexingPrivate *grainspotter_prepare(IndexingMethod *indm, UnitCell *cell,
                                      struct detector *det,
                                      struct beam_params *beam, float *ltl)
{
	struct grainspotter_private *gp;

	if ( cell == NULL ) {
		ERROR("GrainSpotter needs a unit cell.\n");
		return NULL;
	}

	gp = calloc(1, sizeof(*gp));
	if ( gp == NULL ) return NULL;

	/* Flags that GrainSpotter knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_PEAKS
	       | INDEXING_USE_LATTICE_TYPE | INDEXING_USE_CELL_PARAMETERS;

	gp->cell = cell;
	gp->indm = *indm;

	return (IndexingPrivate *)gp;
}


void grainspotter_cleanup(IndexingPrivate *pp)
{
	struct grainspotter_private *p;

	p = (struct grainspotter_private *)pp;
	free(p);
}
