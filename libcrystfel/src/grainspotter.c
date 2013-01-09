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


#define GRAINSPOTTER_VERBOSE 0


struct grainspotter_data {

	/* Low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

};


static int read_matrix(struct image *image, char *filename)
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

	/* One line per grain */
	if ( fgets(line, 1024, fh) == NULL ) {
		ERROR("Failed to read GFF file.\n");
		return 1;
	}

	STATUS("'%s'\n", line);

	r = sscanf(line,
         "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
         &d1, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2,
         &ubi11, &ubi12, &ubi13, &ubi21, &ubi22, &ubi23, &ubi31, &ubi32, &ubi33);

	if ( r != 24 ) {
		ERROR("Only %i parameters in GFF file\n", r);
		return 1;
	}

	fclose(fh);

	image->candidate_cells[0] = cell_new();

//	cell_set_cartesian(image->candidate_cells[0],
//	                   ubi11/1e10, ubi12/1e10, ubi13/1e10,
//	                   ubi21/1e10, ubi22/1e10, ubi23/1e10,
//	                   ubi31/1e10, ubi32/1e10, ubi33/1e10);

	cell_set_cartesian(image->candidate_cells[0],
	                   ubi12/1e10, ubi13/1e10, ubi11/1e10,
	                   ubi22/1e10, ubi23/1e10, ubi21/1e10,
	                   ubi32/1e10, ubi33/1e10, ubi31/1e10);


        image->ncells = 1;
        cell_print(image->candidate_cells[0]);

        return 0;
}


static int grainspotter_readable(struct image *image,
                                 struct grainspotter_data *grainspotter)
{
	int rval;

	rval = read(grainspotter->pty,
	            grainspotter->rbuffer, grainspotter->rbuflen);

	if ( rval == -1 ) {
		ERROR("Read failed: %s\n", strerror(errno));
		return 1;
	}

	grainspotter->rbuffer[rval] = '\0';

	printf("GrainSpotter: '%s'\n", grainspotter->rbuffer);

	return 0;
}


static void write_gve(struct image *image, UnitCell *cell)
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

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
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
		        atan2(f->ry, f->rx), 0.0);   /* eta, omega */

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

	fclose(fh);

	return filename;
}


void run_grainspotter(struct image *image, UnitCell *cell)
{
	unsigned int opts;
	int status;
	int rval;
	struct grainspotter_data *grainspotter;
	char *ini_filename;
	char gff_filename[1024];

	write_gve(image, cell);
	ini_filename = write_ini(image);

	if ( ini_filename == NULL ) {
		ERROR("Failed to write ini file for GrainSpotter.\n");
		return;
	}

	grainspotter = malloc(sizeof(struct grainspotter_data));
	if ( grainspotter == NULL ) {
		ERROR("Couldn't allocate memory for GrainSpotter data.\n");
		return;
	}

	snprintf(gff_filename, 1023, "xfel-%i.gff", image->id);
	remove(gff_filename);

	grainspotter->pid = forkpty(&grainspotter->pty, NULL, NULL, NULL);
	if ( grainspotter->pid == -1 ) {
		ERROR("Failed to fork for GrainSpotter: %s\n", strerror(errno));
		return;
	}
	if ( grainspotter->pid == 0 ) {

		/* Child process: invoke GrainSpotter */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		STATUS("Running GrainSpotter.0.90 '%s'\n", ini_filename);
		execlp("GrainSpotter.0.90", "", ini_filename, (char *)NULL);
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
	}

	if ( read_matrix(image, gff_filename) != 0 ) {
		ERROR("Failed to read matrix\n");
	}

	free(grainspotter);
}
