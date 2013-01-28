/*
 * xds.c
 *
 * Invoke xds crystal autoindexing
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


#define xds_VERBOSE 0


struct xds_data {

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


static int xds_readable(struct image *image,
                                 struct xds_data *xds)
{
	int rval;

	rval = read(xds->pty,
	            xds->rbuffer, xds->rbuflen);

	if ( rval == -1 ) {
		ERROR("Read failed: %s\n", strerror(errno));
		return 1;
	}

	xds->rbuffer[rval] = '\0';

	printf("xds: '%s'\n", xds->rbuffer);

	return 0;
}

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


//void run_xds(struct image *image, UnitCell *cell)
{
	unsigned int opts;
	int status;
	int rval;
	struct xds_data *xds;
	char *ini_filename;
	char gff_filename[1024];

	write_gve(image, cell);
	ini_filename = write_ini(image);

	if ( ini_filename == NULL ) {
		ERROR("Failed to write ini file for xds.\n");
		return;
	}

	xds = malloc(sizeof(struct xds_data));
	if ( xds == NULL ) {
		ERROR("Couldn't allocate memory for xds data.\n");
		return;
	}

	snprintf(gff_filename, 1023, "xfel-%i.gff", image->id);
	remove(gff_filename);

	xds->pid = forkpty(&xds->pty, NULL, NULL, NULL);
	if ( xds->pid == -1 ) {
		ERROR("Failed to fork for xds: %s\n", strerror(errno));
		return;
	}
	if ( xds->pid == 0 ) {

		/* Child process: invoke xds */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		STATUS("Running xds.0.90 '%s'\n", ini_filename);
		execlp("xds.0.90", "", ini_filename, (char *)NULL);
		ERROR("Failed to invoke xds.\n");
		_exit(0);

	}
	free(ini_filename);

	xds->rbuffer = malloc(256);
	xds->rbuflen = 256;
	xds->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(xds->pty, F_GETFL);
	fcntl(xds->pty, F_SETFL, opts | O_NONBLOCK);

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(xds->pty, &fds);

		tv.tv_sec = 20000;
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

			}

		} else if ( sval != 0 ) {
			rval = xds_readable(image, xds);
		} else {
			ERROR("No response from xds..\n");
			rval = 1;
		}

	} while ( !rval );

	close(xds->pty);
	free(xds->rbuffer);
	waitpid(xds->pid, &status, 0);

	if ( status != 0 ) {
		ERROR("xds doesn't seem to be working properly.\n");
	}

	if ( read_matrix(image, gff_filename) != 0 ) {
		ERROR("Failed to read matrix\n");
	}

	free(xds);
}
