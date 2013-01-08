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
#include "grainspotter.h"
#include "utils.h"
#include "peaks.h"


#define GRAINSPOTTER_VERBOSE 0


struct grainspotter_data {

	/* Low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

};


static int grainspotter_readable(struct image *image,
                                 struct grainspotter_data *grainspotter)
{
	int rval;

	rval = read(grainspotter->pty,
	            grainspotter->rbuffer+grainspotter->rbufpos,
	            grainspotter->rbuflen-grainspotter->rbufpos);

	if ( (rval == -1) || (rval == 0) ) return 1;

	/* FIXME! (if needed) */
	//grainspotter->rbufpos += rval;
	//assert(grainspotter->rbufpos <= grainspotter->rbuflen);

	return 0;
}


static void write_gve(struct image *image)
{
	FILE *fh;
	int i;
	char filename[1024];

	snprintf(filename, 1023, "xfel-%i.gve", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}
	fprintf(fh, "%f\n", 0.5);  /* Lie about the wavelength.  */

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		fprintf(fh, "%10f %10f %10f %8f\n",
		        f->rx/1e10, f->ry/1e10, f->rz/1e10, 1.0);

	}
	fclose(fh);
}


void run_grainspotter(struct image *image, UnitCell *cell)
{
	unsigned int opts;
	int status;
	int rval;
	struct grainspotter_data *grainspotter;

	write_gve(image);

	grainspotter = malloc(sizeof(struct grainspotter_data));
	if ( grainspotter == NULL ) {
		ERROR("Couldn't allocate memory for GrainSpotter data.\n");
		return;
	}

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

		execlp("GrainSpotter.0.90", "", (char *)NULL);
		ERROR("Failed to invoke GrainSpotter.\n");
		_exit(0);

	}

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

	free(grainspotter);
}
