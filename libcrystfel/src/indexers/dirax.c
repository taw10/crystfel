/*
 * dirax.c
 *
 * Invoke the DirAx auto-indexing program
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2021 Thomas White <taw@physics.org>
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

#include "image.h"
#include "dirax.h"
#include "utils.h"
#include "peaks.h"
#include "cell-utils.h"


/** \file dirax.h */

#define DIRAX_VERBOSE 0

#define MAX_DIRAX_CELL_CANDIDATES (5)


typedef enum {
	DIRAX_INPUT_NONE,
	DIRAX_INPUT_LINE,
	DIRAX_INPUT_PROMPT,
	DIRAX_INPUT_ACL
} DirAxInputType;


struct dirax_private {
	IndexingMethod          indm;
	UnitCell                *template;
};


struct dirax_data {

	/* DirAx auto-indexing low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

	/* DirAx auto-indexing high-level stuff */
	int                     step;
	int                     finished_ok;
	int                     read_cell;
	int                     best_acl;
	int                     best_acl_nh;
	int                     acls_tried[MAX_DIRAX_CELL_CANDIDATES];
	int                     n_acls_tried;
	int                     done;
	int                     success;

	float                   ax;
	float                   ay;
	float                   az;
	float                   bx;
	float                   by;
	float                   bz;
	float                   cx;
	float                   cy;
	float                   cz;

	struct dirax_private    *dp;

};


static void dirax_parseline(const char *line, struct image *image,
                            struct dirax_data *dirax)
{
	int rf, i, di, acl, acl_nh;
	float d;

	#if DIRAX_VERBOSE
	char *copy;

	copy = cfstrdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='r';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("DirAx: %s\n", copy);
	cffree(copy);
	#endif

	if ( strstr(line, "reflections from file") ) {
		ERROR("DirAx can't understand this data.\n");
		return;
	}

	/* Is this the first line of a unit cell specification? */
	rf = 0; i = 0;
	while ( (i<strlen(line)) && ((line[i] == 'R')
				|| (line[i] == 'D') || (line[i] == ' ')) ) {
		if ( line[i] == 'R' ) rf = 1;
		if ( (line[i] == 'D') && rf ) {
			dirax->read_cell = 1;
			return;
		}
		i++;
	}

	/* Parse unit cell vectors as appropriate */
	if ( dirax->read_cell == 1 ) {

		/* First row of unit cell values */
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &dirax->ax, &dirax->ay, &dirax->az);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line:\n");
			ERROR("'%s'\n", line);
			dirax->read_cell = 0;
			return;
		}
		dirax->ax *= 1e-10;
		dirax->ay *= 1e-10;
		dirax->az *= 1e-10;
		dirax->read_cell++;
		return;

	} else if ( dirax->read_cell == 2 ) {

		/* Second row of unit cell values */
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &dirax->bx, &dirax->by, &dirax->bz);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line:\n");
			ERROR("'%s'\n", line);
			dirax->read_cell = 0;
			return;
		}
		dirax->bx *= 1e-10;
		dirax->by *= 1e-10;
		dirax->bz *= 1e-10;
		dirax->read_cell++;
		return;

	} else if ( dirax->read_cell == 3 ) {

		/* Third row of unit cell values */
		int r;
		UnitCell *cell;
		Crystal *cr;

		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &dirax->cx, &dirax->cy, &dirax->cz);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line:\n");
			ERROR("'%s'\n", line);
			dirax->read_cell = 0;
			return;
		}
		dirax->cx *= 1e-10;
		dirax->cy *= 1e-10;
		dirax->cz *= 1e-10;
		dirax->read_cell = 0;

		cell = cell_new();
		cell_set_cartesian(cell, dirax->ax, dirax->ay, dirax->az,
		                         dirax->bx, dirax->by, dirax->bz,
		                         dirax->cx, dirax->cy, dirax->cz);

		/* Finished reading a cell. */

		cr = crystal_new();
		if ( cr == NULL ) {
			ERROR("Failed to allocate crystal.\n");
			return;
		}
		crystal_set_cell(cr, cell);
		image_add_crystal(image, cr);
		dirax->done = 1;
		dirax->success = 1;

		return;

	}

	dirax->read_cell = 0;

	if ( sscanf(line, "%i %i %f %f %f %f %f %f %i", &acl, &acl_nh,
	                               &d, &d, &d, &d, &d, &d, &di) == 9 ) {
		if ( acl_nh > dirax->best_acl_nh ) {

			int i_acl, found = 0;

			for ( i_acl=0; i_acl<dirax->n_acls_tried; i_acl++ ) {
				if ( dirax->acls_tried[i_acl] == acl ) found = 1;
			}

			if ( !found ) {
				dirax->best_acl_nh = acl_nh;
				dirax->best_acl = acl;
			}

		}
	}
}


static void dirax_sendline(const char *line, struct dirax_data *dirax)
{
	#if DIRAX_VERBOSE
	char *copy;
	int i;

	copy = cfstrdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='\0';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("To DirAx: '%s'\n", copy);
	cffree(copy);
	#endif

	if ( write(dirax->pty, line, strlen(line)) == -1 ) {
		ERROR("write() To dirax failed: %s\n", strerror(errno));
	}
}


static void dirax_send_next(struct image *image, struct dirax_data *dirax)
{
	char tmp[32];

	switch ( dirax->step ) {

		case 1 :
		dirax_sendline("\\echo off\n", dirax);
		break;

		case 2 :
		snprintf(tmp, 31, "read xfel.drx\n");
		dirax_sendline(tmp, dirax);
		break;

		case 3 :
		dirax_sendline("dmax 1000\n", dirax);
		break;

		case 4 :
		dirax_sendline("indexfit 2\n", dirax);
		break;

		case 5 :
		dirax_sendline("levelfit 1000\n", dirax);
		break;

		case 6 :
		dirax_sendline("go\n", dirax);
		dirax->finished_ok = 1;
		break;

		case 7 :
		dirax_sendline("acl\n", dirax);
		break;

		case 8 :
		if ( dirax->best_acl_nh == 0 ) {
			/* At this point, DirAx is presenting its ACL prompt
			 * and waiting for a single number.  Use an extra
			 * newline to choose automatic ACL selection before
			 * exiting. */
			dirax_sendline("\nexit\n", dirax);
			break;
		}
		snprintf(tmp, 31, "%i\n", dirax->best_acl);
		dirax->acls_tried[dirax->n_acls_tried++] = dirax->best_acl;
		dirax_sendline(tmp, dirax);
		break;

		case 9 :
		dirax_sendline("cell\n", dirax);
		break;

		case 10 :
		if ( dirax->n_acls_tried == MAX_DIRAX_CELL_CANDIDATES ) {
			dirax_sendline("exit\n", dirax);
		} else {
			/* Go back round for another cell */
			dirax->best_acl_nh = 0;
			dirax->step = 7;
			dirax_sendline("acl\n", dirax);
		}
		break;

		default:
		dirax_sendline("exit\n", dirax);
		return;

	}

	dirax->step++;
}


static int dirax_readable(struct image *image, struct dirax_data *dirax)
{
	int rval;
	int no_string = 0;

	rval = read(dirax->pty, dirax->rbuffer+dirax->rbufpos,
	            dirax->rbuflen-dirax->rbufpos);

	if ( (rval == -1) || (rval == 0) ) return 1;

	dirax->rbufpos += rval;
	assert(dirax->rbufpos <= dirax->rbuflen);

	while ( (!no_string) && (dirax->rbufpos > 0) ) {

		int i;
		int block_ready = 0;
		DirAxInputType type = DIRAX_INPUT_NONE;

		/* See if there's a full line in the buffer yet */
		for ( i=0; i<dirax->rbufpos-1; i++ ) {
			/* Means the last value looked at is rbufpos-2 */

			/* Is there a prompt in the buffer? */
			if ( (i+7 <= dirax->rbufpos)
			  && (!strncmp(dirax->rbuffer+i, "Dirax> ", 7)) ) {
				block_ready = 1;
				type = DIRAX_INPUT_PROMPT;
				break;
			}

			/* How about an ACL prompt? */
			if ( (i+10 <= dirax->rbufpos)
			  && (!strncmp(dirax->rbuffer+i, "acl/auto [", 10)) ) {
				block_ready = 1;
				type = DIRAX_INPUT_ACL;
				break;
			}

			if ( (dirax->rbuffer[i] == '\r')
			  && (dirax->rbuffer[i+1] == '\n') ) {
				block_ready = 1;
				type = DIRAX_INPUT_LINE;
				break;
			}

		}

		if ( block_ready ) {

			unsigned int new_rbuflen;
			unsigned int endbit_length;
			char *block_buffer = NULL;

			switch ( type ) {

				case DIRAX_INPUT_LINE :
				/* Make buffer a bit too big to keep Valgrind
				 * quiet about alignment errors */
				block_buffer = cfmalloc(i+4);
				memcpy(block_buffer, dirax->rbuffer, i);
				block_buffer[i] = '\0';

				if ( block_buffer[0] == '\r' ) {
					memmove(block_buffer, block_buffer+1, i);
				}

				dirax_parseline(block_buffer, image, dirax);
				cffree(block_buffer);
				endbit_length = i+2;
				break;

				case DIRAX_INPUT_PROMPT :
				dirax_send_next(image, dirax);
				endbit_length = i+7;
				break;

				case DIRAX_INPUT_ACL :
				dirax_send_next(image,dirax);
				endbit_length = i+10;
				break;

				default :
				/* Obviously, this never happens :) */
				ERROR("Unrecognised DirAx input mode! "
				      "I don't know how to understand DirAx\n");
				return 1;

			}

			/* Now the block's been parsed, it should be
			 * forgotten about */
			memmove(dirax->rbuffer,
			        dirax->rbuffer + endbit_length,
			        dirax->rbuflen - endbit_length);

			/* Subtract the number of bytes removed */
			dirax->rbufpos = dirax->rbufpos
			                       - endbit_length;
			new_rbuflen = dirax->rbuflen - endbit_length;
			if ( new_rbuflen == 0 ) new_rbuflen = 256;
			dirax->rbuffer = cfrealloc(dirax->rbuffer, new_rbuflen);
			dirax->rbuflen = new_rbuflen;

		} else {

			if ( dirax->rbufpos == dirax->rbuflen ) {

				/* More buffer space is needed */
				dirax->rbuffer = cfrealloc(dirax->rbuffer,
				                           dirax->rbuflen + 256);
				dirax->rbuflen = dirax->rbuflen + 256;
				/* The new space gets used at the next
				 * read, shortly... */

			}
			no_string = 1;

		}

	}

	return 0;
}


static void write_drx(struct image *image)
{
	FILE *fh;
	int i;
	char filename[1024];

	snprintf(filename, 1023, "xfel.drx");

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		return;
	}
	fprintf(fh, "%f\n", 0.5);  /* Lie about the wavelength.  */

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double r[3];

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		detgeom_transform_coords(&image->detgeom->panels[f->pn],
		                         f->fs, f->ss, image->lambda,
		                         0.0, 0.0, r);

		fprintf(fh, "%10f %10f %10f %8f\n",
		        r[0]/1e10, r[1]/1e10, r[2]/1e10, 1.0);

	}
	fclose(fh);
}


int run_dirax(struct image *image, void *ipriv)
{
	unsigned int opts;
	int status;
	int rval;
	struct dirax_data *dirax;

	write_drx(image);

	dirax = cfmalloc(sizeof(struct dirax_data));
	if ( dirax == NULL ) {
		ERROR("Couldn't allocate memory for DirAx data.\n");
		return 0;
	}

	dirax->pid = forkpty(&dirax->pty, NULL, NULL, NULL);
	if ( dirax->pid == -1 ) {
		ERROR("Failed to fork for DirAx: %s\n", strerror(errno));
		return 0;
	}
	if ( dirax->pid == 0 ) {

		/* Child process: invoke DirAx */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("dirax", "dirax", (char *)NULL);
		ERROR("Failed to invoke DirAx.\n");
		_exit(0);

	}

	dirax->rbuffer = cfmalloc(256);
	dirax->rbuflen = 256;
	dirax->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(dirax->pty, F_GETFL);
	fcntl(dirax->pty, F_SETFL, opts | O_NONBLOCK);

	dirax->step = 1;	/* This starts the "initialisation" procedure */
	dirax->finished_ok = 0;
	dirax->read_cell = 0;
	dirax->n_acls_tried = 0;
	dirax->best_acl_nh = 0;
	dirax->done = 0;
	dirax->success = 0;
	dirax->dp = (struct dirax_private *)ipriv;

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(dirax->pty, &fds);

		tv.tv_sec = 30;
		tv.tv_usec = 0;

		sval = select(dirax->pty+1, &fds, NULL, NULL, &tv);

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
			rval = dirax_readable(image, dirax);
		} else {
			ERROR("No response from DirAx..\n");
			rval = 1;
		}

	} while ( !rval && !dirax->success );

	close(dirax->pty);
	cffree(dirax->rbuffer);
	waitpid(dirax->pid, &status, 0);

	if ( dirax->finished_ok == 0 ) {
		ERROR("DirAx doesn't seem to be working properly.\n");
	}

	rval = dirax->success;
	cffree(dirax);
	return rval;
}


void *dirax_prepare(IndexingMethod *indm, UnitCell *cell)
{
	struct dirax_private *dp;

	if ( dirax_probe(cell) == NULL ) {
		ERROR("DirAx does not appear to run properly.\n");
		ERROR("Please check your DirAx installation.\n");
		return NULL;
	}

	/* Flags that DirAx knows about */
	*indm &= INDEXING_METHOD_MASK;

	dp = cfmalloc(sizeof(struct dirax_private));
	if ( dp == NULL ) return NULL;

	dp->template = cell;
	dp->indm = *indm;

	return (IndexingPrivate *)dp;
}


void dirax_cleanup(void *pp)
{
	struct dirax_private *p;
	p = (struct dirax_private *)pp;
	cffree(p);
}


const char *dirax_probe(UnitCell *cell)
{
	pid_t pid;
	int pty;
	int status;
	FILE *fh;
	char line[1024];
	int ok = 0;

	pid = forkpty(&pty, NULL, NULL, NULL);
	if ( pid == -1 ) {
		return NULL;
	}
	if ( pid == 0 ) {

		/* Child process: invoke DirAx */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("dirax", "dirax", (char *)NULL);
		_exit(1);

	}

	fh = fdopen(pty, "r");
	if ( fgets(line, 1024, fh) != NULL ) {
		if ( strncmp(line, "dirax", 5) == 0 ) {
			ok = 1;
		}
	}

	fclose(fh);
	close(pty);
	waitpid(pid, &status, 0);

	if ( ok ) return "dirax";
	return NULL;
}
