/*
 * dirax.c
 *
 * Invoke the DirAx auto-indexing program
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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
#include "dirax.h"
#include "utils.h"
#include "sfac.h"
#include "peaks.h"


#define DIRAX_VERBOSE 0

#define MAX_DIRAX_CELL_CANDIDATES (5)


typedef enum {
	DIRAX_INPUT_NONE,
	DIRAX_INPUT_LINE,
	DIRAX_INPUT_PROMPT,
	DIRAX_INPUT_ACL
} DirAxInputType;


struct dirax_data {

	/* DirAx auto-indexing low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

	/* DirAx auto-indexing high-level stuff */
	int                     step;
	int                     read_cell;
	int                     best_acl;
	int                     best_acl_nh;
	int                     acls_tried[MAX_CELL_CANDIDATES];
	int                     n_acls_tried;

};


static void dirax_parseline(const char *line, struct image *image,
                            struct dirax_data *dirax)
{
	int rf, i, di, acl, acl_nh;
	float d;

	#if DIRAX_VERBOSE
	char *copy;

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='r';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("DirAx: %s\n", copy);
	free(copy);
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
			image->candidate_cells[image->ncells] = cell_new();
			return;
		}
		i++;
	}

	/* Parse unit cell vectors as appropriate */
	if ( dirax->read_cell == 1 ) {
		/* First row of unit cell values */
		float ax, ay, az;
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &ax, &ay, &az);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line\n");
			dirax->read_cell = 0;
			free(image->candidate_cells[image->ncells]);
			return;
		}
		cell_set_cartesian_a(image->candidate_cells[image->ncells],
		                     ax*1e-10, ay*1e-10, az*1e-10);
		dirax->read_cell++;
		return;
	} else if ( dirax->read_cell == 2 ) {
		/* Second row of unit cell values */
		float bx, by, bz;
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &bx, &by, &bz);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line\n");
			dirax->read_cell = 0;
			free(image->candidate_cells[image->ncells]);
			return;
		}
		cell_set_cartesian_b(image->candidate_cells[image->ncells],
		                     bx*1e-10, by*1e-10, bz*1e-10);
		dirax->read_cell++;
		return;
	} else if ( dirax->read_cell == 3 ) {
		/* Third row of unit cell values */
		float cx, cy, cz;
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &cx, &cy, &cz);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line\n");
			dirax->read_cell = 0;
			free(image->candidate_cells[image->ncells]);
			return;
		}
		cell_set_cartesian_c(image->candidate_cells[image->ncells++],
		                     cx*1e-10, cy*1e-10, cz*1e-10);
		dirax->read_cell = 0;
		return;
	}

	dirax->read_cell = 0;

	if ( sscanf(line, "%i %i %f %f %f %f %f %f %i", &acl, &acl_nh,
	                               &d, &d, &d, &d, &d, &d, &di) == 9 ) {
		if ( acl_nh > dirax->best_acl_nh ) {

			int i, found = 0;

			for ( i=0; i<dirax->n_acls_tried; i++ ) {
				if ( dirax->acls_tried[i] == acl ) found = 1;
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

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='\0';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("To DirAx: '%s'\n", copy);
	free(copy);
	#endif

	write(dirax->pty, line, strlen(line));
}


static void dirax_send_next(struct image *image, struct dirax_data *dirax)
{
	char tmp[32];

	switch ( dirax->step ) {

	case 1 :
		dirax_sendline("\\echo off\n", dirax);
		break;

	case 2 :
		snprintf(tmp, 31, "read xfel-%i.drx\n", image->id);
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
		break;

	case 7 :
		dirax_sendline("acl\n", dirax);
		break;

	case 8 :
		if ( dirax->best_acl_nh == 0 ) {
			STATUS("No more cells to try.\n");
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

				block_buffer = malloc(i+1);
				memcpy(block_buffer, dirax->rbuffer, i);
				block_buffer[i] = '\0';

				if ( block_buffer[0] == '\r' ) {
					memmove(block_buffer, block_buffer+1, i);
				}

				dirax_parseline(block_buffer, image, dirax);
				free(block_buffer);
				endbit_length = i+2;
				break;

			case DIRAX_INPUT_PROMPT :

				dirax_send_next(image, dirax);
				endbit_length = i+7;
				break;

			case DIRAX_INPUT_ACL :

				dirax_send_next(image,dirax );
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
			dirax->rbuffer = realloc(dirax->rbuffer,
			                               new_rbuflen);
			dirax->rbuflen = new_rbuflen;

		} else {

			if ( dirax->rbufpos==dirax->rbuflen ) {

				/* More buffer space is needed */
				dirax->rbuffer = realloc(
				                    dirax->rbuffer,
				                    dirax->rbuflen + 256);
				dirax->rbuflen = dirax->rbuflen
				                                          + 256;
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

	snprintf(filename, 1023, "xfel-%i.drx", image->id);

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file xfel.drx\n");
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



void run_dirax(struct image *image)
{
	unsigned int opts;
	int status;
	int rval;
	struct dirax_data *dirax;

	write_drx(image);

	dirax = malloc(sizeof(struct dirax_data));

	dirax->pid = forkpty(&dirax->pty, NULL, NULL, NULL);
	if ( dirax->pid == -1 ) {
		ERROR("Failed to fork for DirAx\n");
		return;
	}
	if ( dirax->pid == 0 ) {

		/* Child process: invoke DirAx */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		execlp("dirax", "", (char *)NULL);
		ERROR("Failed to invoke DirAx.\n");
		_exit(0);

	}

	dirax->rbuffer = malloc(256);
	dirax->rbuflen = 256;
	dirax->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(dirax->pty, F_GETFL);
	fcntl(dirax->pty, F_SETFL, opts | O_NONBLOCK);

	dirax->step = 1;	/* This starts the "initialisation" procedure */
	dirax->read_cell = 0;
	dirax->n_acls_tried = 0;
	dirax->best_acl_nh = 0;

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(dirax->pty, &fds);

		tv.tv_sec = 10;
		tv.tv_usec = 0;

		sval = select(dirax->pty+1, &fds, NULL, NULL, &tv);

		if ( sval == -1 ) {
			int err = errno;
			ERROR("select() failed: %s\n", strerror(err));
			rval = 1;
		} else if ( sval != 0 ) {
			rval = dirax_readable(image, dirax);
		} else {
			ERROR("No response from DirAx..\n");
			rval = 1;
		}

	} while ( !rval );

	close(dirax->pty);
	free(dirax->rbuffer);
	wait(&status);

	free(dirax);

	return;
}
