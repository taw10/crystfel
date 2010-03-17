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
	DIRAX_INPUT_PROMPT
} DirAxInputType;


static void dirax_parseline(const char *line, struct image *image)
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
			image->dirax_read_cell = 1;
			image->candidate_cells[image->ncells] = cell_new();
			return;
		}
		i++;
	}

	/* Parse unit cell vectors as appropriate */
	if ( image->dirax_read_cell == 1 ) {
		/* First row of unit cell values */
		float ax, ay, az;
		sscanf(line, "%f %f %f %f %f %f", &d, &d, &d, &ax, &ay, &az);
		cell_set_cartesian_a(image->candidate_cells[image->ncells],
		                     ax*1e-10, ay*1e-10, az*1e-10);
		image->dirax_read_cell++;
		return;
	} else if ( image->dirax_read_cell == 2 ) {
		/* First row of unit cell values */
		float bx, by, bz;
		sscanf(line, "%f %f %f %f %f %f", &d, &d, &d, &bx, &by, &bz);
		cell_set_cartesian_b(image->candidate_cells[image->ncells],
		                     bx*1e-10, by*1e-10, bz*1e-10);
		image->dirax_read_cell++;
		return;
	} else if ( image->dirax_read_cell == 3 ) {
		/* First row of unit cell values */
		float cx, cy, cz;
		sscanf(line, "%f %f %f %f %f %f", &d, &d, &d, &cx, &cy, &cz);
		cell_set_cartesian_c(image->candidate_cells[image->ncells++],
		                     cx*1e-10, cy*1e-10, cz*1e-10);
		image->dirax_read_cell = 0;
		return;
	}

	image->dirax_read_cell = 0;

	if ( sscanf(line, "%i %i %f %f %f %f %f %f %i", &acl, &acl_nh,
	                               &d, &d, &d, &d, &d, &d, &di) == 9 ) {
		if ( acl_nh > image->best_acl_nh ) {

			int i, found = 0;

			for ( i=0; i<image->n_acls_tried; i++ ) {
				if ( image->acls_tried[i] == acl ) found = 1;
			}

			if ( !found ) {
				image->best_acl_nh = acl_nh;
				image->best_acl = acl;
			}

		}
	}
}


static void dirax_sendline(const char *line, struct image *image)
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

	write(image->dirax_pty, line, strlen(line));
}


static void dirax_send_next(struct image *image)
{
	char tmp[32];

	switch ( image->dirax_step ) {

	case 1 :
		dirax_sendline("\\echo off\n", image);
		break;

	case 2 :
		dirax_sendline("read xfel.drx\n", image);
		break;

	case 3 :
		dirax_sendline("dmax 1000\n", image);
		break;

	case 4 :
		dirax_sendline("indexfit 2\n", image);
		break;

	case 5 :
		dirax_sendline("levelfit 1000\n", image);
		break;

	case 6 :
		dirax_sendline("go\n", image);
		break;

	case 7 :
		dirax_sendline("acl\n", image);
		break;

	case 8 :
		if ( image->best_acl_nh == 0 ) {
			STATUS("No more cells to try.\n");
			dirax_sendline("exit", image);
			break;
		}
		snprintf(tmp, 31, "acl %i\n", image->best_acl);
		image->acls_tried[image->n_acls_tried++] = image->best_acl;
		dirax_sendline(tmp, image);
		break;

	case 9 :
		dirax_sendline("cell\n", image);
		break;

	case 10 :
		if ( image->n_acls_tried == MAX_DIRAX_CELL_CANDIDATES ) {
			STATUS("That's enough cells.\n");
			dirax_sendline("exit\n", image);
		} else {
			/* Go back round for another cell */
			image->best_acl_nh = 0;
			image->dirax_step = 6;
			dirax_sendline("acl\n", image);
		}
		break;

	default:
		dirax_sendline("exit\n", image);
		return;

	}

	image->dirax_step++;
}


static int dirax_readable(struct image *image)
{
	int rval;
	int no_string = 0;

	rval = read(image->dirax_pty, image->dirax_rbuffer+image->dirax_rbufpos,
	            image->dirax_rbuflen-image->dirax_rbufpos);

	if ( (rval == -1) || (rval == 0) ) return 1;

	image->dirax_rbufpos += rval;
	assert(image->dirax_rbufpos <= image->dirax_rbuflen);

	while ( (!no_string) && (image->dirax_rbufpos > 0) ) {

		int i;
		int block_ready = 0;
		DirAxInputType type = DIRAX_INPUT_NONE;

		/* See if there's a full line in the buffer yet */
		for ( i=0; i<image->dirax_rbufpos-1; i++ ) {
			/* Means the last value looked at is rbufpos-2 */

			/* Is there a prompt in the buffer? */
			if ( i+7 <= image->dirax_rbufpos ) {
				if ( (strncmp(image->dirax_rbuffer+i,
						"Dirax> ", 7) == 0)
				  || (strncmp(image->dirax_rbuffer+i,
				                "PROMPT:", 7) == 0)
				  || (strncmp(image->dirax_rbuffer+i,
				                "acl/auto:", 8) == 0) ) {
					block_ready = 1;
					type = DIRAX_INPUT_PROMPT;
					break;
				}
			}

			if ( (image->dirax_rbuffer[i] == '\r')
			  && (image->dirax_rbuffer[i+1] == '\n') ) {
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
				memcpy(block_buffer, image->dirax_rbuffer, i);
				block_buffer[i] = '\0';

				if ( block_buffer[0] == '\r' ) {
					memmove(block_buffer, block_buffer+1, i);
				}

				dirax_parseline(block_buffer, image);
				free(block_buffer);
				endbit_length = i+2;
				break;

			case DIRAX_INPUT_PROMPT :

				dirax_send_next(image);
				endbit_length = i+7;
				break;

			default :

				/* Obviously, this never happens :) */
				ERROR("Unrecognised DirAx input mode! "
				      "I don't know how to understand DirAx\n");
				return 1;

			}

			/* Now the block's been parsed, it should be
			 * forgotten about */
			memmove(image->dirax_rbuffer,
			        image->dirax_rbuffer + endbit_length,
			        image->dirax_rbuflen - endbit_length);

			/* Subtract the number of bytes removed */
			image->dirax_rbufpos = image->dirax_rbufpos
			                       - endbit_length;
			new_rbuflen = image->dirax_rbuflen - endbit_length;
			if ( new_rbuflen == 0 ) new_rbuflen = 256;
			image->dirax_rbuffer = realloc(image->dirax_rbuffer,
			                               new_rbuflen);
			image->dirax_rbuflen = new_rbuflen;

		} else {

			if ( image->dirax_rbufpos==image->dirax_rbuflen ) {

				/* More buffer space is needed */
				image->dirax_rbuffer = realloc(
				                    image->dirax_rbuffer,
				                    image->dirax_rbuflen + 256);
				image->dirax_rbuflen = image->dirax_rbuflen
				                                          + 256;
				/* The new space gets used at the next
				 * read, shortly... */

			}
			no_string = 1;

		}

	}

	return 0;
}


void run_dirax(struct image *image)
{
	unsigned int opts;
	int saved_stderr;
	int status;
	int rval;

	saved_stderr = dup(STDERR_FILENO);
	image->dirax_pid = forkpty(&image->dirax_pty, NULL, NULL, NULL);
	if ( image->dirax_pid == -1 ) {
		ERROR("Failed to fork for DirAx\n");
		return;
	}
	if ( image->dirax_pid == 0 ) {

		/* Child process: invoke DirAx */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		/* Reconnect stderr */
		dup2(saved_stderr, STDERR_FILENO);

		execlp("dirax", "", (char *)NULL);
		ERROR("Failed to invoke DirAx.\n");
		_exit(0);

	}

	image->dirax_rbuffer = malloc(256);
	image->dirax_rbuflen = 256;
	image->dirax_rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(image->dirax_pty, F_GETFL);
	fcntl(image->dirax_pty, F_SETFL, opts | O_NONBLOCK);

	image->dirax_step = 1;	/* This starts the "initialisation" procedure */
	image->dirax_read_cell = 0;
	image->n_acls_tried = 0;
	image->best_acl_nh = 0;

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(image->dirax_pty, &fds);

		tv.tv_sec = 10;
		tv.tv_usec = 0;

		sval = select(image->dirax_pty+1, &fds, NULL, NULL, &tv);

		if ( sval == -1 ) {
			ERROR("select() failed.\n");
			rval = 1;
		} else if ( sval != 0 ) {
			rval = dirax_readable(image);
		} else {
			ERROR("No response from DirAx..\n");
			rval = 1;
		}

	} while ( !rval );

	close(image->dirax_pty);
	close(saved_stderr);
	wait(&status);

	return;
}
