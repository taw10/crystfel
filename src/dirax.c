/*
 * dirax.c
 *
 * Invoke the DirAx auto-indexing program
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <glib.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <pty.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/ioctl.h>
#include <termio.h>
#include <sgtty.h>


#include "image.h"
#include "dirax.h"
#include "utils.h"


typedef enum {
	DIRAX_INPUT_NONE,
	DIRAX_INPUT_LINE,
	DIRAX_INPUT_PROMPT
} DirAxInputType;


static void dirax_parseline(const char *line, struct image *image)
{
	int i, rf;
	char *copy;

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='r';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("DirAx: %s\n", copy);
	free(copy);

	if ( strstr(line, "reflections from file") ) {
		ERROR("DirAx can't understand this data.");
		return;
	}

	/* Is this the first line of a unit cell specification? */
	rf = 0; i = 0;
	while ( (i<strlen(line)) && ((line[i] == 'R')
				|| (line[i] == 'D') || (line[i] == ' ')) ) {
		if ( line[i] == 'R' ) rf = 1;
		if ( (line[i] == 'D') && rf ) {
			image->dirax_read_cell = 1;
			if ( image->cell ) {
				free(image->cell);
			}
			image->cell = cell_new();
			return;
		}
		i++;
	}

	/* Parse unit cell vectors as appropriate */
	if ( image->dirax_read_cell == 1 ) {
		/* First row of unit cell values */
		float x1, x2, x3;
		sscanf(line, "%f %f %f", &x1, &x2, &x3);
		cell_set_cartesian_x(image->cell, x1*1e10, x2*1e10, x3*1e10);
		image->dirax_read_cell++;
		return;
	} else if ( image->dirax_read_cell == 2 ) {
		/* First row of unit cell values */
		float y1, y2, y3;
		sscanf(line, "%f %f %f", &y1, &y2, &y3);
		cell_set_cartesian_y(image->cell, y1*1e10, y2*1e10, y3*1e10);
		image->dirax_read_cell++;
		return;
	} else if ( image->dirax_read_cell == 3 ) {
		/* First row of unit cell values */
		float z1, z2, z3;
		sscanf(line, "%f %f %f", &z1, &z2, &z3);
		cell_set_cartesian_z(image->cell, z1*1e10, z2*1e10, z3*1e10);
		STATUS("Read a reciprocal unit cell from DirAx\n");
		/* FIXME: Do something */
		image->dirax_read_cell = 0;
		return;
	}

	image->dirax_read_cell = 0;
}


static void dirax_sendline(const char *line, struct image *image)
{
	char *copy;
	int i;

	write(image->dirax_pty, line, strlen(line));

	copy = strdup(line);
	for ( i=0; i<strlen(copy); i++ ) {
		if ( copy[i] == '\r' ) copy[i]='\0';
		if ( copy[i] == '\n' ) copy[i]='\0';
	}
	STATUS("To DirAx: '%s'\n", copy);
	free(copy);
}


static void dirax_send_next(struct image *image)
{
	switch ( image->dirax_step ) {

		case 1 : {
			dirax_sendline("\\echo off\n", image);
			image->dirax_step++;
			break;
		}

		case 2 : {
			dirax_sendline("read xfel.drx\n", image);
			image->dirax_step++;
			break;
		}

		case 3 : {
			dirax_sendline("dmax 10\n", image);
			image->dirax_step++;
			break;
		}

		case 4 : {
			dirax_sendline("indexfit 2\n", image);
			image->dirax_step++;
			break;
		}

		case 5 : {
			dirax_sendline("levelfit 200\n", image);
			image->dirax_step++;
			break;
		}

		case 6 : {
			dirax_sendline("go\n", image);
			image->dirax_step++;
			break;
		}

		case 7 : {
			dirax_sendline("cell\n", image);
			image->dirax_step++;
			break;
		}

		default: {
			image->dirax_step = 0;
			STATUS("DirAx is idle\n");
		}

	}
}


static gboolean dirax_readable(GIOChannel *dirax, GIOCondition condition,
                               struct image *image)
{
	int rval;

	rval = read(image->dirax_pty, image->dirax_rbuffer+image->dirax_rbufpos,
			image->dirax_rbuflen-image->dirax_rbufpos);

	if ( (rval == -1) || (rval == 0) ) {

		ERROR("Lost connection to DirAx\n");
		waitpid(image->dirax_pid, NULL, 0);
		g_io_channel_shutdown(image->dirax, FALSE, NULL);
		image->dirax = NULL;
		return FALSE;

	} else {

		int no_string = 0;

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
					  		"PROMPT:", 7) == 0) ) {
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

				switch ( type ) {

					case DIRAX_INPUT_LINE : {

						char *block_buffer = NULL;

						block_buffer = malloc(i+1);
						memcpy(block_buffer,
							image->dirax_rbuffer, i);
						block_buffer[i] = '\0';

						if ( block_buffer[0] == '\r' ) {
							memmove(block_buffer,
							  block_buffer+1, i);
						}

						dirax_parseline(block_buffer,
									image);
						free(block_buffer);
						endbit_length = i+2;

						break;

					}

					case DIRAX_INPUT_PROMPT : {

						dirax_send_next(image);
						endbit_length = i+7;
						break;

					}

					default : {
						ERROR(
			" Unrecognised input mode (this never happens!)\n");
						abort();
					}

				}

				/* Now the block's been parsed, it should be
				 * forgotten about */
				memmove(image->dirax_rbuffer,
					image->dirax_rbuffer + endbit_length,
					image->dirax_rbuflen - endbit_length);

				/* Subtract the number of bytes removed */
				image->dirax_rbufpos = image->dirax_rbufpos
								- endbit_length;
				new_rbuflen = image->dirax_rbuflen
								- endbit_length;
				if ( new_rbuflen == 0 ) {
					new_rbuflen = 256;
				}
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

	}

	return TRUE;
}


static int map_position(struct image *image, double x, double y,
                        double *rx, double *ry, double *rz)
{
	/* "Input" space */
	double d;

	/* Angular description of reflection */
	double theta, psi, k;

	x -= image->x_centre;
	y -= image->y_centre;
	k = 1.0 / image->lambda;

	/* FIXME: Don't process lower CCD for now */
	if ( y < 0 ) return 0;

	if ( image->fmode == FORMULATION_CLEN ) {

		/* Convert pixels to metres */
		x /= image->resolution;
		y /= image->resolution;	/* Convert pixels to metres */
		d = sqrt((x*x) + (y*y));
		theta = atan2(d, image->camera_len);

	} else if (image->fmode == FORMULATION_PIXELSIZE ) {

		/* Convert pixels to metres^-1 */
		x *= image->pixel_size;
		y *= image->pixel_size;	/* Convert pixels to metres^-1 */
		d = sqrt((x*x) + (y*y));
		theta = atan2(d, k);

	} else {
		ERROR("Unrecognised formulation mode in mapping_scale.\n");
		return -1;
	}

	psi = atan2(y, x);

	*rx = k*sin(theta)*cos(psi);
	*ry = k*sin(theta)*sin(psi);
	*rz = k - k*cos(theta);

	return 0;
}


#define PEAK_WINDOW_SIZE (10)

static void search_peaks(struct image *image)
{
	FILE *fh;
	int x, y, width, height;
	int16_t *data;

	fh = fopen("xfel.drx", "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file xfel.drx\n");
		return;
	}
	fprintf(fh, "%f\n", 0.5);  /* Lie about the wavelength.  */

	data = image->data;
	width = image->width;
	height = image->height;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	for ( x=1; x<image->width-1; x++ ) {
	for ( y=1; y<image->height-1; y++ ) {

		double dx1, dx2, dy1, dy2;
		double dxs, dys;
		double grad;
		int mask_x, mask_y;
		int sx, sy;
		double max;
		unsigned int did_something = 1;

		/* Overall threshold */
		if ( data[x+width*y] < 800 ) continue;

		/* Ignore streak */
		if ( abs(x-image->x_centre) < 15 ) continue;

		/* Get gradients */
		dx1 = data[x+width*y] - data[(x+1)+width*y];
		dx2 = data[(x-1)+width*y] - data[x+width*y];
		dy1 = data[x+width*y] - data[(x+1)+width*(y+1)];
		dy2 = data[x+width*(y-1)] - data[x+width*y];

		/* Average gradient measurements from both sides */
		dxs = ((dx1*dx1) + (dx2*dx2)) / 2;
		dys = ((dy1*dy1) + (dy2*dy2)) / 2;

		/* Calculate overall gradient */
		grad = dxs + dys;

		if ( grad < 2000000 ) continue;

		mask_x = x;
		mask_y = y;

		while ( (did_something)
		     && (distance(mask_x, mask_y, x, y)<50) ) {

			max = data[mask_x+width*mask_y];
			did_something = 0;

			for ( sy=biggest(mask_y-PEAK_WINDOW_SIZE/2, 0);
			      sy<smallest(mask_y+PEAK_WINDOW_SIZE/2, height);
			      sy++ ) {
			for ( sx=biggest(mask_x-PEAK_WINDOW_SIZE/2, 0);
			      sx<smallest(mask_x+PEAK_WINDOW_SIZE/2, width);
			      sx++ ) {

				if ( data[sx+width*sy] > max ) {
					max = data[sx+width*sy];
					mask_x = sx;
					mask_y = sy;
					did_something = 1;
				}

			}
			}

		}

		if ( !did_something ) {

			double d;
			int idx;

			assert(mask_x<image->width);
			assert(mask_y<image->height);
			assert(mask_x>=0);
			assert(mask_y>=0);

			/* Too far from foot point? */
			if ( distance(mask_x, mask_y, x, y) > 50.0 )  continue;

			/* Check for a feature at exactly the
			 * same coordinates */
			image_feature_closest(image->features, mask_x, mask_y,
			                      &d, &idx);

			if ( d > 1.0 ) {

				double rx = 0.0;
				double ry = 0.0;
				double rz = 0.0;

				/* Map and record reflection */
				printf("%i %i\n", x, y);

				image_add_feature(image->features,
				                  mask_x, mask_y, image, 1.0);
				map_position(image, x, y, &rx, &ry, &rz);
				fprintf(fh, "%10f %10f %10f %8f\n",
				        rx/1e10, ry/1e10, rz/1e10, 1.0);
			}

		}


	}
	}

	fclose(fh);
}


void index_pattern(struct image *image, int no_index)
{
	unsigned int opts;
	int saved_stderr;
	GMainLoop *ml;

	/* Do peak search and splurge out 'xfel.drx' */
	search_peaks(image);

	if ( no_index ) return;

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

	image->dirax = g_io_channel_unix_new(image->dirax_pty);
	g_io_add_watch(image->dirax, G_IO_IN | G_IO_HUP,
	               (GIOFunc)dirax_readable, image);

	ml = g_main_loop_new(NULL, FALSE);
	g_main_loop_run(ml);

	return;
}
