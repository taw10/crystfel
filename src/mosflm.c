/*
 * mosflm.c
 *
 * Invoke the DPS auto-indexing algorithm through MOSFLM
 *
 * This will actuaully run DirAx for now... will be fixed soon...
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
#include "mosflm.h"
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
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &ax, &ay, &az);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line\n");
			image->dirax_read_cell = 0;
			free(image->candidate_cells[image->ncells]);
			return;
		}
		cell_set_cartesian_a(image->candidate_cells[image->ncells],
		                     ax*1e-10, ay*1e-10, az*1e-10);
		image->dirax_read_cell++;
		return;
	} else if ( image->dirax_read_cell == 2 ) {
		/* Second row of unit cell values */
		float bx, by, bz;
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &bx, &by, &bz);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line\n");
			image->dirax_read_cell = 0;
			free(image->candidate_cells[image->ncells]);
			return;
		}
		cell_set_cartesian_b(image->candidate_cells[image->ncells],
		                     bx*1e-10, by*1e-10, bz*1e-10);
		image->dirax_read_cell++;
		return;
	} else if ( image->dirax_read_cell == 3 ) {
		/* Third row of unit cell values */
		float cx, cy, cz;
		int r;
		r = sscanf(line, "%f %f %f %f %f %f",
		           &d, &d, &d, &cx, &cy, &cz);
		if ( r != 6 ) {
			ERROR("Couldn't understand cell line\n");
			image->dirax_read_cell = 0;
			free(image->candidate_cells[image->ncells]);
			return;
		}
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
		snprintf(tmp, 31, "read xfel-%i.drx\n", image->id);
		dirax_sendline(tmp, image);
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
			/* At this point, DirAx is presenting its ACL prompt
			 * and waiting for a single number.  Use an extra
			 * newline to choose automatic ACL selection before
			 * exiting. */
			dirax_sendline("\nexit\n", image);
			break;
		}
		snprintf(tmp, 31, "%i\n", image->best_acl);
		image->acls_tried[image->n_acls_tried++] = image->best_acl;
		dirax_sendline(tmp, image);
		break;

	case 9 :
		dirax_sendline("cell\n", image);
		break;

	case 10 :
		if ( image->n_acls_tried == MAX_DIRAX_CELL_CANDIDATES ) {
			dirax_sendline("exit\n", image);
		} else {
			/* Go back round for another cell */
			image->best_acl_nh = 0;
			image->dirax_step = 7;
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
			if ( (i+7 <= image->dirax_rbufpos)
			  && (!strncmp(image->dirax_rbuffer+i, "Dirax> ", 7)) ) {
				block_ready = 1;
				type = DIRAX_INPUT_PROMPT;
				break;
			}

			/* How about an ACL prompt? */
			if ( (i+10 <= image->dirax_rbufpos)
			  && (!strncmp(image->dirax_rbuffer+i, "acl/auto [", 10)) ) {
				block_ready = 1;
				type = DIRAX_INPUT_ACL;
				break;
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

			case DIRAX_INPUT_ACL :

				dirax_send_next(image);
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


static int read_newmat(const char * filename, struct image *image)
{
	FILE * fh;
	float asx, asy, asz;
	float bsx, bsy, bsz;
	float csx, csy, csz;
	int n;
	double c;
	
	fh = fopen(filename,"r");
	if (fh == NULL){ 
		return 1;
	}
	n  = fscanf(fh,"%f %f %f\n",&asx,&bsx,&csx);
	n += fscanf(fh,"%f %f %f\n",&asy,&bsy,&csy);
	n += fscanf(fh,"%f %f %f\n",&asz,&bsz,&csz);
	if (n != 9) {
		return 1;
	}
	fclose(fh);
	
	/* mosflm A matrix is multiplied by lambda, so fix this */
	c = 1e-10/image->lambda;
	
	cell_set_reciprocal(image->candidate_cells[0],
                            asz*c, asy*c, asx*c,
                            bsz*c, bsy*c, bsx*c,
                            csz*c, csy*c, csx*c);
               
        image->ncells = 1;
        
        return 0;
}


void run_mosflm(struct image *image, UnitCell *cell)
{
	unsigned int opts;
	int status;
	int rval;
	int i,j;
	char mos_cmd[1024];
	char symm[64];
	const char *sg;
	double a,b,c,alpha,beta,gamma;
	double wavelength; /* angstrom */
	const char newmatfile[128];
	int fail;
	
	printf("Mosflm is not fully implemented.  Using DirAx for now.\n");
	
	wavelength = image->lambda*1e10;
	cell_get_parameters(cell, &a, &b, &c, &alpha, &beta, &gamma);
	sg = cell_get_spacegroup(cell);
	sprintf(newmatfile,"xfel-%i.newmat",image->id);
	
	/* need to remove white space from spacegroup... */
	j = 0;
	for(i = 0; i < strlen(sg);i++)
	{
		if (sg[i] != ' ') {
			symm[j] = sg[i];
			j++;
		}
	}
	symm[j] = '\0';	
	
	
	/* build a script to run mosflm */
	sprintf(mos_cmd,"%s","ipmosflm << eof-mosflm >> /dev/null\n");
	sprintf(mos_cmd,"%s%s",mos_cmd,
	                     "DETECTOR ROTATION HORIZONTAL ANTICLOCKWISE"
	                     " ORIGIN LL FAST HORIZONTAL RECTANGULAR\n");
	sprintf(mos_cmd,"%sCELL %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
	                   mos_cmd,
	                   a*1e10,b*1e10,c*1e10,
	                   rad2deg(alpha),rad2deg(beta),rad2deg(gamma));
	sprintf(mos_cmd,"%sSYMM %s\n",mos_cmd,symm);
	sprintf(mos_cmd,"%sDISTANCE %8.4f\n",mos_cmd,67.8);
	sprintf(mos_cmd,"%sBEAM %8.4f %8.4f\n",mos_cmd,0.0,0.0);
	sprintf(mos_cmd,"%sWAVELENGTH %10.5f\n",mos_cmd,wavelength);
	sprintf(mos_cmd,"%sNEWMAT %s\n",mos_cmd,newmatfile);
	sprintf(mos_cmd,"%sIMAGE xfel-%i_001.img phi 0 0\n",mos_cmd,image->id);
	sprintf(mos_cmd,"%sAUTOINDEX DPS FILE xfel-%i.spt IMAGE 1\n",
	                      mos_cmd,image->id);
	sprintf(mos_cmd,"%sGO\n",mos_cmd);
	sprintf(mos_cmd,"%s%s",mos_cmd,"eof-mosflm\n");

	/* Run the mosflm script */
	fail = system(mos_cmd);
	if (fail) { 
		ERROR("mosflm execution failed.\n");
		return; 
	}

	/* Read the mosflm NEWMAT file and set cell candidate */ 
	/* Existence of this file means possible success. Pretty shady. */
	fail = read_newmat(newmatfile,image);
	if (fail) {
		printf("Failed to read mosflm NEWMAT file.\n"); 
		return;
	}

	/* remove the mosflm NEWMAT file */
	remove(newmatfile);
	
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
	free(image->dirax_rbuffer);
	wait(&status);

	return;
}
