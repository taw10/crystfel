/*
 * felix.c
 *
 * Invoke Felix for multi-crystal autoindexing.
 *
 * Copyright © 2015-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2015 Thomas White <taw@physics.org>
 *   2015 Kenneth Beyerlein <kenneth.beyerlein@desy.de>
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
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
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
#include "felix.h"


#define FELIX_VERBOSE 0
#define MAX_MEAN_IA 0.25
#define MIN_NGV 15
#define MIN_NGV_UNIQUE 10


/* Global private data, prepared once */
struct felix_private
{
	IndexingMethod indm;
	UnitCell *cell;

	/* Options specific to Felix */
	int spacegroup;
	float tthrange_min;
	float tthrange_max;
	float etarange_min;
	float etarange_max;
	float domega;
	float omegarange_min;
	float omegarange_max;
	int min_measurements;   /* related to "cuts" */
	float min_completeness; /* related to "cuts" */
	float min_uniqueness;   /* related to "cuts" */
	int n_voxels;           /* related to "frustsumsize" */
	float test_fraction;    /* related to "frustsumsize" */
	float sigma_tth;        /* related to "uncertainties" */
	float sigma_eta;        /* related to "uncertainties" */
	float sigma_omega;      /* related to "uncertainties" */
	int n_sigmas;
	int force4frustums;

	/*Felix v0.3 options*/
	int orispace_frustum;
	int orispace_octa;
	char *readhkl_file;
	float maxtime;

};


/* Data needed to call Felix */
struct felix_data {

	struct felix_private *gp;

	/* Low-level stuff */
	int                     pty;
	pid_t                   pid;
	char                    *rbuffer;
	int                     rbufpos;
	int                     rbuflen;

};


static int read_felix(struct felix_private *gp, struct image *image,
                       char *filename)
{
	FILE *fh;
	int d1;
	float d2;
	float ubi11, ubi12, ubi13;
	float ubi21, ubi22, ubi23;
	float ubi31, ubi32, ubi33;
	float mean_ia;
	int ngv;
	char line[1024];
	int r;
	int n_crystals = 0;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Can't open '%s'\n", filename);
		return 0;
	}

	/* Read and discard first line */
	if ( fgets( line, 1024, fh ) == NULL ) {
		ERROR("Failed to read *.felix file.\n");
		return 0;
	}

	do {

		Crystal *cr;
		UnitCell *cell;

		/* One line per grain */
		if ( fgets( line, 1024, fh ) == NULL ) {
			break;
		}
		/* File format of the .felix files
		 * version 0.1 - 0.2
		 *
		 * r = sscanf(line, "%i %f %i %i %f %f %f %f %f %f %f %f %f"
		 *	         "%f %f %f %f %f %f %f %f %f %f %f %f",
		 *	         &d1, &mean_ia, &ngv, &ngv_unique, &d2, &d2, &d2,
		 *	         &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2, &d2,
		 *         &ubi11, &ubi12, &ubi13,
		 *	         &ubi21, &ubi22, &ubi23,
		 *	         &ubi31, &ubi32, &ubi33);
		 *
		 * if ( r != 25 ) {
		 *		ERROR("Only %i parameters in .felix file\n", r);
		 *		return 1;
		 * }
		 */

		/* version 0.3 - present */
		r = sscanf(line, "%i %f %i %f %f %f %f %f %f"
					 "%f %f %f %f %f %f %f %f %f %f %f %f",
					 &d1, &mean_ia, &ngv, &d2, &d2,
					 &d2, &d2, &d2, &d2, &d2, &d2, &d2,
					 &ubi11, &ubi12, &ubi13,
					 &ubi21, &ubi22, &ubi23,
					 &ubi31, &ubi32, &ubi33);

		if ( r != 21 ) {
			ERROR("Only %i parameters in .felix file, "
			      "check version and format.\n", r);
			return -1;
		}

		cell = cell_new();

		cell_set_cartesian(cell, ubi12/1e10, ubi13/1e10, ubi11/1e10,
			                 ubi22/1e10, ubi23/1e10, ubi21/1e10,
			                 ubi32/1e10, ubi33/1e10, ubi31/1e10);
		cell_set_lattice_type(cell, cell_get_lattice_type(gp->cell));
		cell_set_centering(cell, cell_get_centering(gp->cell));
		cell_set_unique_axis(cell, cell_get_unique_axis(gp->cell));

		cr = crystal_new();
		if ( cr == NULL ) {
			ERROR( "Failed to allocate crystal.\n" );
			return 0;
		}

		crystal_set_cell(cr, cell);

		/* Poor indexing criterion for Felix v0.1
		 *
		 * if (mean_ia > MAX_MEAN_IA || ngv < MIN_NGV ||
		 * 		ngv_unique < MIN_NGV_UNIQUE ){
		 * 		crystal_set_user_flag(cr, 1);
		 */

		/* Poor indexing criterion for Felix v0.3 */

		if ( mean_ia > MAX_MEAN_IA || ngv < MIN_NGV ){
			crystal_set_user_flag(cr, 1);
		}

		/* All crystals are saved to the image,
		 * but only good ones will be later written to the stream file.
		 */

		image_add_crystal(image, cr);
		n_crystals++;

	} while ( !feof(fh) );

	fclose(fh);

	return n_crystals;
}


static void gs_parseline(char *line, struct image *image,
                         struct felix_data *gs)
{
	#if FELIX_VERBOSE
	STATUS("%s\n", line);
	#endif
}


static int felix_readable(struct image *image, struct felix_data *gs)
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
			  && ( gs->rbuffer[i+1] == '\n' ) ) {
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

			if ( gs->rbufpos == gs->rbuflen ) {

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


static void write_gve(struct image *image, struct felix_private *gp)
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


static char *write_ini(struct image *image, struct felix_private *gp)
{
	FILE *fh;
	char *filename;
	char gveFilename[1024];
	char logFilename[1024];
	double tt;

	filename = malloc(1024);
	if ( filename == NULL ) return NULL;

	snprintf(filename, 1023, "xfel-%i.ini", image->id);
	snprintf(gveFilename, 1023, "xfel-%i.gve", image->id);
	snprintf(logFilename, 1023, "xfel-%i.log", image->id);

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

	fprintf(fh, "spacegroup %i\n", gp->spacegroup);
	fprintf(fh, "tthrange %f %f\n", rad2deg(gp->tthrange_min),
	                                rad2deg(gp->tthrange_max));
	fprintf(fh, "etarange %f %f\n", gp->etarange_min, gp->etarange_max);
	fprintf(fh, "domega %f\n", gp->domega);
	fprintf(fh, "omegarange %f %f\n", gp->omegarange_min, gp->omegarange_max);
	fprintf(fh, "filespecs %s %s\n", gveFilename, logFilename);
	fprintf(fh, "cuts %i %f %f\n", gp->min_measurements, gp->min_completeness,
			gp->min_uniqueness);
	fprintf(fh, "frustumsize %i %f\n", gp->n_voxels,
			gp->test_fraction);
	fprintf(fh, "uncertainties %f %f %f\n", gp->sigma_tth,
			gp->sigma_eta, gp->sigma_omega);
	fprintf(fh, "nsigmas %i\n", gp->n_sigmas);

	if ( gp->force4frustums == 1 ){
		fprintf(fh, "force4frustums\n");
	}

	if ( gp->orispace_frustum == 1 ){
		fprintf(fh, "orispace frustum\n");
	} else if ( gp->orispace_octa ==1 ){
		fprintf(fh, "orispace octa\n");
	} else{
		ERROR("No felix supported orispace specified.\n");
		free(filename);
		filename = NULL;
	}

	/* If an hkl file is not specified, generate the peak list. */
	if ( gp->readhkl_file != NULL ){
		fprintf(fh, "readhkl %s\n", gp->readhkl_file);
	} else{
		fprintf(fh, "genhkl\n");
	}

	fprintf(fh, "maxtime %f\n", gp->maxtime);

	fclose(fh);

	return filename;
}


int felix_index(struct image *image, IndexingPrivate *ipriv)
{
	unsigned int opts;
	int status;
	int rval;
	struct felix_data *felix;
	struct felix_private *gp = (struct felix_private *) ipriv;
	char *ini_filename;
	char gff_filename[1024];

	write_gve (image, gp);
	ini_filename = write_ini (image, gp);

	if ( ini_filename == NULL ) {
		ERROR("Failed to write ini file for Felix.\n");
		return 0;
	}

	felix = malloc(sizeof(struct felix_data));
	if ( felix == NULL ) {
		ERROR("Couldn't allocate memory for Felix data.\n");
		return 0;
	}

	felix->gp = gp;

	snprintf(gff_filename, 1023, "xfel-%i.felix", image->id);
	remove(gff_filename);

	felix->pid = forkpty(&felix->pty, NULL, NULL, NULL);
	if ( felix->pid == -1 ) {
		ERROR("Failed to fork for Felix: %s\n", strerror(errno));
		return 0;
	}
	if ( felix->pid == 0 ) {

		/* Child process: invoke Felix */
		struct termios t;

		/* Turn echo off */
		tcgetattr(STDIN_FILENO, &t);
		t.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
		tcsetattr(STDIN_FILENO, TCSANOW, &t);

		STATUS("Running Felix '%s'\n", ini_filename);
		execlp("Felix", "Felix", ini_filename, (char *)NULL);
		ERROR("Failed to invoke Felix.\n");
		_exit(0);

	}

	free(ini_filename);

	felix->rbuffer = malloc(256);
	felix->rbuflen = 256;
	felix->rbufpos = 0;

	/* Set non-blocking */
	opts = fcntl(felix->pty, F_GETFL);
	fcntl(felix->pty, F_SETFL, opts | O_NONBLOCK);

	do {

		fd_set fds;
		struct timeval tv;
		int sval;

		FD_ZERO(&fds);
		FD_SET(felix->pty, &fds);

		tv.tv_sec = gp->maxtime + 1.0;
		tv.tv_usec = 0;

		sval = select(felix->pty+1, &fds, NULL, NULL, &tv);

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
			rval = felix_readable(image, felix);
		} else {
			ERROR("No response from Felix..\n");
			rval = 1;
		}

	} while ( !rval );

	close(felix->pty);
	free(felix->rbuffer);
	waitpid(felix->pid, &status, 0);

	if ( status != 0 ) {
		ERROR("Felix either timed out, or is not working properly.\n");
		free(felix);
		return 0;
	}

	rval = read_felix(gp, image, gff_filename);

	free(felix);
	return rval;

}


static int sg_number_for_cell(UnitCell *cell)
{
	LatticeType lattice = cell_get_lattice_type(cell);
	char cen = cell_get_centering(cell);

	switch (lattice)
	{
		case L_TRICLINIC:
		return 1;  /* P1 */

		case L_MONOCLINIC:
		switch ( cen ) {
			case 'P' : return 3;  /* P2 */
			case 'C' : return 5;  /* C2 */
			default : return 0;
		}

		case L_ORTHORHOMBIC:
		switch ( cen ) {
			case 'P' : return 16;  /* P222 */
			case 'C' : return 21;  /* C222 */
			case 'F' : return 22;  /* F222 */
			case 'I' : return 23;  /* I222 */
			case 'A' : return 38;  /* Amm2 */
			default : return 0;
		}

		case L_TETRAGONAL:
		switch ( cen ) {
			case 'P' : return 89;  /* P422 */
			case 'I' : return 97;  /* I422 */
			default : return 0;
		}

		case L_RHOMBOHEDRAL:
		return 155;  /* R32 */

		case L_HEXAGONAL:
		switch ( cen ) {
			case 'P' : return 177;  /* P622 */
			case 'H' : return 143;  /* P3 */
			default : return 0;
		}

		case L_CUBIC:
		switch ( cen ) {
			case 'P' : return 207;  /* P432 */
			case 'F' : return 209;  /* F432 */
			case 'I' : return 211;  /* I432 */
			default : return 0;
		}

		default:
		return 0;
	}
}


void *felix_prepare(IndexingMethod *indm, UnitCell *cell,
                    struct felix_options *opts)
{
	struct felix_private *gp;

	if ( felix_probe(cell) == NULL ) {
		ERROR("Felix does not appear to run properly.\n");
		ERROR("Please check your Felix installation.\n");
		return NULL;
	}

	if ( !cell_has_parameters(cell) ) {
		ERROR("Felix needs a unit cell.\n");
		return NULL;
	}

	gp = calloc(1, sizeof(*gp));
	if ( gp == NULL ) return NULL;

	/* Flags that Felix knows about */
	*indm &= INDEXING_METHOD_MASK
	       | INDEXING_USE_LATTICE_TYPE | INDEXING_USE_CELL_PARAMETERS;

	gp->cell = cell;
	gp->indm = *indm;

	/* Default values of felix options */
	gp->spacegroup = sg_number_for_cell(cell);
	if ( gp->spacegroup == 0 ) {
		ERROR("Couldn't determine representative space group for your cell.\n");
		ERROR("Try again with a more conventional cell.\n");
		return NULL;
	}

	/* Default parameters */
	gp->n_voxels = 100;
	gp->etarange_min = 0;
	gp->etarange_max = 360;
	gp->domega = 2;
	gp->omegarange_min = -1.0;
	gp->omegarange_max = 1.0;
	gp->min_measurements = 15;
	gp->min_completeness = 0.001;
	gp->min_uniqueness = 0.5;
	gp->test_fraction = 0.75;
	gp->sigma_tth = 0.2;
	gp->sigma_eta = 0.2;
	gp->sigma_omega = 0.2;
	gp->n_sigmas = 1;
	gp->force4frustums = 0;
	gp->orispace_frustum = 1;
	gp->orispace_octa = 0;
	gp->readhkl_file = NULL;
	gp->maxtime = 120.0;
	gp->tthrange_min = deg2rad(0.0);
	gp->tthrange_max = deg2rad(30.0);

	if ( opts->ttmin > 0.0 ) {
		gp->tthrange_min = opts->ttmin;
	}
	if ( opts->ttmax > 0.0 ) {
		gp->tthrange_max = opts->ttmax;
	}
	if ( opts->min_measurements > 0 ) {
		gp->min_measurements = opts->min_measurements;
	}
	if ( opts->min_completeness > 0.0 ) {
		gp->min_completeness = opts->min_completeness;
	}
	if ( opts->min_uniqueness > 0.0 ) {
		gp->min_uniqueness = opts->min_uniqueness;
	}
	if ( opts->n_voxels > 0 ) {
		gp->n_voxels = opts->n_voxels;
	}
	if ( opts->test_fraction > 0.0 ) {
		gp->test_fraction = opts->test_fraction;
	}
	if ( opts->sigma > 0.0 ) {
		gp->sigma_tth = opts->sigma;
		gp->sigma_eta = opts->sigma;
		gp->sigma_omega = opts->sigma;
	}
	if (opts->domega > 0.0 ) {
		gp->domega = opts -> domega;
	}

	return (IndexingPrivate *)gp;
}


void felix_cleanup(IndexingPrivate *pp)
{
	struct felix_private *p;

	p = (struct felix_private *) pp;
	free(p->readhkl_file);
	free(p);
}


const char *felix_probe(UnitCell *cell)
{
	pid_t pid;
	int pty;
	int status;
	FILE *fh;
	char line[1024];
	int ok = 0;

	if ( !cell_has_parameters(cell) ) {
		return NULL;
	}

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

		execlp("Felix", "Felix", (char *)NULL);
		_exit(1);

	}

	fh = fdopen(pty, "r");
	if ( fgets(line, 1024, fh) == NULL ) {
		ok = 0;
	} else {
		if ( strncmp(line, "Felix", 5) == 0 ) {
			ok = 1;
		}
	}

	fclose(fh);
	close(pty);
	waitpid(pid, &status, 0);

	if ( ok ) return "felix";
	return NULL;
}
