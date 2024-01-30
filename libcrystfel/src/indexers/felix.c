/*
 * felix.c
 *
 * Invoke Felix for multi-crystal autoindexing.
 *
 * Copyright © 2015-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2015-2021 Thomas White <taw@physics.org>
 *   2015      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
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

/** \file felix.h */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include <errno.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#ifdef HAVE_FORKPTY_PTY_H
#include <pty.h>
#endif
#ifdef HAVE_FORKPTY_UTIL_H
#include <util.h>
#endif

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell.h"
#include "cell-utils.h"
#include "felix.h"
#include "index.h"


#define FELIX_VERBOSE 0


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
	int min_visits;   /* related to "cuts" */
	float min_completeness; /* related to "cuts" */
	float max_uniqueness;   /* related to "cuts" */
	int n_voxels;           /* related to "frustsumsize" */
	float fraction_max_visits;    /* related to "frustsumsize" */
	float sigma_tth;        /* related to "uncertainties" */
	float sigma_eta;        /* related to "uncertainties" */
	float sigma_omega;      /* related to "uncertainties" */
	int n_sigmas;
	int force4frustums;
	float max_internal_angle;

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
	int n_crystals = 0;

	fh = fopen(filename, "r");
	if ( fh == NULL ) {
		ERROR("Can't open '%s'\n", filename);
		return 0;
	}

	/* Read and discard first line */
	if ( fgets( line, 1024, fh ) == NULL ) {
		ERROR("Failed to read *.felix file.\n");
		fclose(fh);
		return 0;
	}

	do {

		Crystal *cr;
		UnitCell *cell;
		int r;

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

		if ( mean_ia > gp->max_internal_angle ){
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

			block_buffer = cfmalloc(i+1);
			memcpy(block_buffer, gs->rbuffer, i);
			block_buffer[i] = '\0';

			if ( block_buffer[0] == '\r' ) {
				memmove(block_buffer, block_buffer+1, i);
			}

			gs_parseline(block_buffer, image, gs);
			cffree(block_buffer);
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
			gs->rbuffer = cfrealloc(gs->rbuffer, new_rbuflen);
			gs->rbuflen = new_rbuflen;

		} else {

			if ( gs->rbufpos == gs->rbuflen ) {

				/* More buffer space is needed */
				gs->rbuffer = cfrealloc(gs->rbuffer,
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
	snprintf(filename, 1023, "xfel.gve");
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
		double r[3];

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		detgeom_transform_coords(&image->detgeom->panels[f->pn],
		                         f->fs, f->ss, image->lambda,
		                         0.0, 0.0, r);

		fprintf(fh, "%.6f %.6f %.6f 0 0 %.6f %.6f %.6f 0\n",
		        r[2]/1e10, r[0]/1e10, r[1]/1e10,
		        modulus(r[0], r[1], r[2])/1e10, /* dstar */
		        rad2deg(atan2(r[1], r[0])), 0.0);   /* eta, omega */

	}
	fclose(fh);
}


static char *write_ini(struct image *image, struct felix_private *gp)
{
	FILE *fh;
	char *filename;
	char gveFilename[1024];
	char logFilename[1024];

	filename = cfmalloc(1024);
	if ( filename == NULL ) return NULL;

	snprintf(filename, 1023, "xfel.ini");
	snprintf(gveFilename, 1023, "xfel.gve");
	snprintf(logFilename, 1023, "xfel.log");

	fh = fopen(filename, "w");
	if ( !fh ) {
		ERROR("Couldn't open temporary file '%s'\n", filename);
		cffree(filename);
		return NULL;
	}

	fprintf(fh, "spacegroup %i\n", gp->spacegroup);
	fprintf(fh, "tthrange %f %f\n", rad2deg(gp->tthrange_min),
	                                rad2deg(gp->tthrange_max));
	fprintf(fh, "etarange %f %f\n", gp->etarange_min, gp->etarange_max);
	fprintf(fh, "domega %f\n", gp->domega);
	fprintf(fh, "omegarange %f %f\n", gp->omegarange_min, gp->omegarange_max);
	fprintf(fh, "filespecs %s %s\n", gveFilename, logFilename);
	fprintf(fh, "cuts %i %f %f\n", gp->min_visits, gp->min_completeness,
			gp->max_uniqueness);
	fprintf(fh, "frustumsize %i %f\n", gp->n_voxels,
			gp->fraction_max_visits);
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
		cffree(filename);
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

	felix = cfmalloc(sizeof(struct felix_data));
	if ( felix == NULL ) {
		ERROR("Couldn't allocate memory for Felix data.\n");
		return 0;
	}

	felix->gp = gp;

	snprintf(gff_filename, 1023, "xfel.felix");
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

	cffree(ini_filename);

	felix->rbuffer = cfmalloc(256);
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
	cffree(felix->rbuffer);
	waitpid(felix->pid, &status, 0);

	if ( status != 0 ) {
		ERROR("Felix either timed out, or is not working properly.\n");
		cffree(felix);
		return 0;
	}

	rval = read_felix(gp, image, gff_filename);

	cffree(felix);
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

	if ( !cell_has_parameters(cell) ) {
		ERROR("Felix needs a unit cell.\n");
		return NULL;
	}

	if ( felix_probe(cell) == NULL ) {
		ERROR("Felix does not appear to run properly.\n");
		ERROR("Please check your Felix installation.\n");
		return NULL;
	}

	gp = cfcalloc(1, sizeof(*gp));
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
		cffree(gp);
		return NULL;
	}

	/* Default parameters */
	gp->n_voxels = 100;
	gp->etarange_min = 0;
	gp->etarange_max = 360;
	gp->domega = 2;
	gp->omegarange_min = -1.0;
	gp->omegarange_max = 1.0;
	gp->min_visits = 15;
	gp->min_completeness = 0.001;
	gp->max_uniqueness = 0.5;
	gp->fraction_max_visits = 0.75;
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
	gp->max_internal_angle = 0.25;

	if ( opts->ttmin > 0.0 ) {
		gp->tthrange_min = opts->ttmin;
	}
	if ( opts->ttmax > 0.0 ) {
		gp->tthrange_max = opts->ttmax;
	}
	if ( opts->min_visits > 0 ) {
		gp->min_visits = opts->min_visits;
	}
	if ( opts->min_completeness > 0.0 ) {
		gp->min_completeness = opts->min_completeness;
	}
	if ( opts->max_uniqueness > 0.0 ) {
		gp->max_uniqueness = opts->max_uniqueness;
	}
	if ( opts->n_voxels > 0 ) {
		gp->n_voxels = opts->n_voxels;
	}
	if ( opts->fraction_max_visits > 0.0 ) {
		gp->fraction_max_visits = opts->fraction_max_visits;
	}
	if ( opts->sigma > 0.0 ) {
		gp->sigma_tth = opts->sigma;
		gp->sigma_eta = opts->sigma;
		gp->sigma_omega = opts->sigma;
	}
	if ( opts->domega > 0.0 ) {
		gp->domega = opts->domega;
	}
	if ( opts->max_internal_angle > 0 ) {
		gp->max_internal_angle = opts->max_internal_angle;
	}

	return (IndexingPrivate *)gp;
}


void felix_cleanup(IndexingPrivate *pp)
{
	struct felix_private *p;

	p = (struct felix_private *) pp;
	cffree(p->readhkl_file);
	cffree(p);
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

	/* Felix will write gmon.out when we test it, which we are
	 * are going to delete afterwards.  Better check the file doesn't exist
	 * first, in case it was important. */
	if ( file_exists("gmon.out") ) {
		ERROR("Please move or delete gmon.out from the working "
		      "directory first.\n");
		exit(1);
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

	unlink("gmon.out");

	if ( ok ) return "felix";
	return NULL;
}


static void felix_show_help()
{
	printf("Parameters for the Felix indexing algorithm:\n"
"     --felix-domega        Degree range of omega (moscaicity) to consider.\n"
"                            Default: 2\n"
"     --felix-fraction-max-visits\n"
"                           Cutoff for minimum fraction of the max visits.\n"
"                            Default: 0.75\n"
"     --felix-max-internal-angle\n"
"                           Cutoff for maximum internal angle between observed\n"
"                            spots and predicted spots. Default: 0.25\n"
"     --felix-max-uniqueness\n"
"                           Cutoff for maximum fraction of found spots which\n"
"                            can belong to other crystallites.  Default: 0.5\n"
"     --felix-min-completeness\n"
"                           Cutoff for minimum fraction of projected spots\n"
"                            found in the pattern. Default: 0.001\n"
"     --felix-min-visits\n"
"                           Cutoff for minimum number of voxel visits.\n"
"                            Default: 15\n"
"     --felix-num-voxels    Number of voxels for Rodrigues space search\n"
"                            Default: 100\n"
"     --felix-sigma         The sigma of the 2theta, eta and omega angles.\n"
"                            Default: 0.2\n"
"     --felix-tthrange-max  Maximum 2theta to consider for indexing (degrees)\n"
"                            Default: 30\n"
"     --felix-tthrange-min  Minimum 2theta to consider for indexing (degrees)\n"
"                            Default: 0\n"
);
}


int felix_default_options(struct felix_options **opts_ptr)
{
	struct felix_options *opts;

	opts = cfmalloc(sizeof(struct felix_options));
	if ( opts == NULL ) return ENOMEM;

	opts->ttmin = -1.0;
	opts->ttmax = -1.0;
	opts->min_visits = 0;
	opts->min_completeness = -1.0;
	opts->max_uniqueness = -1.0;
	opts->n_voxels = 0;
	opts->fraction_max_visits = -1.0;
	opts->sigma = -1.0;
	opts->domega = -1.0;
	opts->max_internal_angle = -1.0;

	*opts_ptr = opts;
	return 0;
}


static error_t felix_parse_arg(int key, char *arg,
                               struct argp_state *state)
{
	struct felix_options **opts_ptr = state->input;
	float tmp;
	int r;

	switch ( key ) {

		case ARGP_KEY_INIT :
		r = felix_default_options(opts_ptr);
		if ( r ) return r;
		break;

		case 1 :
		felix_show_help();
		return EINVAL;

		case 2 :
		if ( sscanf(arg, "%f", &tmp) != 1 ) {
			ERROR("Invalid value for --felix-tthrange-min\n");
			return EINVAL;
		}
		(*opts_ptr)->ttmin = deg2rad(tmp);
		break;

		case 3 :
		if ( sscanf(arg, "%f", &tmp) != 1 ) {
			ERROR("Invalid value for --felix-tthrange-max\n");
			return EINVAL;
		}
		(*opts_ptr)->ttmax = deg2rad(tmp);
		break;

		case 4 :
		if ( sscanf(arg, "%d", &(*opts_ptr)->min_visits) != 1 ) {
			ERROR("Invalid value for --felix-min-visits\n");
			return EINVAL;
		}
		break;

		case 5 :
		if ( sscanf(arg, "%lf", &(*opts_ptr)->min_completeness) != 1 ) {
			ERROR("Invalid value for --felix-min-completeness\n");
			return EINVAL;
		}
		break;

		case 6 :
		if ( sscanf(arg, "%lf", &(*opts_ptr)->max_uniqueness) != 1 ) {
			ERROR("Invalid value for --felix-max-uniqueness\n");
			return EINVAL;
		}
		break;

		case 7 :
		if ( sscanf(arg, "%d", &(*opts_ptr)->n_voxels) != 1 ) {
			ERROR("Invalid value for --felix-num-voxels\n");
			return EINVAL;
		}
		break;

		case 8 :
		if ( sscanf(arg, "%lf", &(*opts_ptr)->fraction_max_visits) != 1 ) {
			ERROR("Invalid value for --felix-fraction-max-visits\n");
			return EINVAL;
		}
		break;

		case 9 :
		if ( sscanf(arg, "%lf", &(*opts_ptr)->sigma) != 1 ) {
			ERROR("Invalid value for --felix-sigma\n");
			return EINVAL;
		}
		break;

		case 10 :
		if ( sscanf(arg, "%lf", &(*opts_ptr)->domega) != 1 ) {
			ERROR("Invalid value for --felix-domega\n");
			return EINVAL;
		}
		break;

		case 11 :
		if ( sscanf(arg, "%lf", &(*opts_ptr)->max_internal_angle) != 1 ) {
			ERROR("Invalid value for --felix-max-internal-angle\n");
			return EINVAL;
		}
		break;

		default :
		return ARGP_ERR_UNKNOWN;

	}

	return 0;
}


static struct argp_option felix_options[] = {

	{"help-felix", 1, NULL, OPTION_NO_USAGE,
	 "Show options for Felix indexing algorithm", 99},
	{"felix-tthrange-min", 2, "2theta", OPTION_HIDDEN, NULL},
	{"felix-tthrange-max", 3, "2theta", OPTION_HIDDEN, NULL},
	{"felix-min-visits", 4, "n", OPTION_HIDDEN, NULL},
	{"felix-min-completeness", 5, "frac", OPTION_HIDDEN, NULL},
	{"felix-max-uniqueness", 6, "n", OPTION_HIDDEN, NULL},
	{"felix-num-voxels", 7, "n", OPTION_HIDDEN, NULL},
	{"felix-fraction-max-visits", 8, "n", OPTION_HIDDEN, NULL},
	{"felix-sigma", 9, "n", OPTION_HIDDEN, NULL},
	{"felix-domega", 10, "n", OPTION_HIDDEN, NULL},
	{"felix-max-internal-angle", 11, "ang", OPTION_HIDDEN, NULL},

	{0}
};


struct argp felix_argp = { felix_options, felix_parse_arg,
                           NULL, NULL, NULL, NULL, NULL };
