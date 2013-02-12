/*
 * stream.c
 *
 * Stream tools
 *
 * Copyright © 2013 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2011      Andrew Aquila
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
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "cell.h"
#include "utils.h"
#include "image.h"
#include "stream.h"
#include "reflist.h"
#include "reflist-utils.h"

#define LATEST_MAJOR_VERSION (2)
#define LATEST_MINOR_VERSION (1)


struct _stream
{
	FILE *fh;

	int major_version;
	int minor_version;
};


static int read_peaks(FILE *fh, struct image *image)
{
	char *rval = NULL;
	int first = 1;

	image->features = image_feature_list_new();

	do {

		char line[1024];
		float x, y, d, intensity;
		int r;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, PEAK_LIST_END_MARKER) == 0 ) return 0;

		r = sscanf(line, "%f %f %f %f", &x, &y, &d, &intensity);
		if ( (r != 4) && (!first) ) {
			ERROR("Failed to parse peak list line.\n");
			ERROR("The failed line was: '%s'\n", line);
			return 1;
		}

		first = 0;
		if ( r == 4 ) {
			image_add_feature(image->features, x, y,
			                  image, intensity, NULL);
		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return 1;
}


static void write_peaks(struct image *image, FILE *ofh)
{
	int i;

	fprintf(ofh, PEAK_LIST_START_MARKER"\n");
	fprintf(ofh, "  fs/px   ss/px  (1/d)/nm^-1   Intensity\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec r;
		double q;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);
		q = modulus(r.u, r.v, r.w);

		fprintf(ofh, "%7.2f %7.2f   %10.2f  %10.2f\n",
		       f->fs, f->ss, q/1.0e9, f->intensity);

	}

	fprintf(ofh, PEAK_LIST_END_MARKER"\n");
}


static void write_crystal(Stream *st, Crystal *cr, int include_reflections)
{
	UnitCell *cell;
	RefList *reflist;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;

	fprintf(st->fh, CRYSTAL_START_MARKER"\n");

	cell = crystal_get_cell(cr);
	assert(cell != NULL);

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	fprintf(st->fh, "Cell parameters %7.5f %7.5f %7.5f nm,"
		        " %7.5f %7.5f %7.5f deg\n",
		        a*1.0e9, b*1.0e9, c*1.0e9,
		        rad2deg(al), rad2deg(be), rad2deg(ga));

	cell_get_reciprocal(cell, &asx, &asy, &asz,
		                  &bsx, &bsy, &bsz,
		                  &csx, &csy, &csz);
	fprintf(st->fh, "astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
		asx/1e9, asy/1e9, asz/1e9);
	fprintf(st->fh, "bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
		bsx/1e9, bsy/1e9, bsz/1e9);
	fprintf(st->fh, "cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	        csx/1e9, csy/1e9, csz/1e9);

	reflist = crystal_get_reflections(cr);
	if ( reflist != NULL ) {

		fprintf(st->fh, "diffraction_resolution_limit"
		                " = %.2f nm^-1 or %.2f A\n",
		                crystal_get_resolution_limit(cr)/1e9,
		                1e9 / crystal_get_resolution_limit(cr));

		fprintf(st->fh, "num_saturated_reflections = %lli\n",
		                crystal_get_num_saturated_reflections(cr));

	}

	if ( include_reflections ) {

		if ( reflist != NULL ) {

			fprintf(st->fh, REFLECTION_START_MARKER"\n");
			write_reflections_to_file(st->fh, reflist);
			fprintf(st->fh, REFLECTION_END_MARKER"\n");

		} else {

			fprintf(st->fh, "No integrated reflections.\n");

		}
	}

	fprintf(st->fh, CRYSTAL_END_MARKER"\n");
}


void write_chunk(Stream *st, struct image *i, struct hdfile *hdfile,
                 int include_peaks, int include_reflections)
{
	int j;

	fprintf(st->fh, CHUNK_START_MARKER"\n");

	fprintf(st->fh, "Image filename: %s\n", i->filename);
	fprintf(st->fh, "indexed_by = %s\n", indexer_str(i->indexed_by));

	if ( i->det != NULL ) {

		int j;

		for ( j=0; j<i->det->n_panels; j++ ) {
			fprintf(st->fh, "camera_length_%s = %f\n",
			        i->det->panels[j].name, i->det->panels[j].clen);
		}

	}

	copy_hdf5_fields(hdfile, i->copyme, st->fh);

	fprintf(st->fh, "num_peaks = %lli\n", i->num_peaks);
	fprintf(st->fh, "num_saturated_peaks = %lli\n", i->num_saturated_peaks);
	if ( include_peaks ) {
		write_peaks(i, st->fh);
	}

	fprintf(st->fh, "photon_energy_eV = %f\n",
	        J_to_eV(ph_lambda_to_en(i->lambda)));

	for ( j=0; j<i->n_crystals; j++ ) {
		write_crystal(st, i->crystals[j], include_reflections);
	}

	fprintf(st->fh, CHUNK_END_MARKER"\n");

	fflush(st->fh);
}


static int find_start_of_chunk(FILE *fh)
{
	char *rval = NULL;
	char line[1024];

	do {

		rval = fgets(line, 1023, fh);

		/* Trouble? */
		if ( rval == NULL ) return 1;

		chomp(line);

	} while ( strcmp(line, CHUNK_START_MARKER) != 0 );

	return 0;
}


void read_crystal(Stream *st, struct image *image)
{
	char line[1024];
	char *rval = NULL;
	struct rvec as, bs, cs;
	int have_as = 0;
	int have_bs = 0;
	int have_cs = 0;
	Crystal *cr;
	int n;
	Crystal **crystals_new;

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal!\n");
		return;
	}

	do {

		float u, v, w;

		rval = fgets(line, 1023, st->fh);

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);
		if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) == 3 ) {
			as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
			have_as = 1;
		}

		if ( sscanf(line, "bstar = %f %f %f", &u, &v, &w) == 3 ) {
			bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
			have_bs = 1;
		}

		if ( sscanf(line, "cstar = %f %f %f", &u, &v, &w) == 3 ) {
			cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
			have_cs = 1;
		}

		if ( have_as && have_bs && have_cs ) {

			UnitCell *cell;

			cell = crystal_get_cell(cr);

			if ( cell != NULL ) {
				ERROR("Duplicate cell found in stream!\n");
				ERROR("I'll use the most recent one.\n");
				cell_free(cell);
			}

			cell = cell_new_from_reciprocal_axes(as, bs, cs);
			crystal_set_cell(cr, cell);

			have_as = 0;  have_bs = 0;  have_cs = 0;

		}

		if ( strncmp(line, "num_saturated_reflections = ", 28) == 0 ) {
			int n = atoi(line+28);
			crystal_set_num_saturated_reflections(cr, n);
		}

		if ( strcmp(line, REFLECTION_START_MARKER) == 0 ) {

			RefList *reflections;

			reflections = read_reflections_from_file(st->fh);

			if ( reflections == NULL ) {
				ERROR("Failed while reading reflections\n");
				break;
			}

			crystal_set_reflections(cr, reflections);

		}

		if ( strcmp(line, CRYSTAL_END_MARKER) == 0 ) break;

	} while ( 1 );

	/* Add crystal to the list for this image */
	n = image->n_crystals+1;
	crystals_new = realloc(image->crystals, n*sizeof(Crystal *));

	if ( crystals_new == NULL ) {
		ERROR("Failed to expand crystal list!\n");
	} else {
		image->crystals = crystals_new;
		image->crystals[image->n_crystals++] = cr;
	}

}


/* Read the next chunk from a stream and fill in 'image' */
int read_chunk(Stream *st, struct image *image)
{
	char line[1024];
	char *rval = NULL;
	int have_filename = 0;
	int have_ev = 0;

	if ( find_start_of_chunk(st->fh) ) return 1;

	image->lambda = -1.0;
	image->features = NULL;
	image->crystals = NULL;
	image->n_crystals = 0;

	do {

		rval = fgets(line, 1023, st->fh);

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);

		if ( strncmp(line, "Image filename: ", 16) == 0 ) {
			image->filename = strdup(line+16);
			have_filename = 1;
		}

		if ( strncmp(line, "photon_energy_eV = ", 19) == 0 ) {
			image->lambda = ph_en_to_lambda(eV_to_J(atof(line+19)));
			have_ev = 1;
		}

		if ( strncmp(line, "camera_length_", 14) == 0 ) {
			if ( image->det != NULL ) {

				int k;
				char name[1024];
				struct panel *p;

				for ( k=0; k<strlen(line)-14; k++ ) {
					char ch = line[k+14];
					name[k] = ch;
					if ( (ch == ' ') || (ch == '=') ) {
						name[k] = '\0';
						break;
					}
				}

				p = find_panel_by_name(image->det, name);
				if ( p == NULL ) {
					ERROR("No panel '%s'\n", name);
				} else {
					p->clen = atof(line+14+k+3);
				}

			}
		}

		if ( strcmp(line, PEAK_LIST_START_MARKER) == 0 ) {
			if ( read_peaks(st->fh, image) ) {
				ERROR("Failed while reading peaks\n");
				return 1;
			}
		}

		if ( strcmp(line, CRYSTAL_START_MARKER) == 0 ) {
			read_crystal(st, image);
		}

		/* A chunk must have at least a filename and a wavelength,
		 * otherwise it's incomplete */
		if ( strcmp(line, CHUNK_END_MARKER) == 0 ) {
			if ( have_filename && have_ev ) return 0;
			ERROR("Incomplete chunk found in input file.\n");
			return 1;
		}

	} while ( 1 );

	if ( !feof(st->fh) ) {
		ERROR("Error reading stream.\n");
	}

	return 1;  /* Either error or EOF, don't care because we will complain
	            * on the terminal if it was an error. */
}


void write_stream_header(FILE *ofh, int argc, char *argv[])
{
	int i;

	fprintf(ofh, "Command line:");
	for ( i=0; i<argc; i++ ) {
		fprintf(ofh, " %s", argv[i]);
	}
	fprintf(ofh, "\n");
	fflush(ofh);
}


Stream *open_stream_for_read(const char *filename)
{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;

	st->fh = fopen(filename, "r");
	if ( st->fh == NULL ) {
		free(st);
		return NULL;
	}

	char line[1024];
	char *rval;

	rval = fgets(line, 1023, st->fh);
	if ( rval == NULL ) {
		ERROR("Failed to read stream version.\n");
		close_stream(st);
		return NULL;
	}

	if ( strncmp(line, "CrystFEL stream format 2.0", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 0;
	} else if ( strncmp(line, "CrystFEL stream format 2.1", 26) == 0 ) {
		st->major_version = 2;
		st->minor_version = 1;
	} else {
		ERROR("Invalid stream, or stream format is too new.\n");
		close_stream(st);
		return NULL;
	}

	return st;
}


Stream *open_stream_fd_for_write(int fd)
{
	Stream *st;

	st = malloc(sizeof(struct _stream));
	if ( st == NULL ) return NULL;

	st->fh = fdopen(fd, "w");
	if ( st->fh == NULL ) {
		free(st);
		return NULL;
	}

	st->major_version = LATEST_MAJOR_VERSION;
	st->minor_version = LATEST_MINOR_VERSION;

	fprintf(st->fh, "CrystFEL stream format %i.%i\n",
	        st->major_version, st->minor_version);

	return st;
}



Stream *open_stream_for_write(const char *filename)
{
	int fd;

	fd = open(filename, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR);
	if ( fd == -1 ) {
		ERROR("Failed to open stream.\n");
		return NULL;
	}

	return open_stream_fd_for_write(fd);
}


void close_stream(Stream *st)
{
	fclose(st->fh);
	free(st);
}


int is_stream(const char *filename)
{
	FILE *fh;
	char line[1024];
	char *rval;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return 0;

	rval = fgets(line, 1023, fh);
	fclose(fh);
	if ( rval == NULL ) return 0;

	if ( strncmp(line, "CrystFEL stream format 2.0", 26) == 0 ) return 1;
	if ( strncmp(line, "CrystFEL stream format 2.1", 26) == 0 ) return 1;

	return 0;
}


void write_line(Stream *st, const char *line)
{
	fprintf(st->fh, "%s\n", line);
}


void write_command(Stream *st, int argc, char *argv[])
{
	int i;

	for ( i=0; i<argc; i++ ) {
		if ( i > 0 ) fprintf(st->fh, " ");
		fprintf(st->fh, "%s", argv[i]);
	}
	fprintf(st->fh, "\n");
}


/**
 * rewind_stream:
 * @st: A %Stream
 *
 * Attempts to set the file pointer for @st to the start of the stream, so that
 * later calls to read_chunk() will repeat the sequence of chunks from the
 * start.
 *
 * Programs must not assume that this operation always succeeds!
 *
 * Returns: non-zero if the stream could not be rewound.
 */
int rewind_stream(Stream *st)
{
	return fseek(st->fh, 0, SEEK_SET);
}
