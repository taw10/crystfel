/*
 * stream.c
 *
 * Stream tools
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 * (c) 2011 Rick Kirian <rkirian@asu.edu>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cell.h"
#include "utils.h"
#include "image.h"
#include "stream.h"
#include "reflist.h"
#include "reflist-utils.h"


#define CHUNK_START_MARKER "----- Begin chunk -----"
#define CHUNK_END_MARKER "----- End chunk -----"
#define PEAK_LIST_START_MARKER "Peaks from peak search"
#define PEAK_LIST_END_MARKER "End of peak list"
#define REFLECTION_START_MARKER "Reflections measured after indexing"
/* REFLECTION_END_MARKER is over in reflist-utils.h because it is also
 * used to terminate a standalone list of reflections */

static void exclusive(const char *a, const char *b)
{
	ERROR("The stream options '%s' and '%s' are mutually exclusive.\n",
	      a, b);
}


int parse_stream_flags(const char *a)
{
	int n, i;
	int ret = STREAM_NONE;
	char **flags;

	n = assplode(a, ",", &flags, ASSPLODE_NONE);

	for ( i=0; i<n; i++ ) {

		if ( strcmp(flags[i], "pixels") == 0) {
			if ( ret & STREAM_INTEGRATED ) {
				exclusive("pixels", "integrated");
				return -1;
			}
			ret |= STREAM_PIXELS;

		} else if ( strcmp(flags[i], "integrated") == 0) {
			if ( ret & STREAM_PIXELS ) {
				exclusive("pixels", "integrated");
				return -1;
			}
			ret |= STREAM_INTEGRATED;

		} else if ( strcmp(flags[i], "peaks") == 0) {
			if ( ret & STREAM_PEAKS_IF_INDEXED ) {
				exclusive("peaks", "peaksifindexed");
				return -1;
			}
			if ( ret & STREAM_PEAKS_IF_NOT_INDEXED ) {
				exclusive("peaks", "peaksifnotindexed");
				return -1;
			}
			ret |= STREAM_PEAKS;

		} else if ( strcmp(flags[i], "peaksifindexed") == 0) {
			if ( ret & STREAM_PEAKS ) {
				exclusive("peaks", "peaksifindexed");
				return -1;
			}
			if ( ret & STREAM_PEAKS_IF_NOT_INDEXED ) {
				exclusive("peaksifnotindexed",
				          "peaksifindexed");
				return -1;
			}
			ret |= STREAM_PEAKS_IF_INDEXED;

		} else if ( strcmp(flags[i], "peaksifnotindexed") == 0) {
			if ( ret & STREAM_PEAKS ) {
				exclusive("peaks", "peaksifnotindexed");
				return -1;
			}
			if ( ret & STREAM_PEAKS_IF_INDEXED ) {
				exclusive("peaksifnotindexed",
				          "peaksifindexed");
				return -1;
			}
			ret |= STREAM_PEAKS_IF_NOT_INDEXED;

		} else {
			ERROR("Unrecognised stream flag '%s'\n", flags[i]);
			return -1;
		}

		free(flags[i]);

	}
	free(flags);

	return ret;
}


int count_patterns(FILE *fh)
{
	char *rval;

	int n_total_patterns = 0;
	do {
		char line[1024];

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);
		if ( strcmp(line, CHUNK_END_MARKER) == 0 ) n_total_patterns++;

	} while ( rval != NULL );

	if ( ferror(fh) ) {
		ERROR("Read error while counting patterns.\n");
		return 0;
	}

	return n_total_patterns;
}


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


void write_chunk(FILE *ofh, struct image *i, int f)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;

	fprintf(ofh, CHUNK_START_MARKER"\n");

	fprintf(ofh, "Image filename: %s\n", i->filename);

	if ( i->indexed_cell != NULL ) {

		cell_get_parameters(i->indexed_cell, &a, &b, &c,
		                                         &al, &be, &ga);
		fprintf(ofh, "Cell parameters %7.5f %7.5f %7.5f nm,"
			     " %7.5f %7.5f %7.5f deg\n",
			     a*1.0e9, b*1.0e9, c*1.0e9,
			     rad2deg(al), rad2deg(be), rad2deg(ga));

		cell_get_reciprocal(i->indexed_cell, &asx, &asy, &asz,
			                  &bsx, &bsy, &bsz,
			                  &csx, &csy, &csz);
		fprintf(ofh, "astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
			asx/1e9, asy/1e9, asz/1e9);
		fprintf(ofh, "bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
			bsx/1e9, bsy/1e9, bsz/1e9);
		fprintf(ofh, "cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
		        csx/1e9, csy/1e9, csz/1e9);

	} else {

		fprintf(ofh, "No unit cell from indexing.\n");

	}

	if ( i->i0_available ) {
		fprintf(ofh, "I0 = %7.5f (arbitrary units)\n", i->i0);
	} else {
		fprintf(ofh, "I0 = invalid\n");
	}

	fprintf(ofh, "photon_energy_eV = %f\n",
	        J_to_eV(ph_lambda_to_en(i->lambda)));

	if ( i->det != NULL ) {

		int j;

		for ( j=0; j<i->det->n_panels; j++ ) {
			fprintf(ofh, "camera_length_%s = %f\n",
			        i->det->panels[j].name, i->det->panels[j].clen);
		}

	}

	if ( (f & STREAM_PEAKS)
	  || ((f & STREAM_PEAKS_IF_INDEXED) && (i->indexed_cell != NULL))
	  || ((f & STREAM_PEAKS_IF_NOT_INDEXED) && (i->indexed_cell == NULL)) )
	{
		fprintf(ofh, "\n");
		write_peaks(i, ofh);
	}

	if ( (f & STREAM_PIXELS) || (f & STREAM_INTEGRATED) ) {

		if ( i->reflections != NULL ) {

			fprintf(ofh, "\n");
			fprintf(ofh, REFLECTION_START_MARKER"\n");
			write_reflections_to_file(ofh, i->reflections,
			                               i->indexed_cell);
			fprintf(ofh, REFLECTION_END_MARKER"\n");

		}
	}

	fprintf(ofh, CHUNK_END_MARKER"\n\n");
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


/* Read the next chunk from a stream and fill in 'image' */
int read_chunk(FILE *fh, struct image *image)
{
	char line[1024];
	char *rval = NULL;
	struct rvec as, bs, cs;
	int have_as = 0;
	int have_bs = 0;
	int have_cs = 0;
	int have_filename = 0;
	int have_cell = 0;
	int have_ev = 0;

	if ( find_start_of_chunk(fh) ) return 1;

	image->i0_available = 0;
	image->i0 = 1.0;
	image->lambda = -1.0;
	image->features = NULL;
	image->reflections = NULL;
	image->indexed_cell = NULL;

	do {

		float u, v, w;

		rval = fgets(line, 1023, fh);

		/* Trouble? */
		if ( rval == NULL ) break;

		chomp(line);

		if ( strncmp(line, "Image filename: ", 16) == 0 ) {
			image->filename = strdup(line+16);
			have_filename = 1;
		}

		if ( strncmp(line, "camera_length_", 14) == 0 ) {
			if ( image->det == NULL ) {
				ERROR("Stream had a camera length, but "
				      "geometry is not currently loaded.\n");
			} else {

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

		if ( strncmp(line, "I0 = ", 5) == 0 ) {
			image->i0 = atof(line+5);
			image->i0_available = 1;
		}

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
			if ( image->indexed_cell != NULL ) {
				ERROR("Duplicate cell found in stream!\n");
				cell_free(image->indexed_cell);
			}
			image->indexed_cell = cell_new_from_reciprocal_axes(as,
			                                                    bs,
			                                                    cs);
			have_cell = 1;
			have_as = 0;  have_bs = 0;  have_cs = 0;
		}

		if ( strncmp(line, "photon_energy_eV = ", 19) == 0 ) {
			image->lambda = ph_en_to_lambda(eV_to_J(atof(line+19)));
			have_ev = 1;
		}

		if ( strcmp(line, PEAK_LIST_START_MARKER) == 0 ) {
			if ( read_peaks(fh, image) ) {
				ERROR("Failed while reading peaks\n");
				return 1;
			}
		}

		if ( strcmp(line, REFLECTION_START_MARKER) == 0 ) {
			image->reflections = read_reflections_from_file(fh);
			if ( image->reflections == NULL ) {
				ERROR("Failed while reading reflections\n");
				return 1;
			}
		}

		if ( strcmp(line, CHUNK_END_MARKER) == 0 ) {
			if ( have_filename && have_ev ) return 0;
			ERROR("Incomplete chunk found in input file.\n");
			return 1;
		}

	} while ( 1 );

	return 1;  /* Either error or EOF, don't care because we will complain
	            * on the terminal if it was an error. */
}


void write_stream_header(FILE *ofh, int argc, char *argv[])
{
	int i;

	fprintf(ofh, "CrystFEL stream format 2.0\n");
	fprintf(ofh, "Command line:");
	for ( i=0; i<argc; i++ ) {
		fprintf(ofh, " %s", argv[i]);
	}
	fprintf(ofh, "\n");
	fflush(ofh);
}


int skip_some_files(FILE *fh, int n)
{
	char *rval = NULL;
	int n_patterns = 0;

	do {

		char line[1024];

		if ( n_patterns == n ) return 0;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		if ( strcmp(line, CHUNK_END_MARKER) == 0 ) n_patterns++;

	} while ( rval != NULL );

	return 1;
}

int is_stream(const char *filename) {
	FILE *fh;
	char line[1024];
	char *rval = NULL;
	fh = fopen(filename, "r");
	rval = fgets(line, 1023, fh);
	fclose(fh);
	if ( rval == NULL ) {
		return -1;
	}
	if ( strncmp(line, "CrystFEL stream format 2.0", 26) == 0 ) {
		return 1;
	}
	else {
		return 0;
	}
	return -1;
}
