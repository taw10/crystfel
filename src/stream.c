/*
 * stream.c
 *
 * Indexed stream tools
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
#include <string.h>

#include "cell.h"
#include "utils.h"
#include "image.h"


int count_patterns(FILE *fh)
{
	char *rval;

	int n_total_patterns = 0;
	do {
		char line[1024];

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		if ( (strncmp(line, "Reflections from indexing", 25) == 0)
		    || (strncmp(line, "New pattern", 11) == 0) ) {
		    n_total_patterns++;
		}
	} while ( rval != NULL );

	return n_total_patterns;
}


static UnitCell *read_orientation_matrix(FILE *fh)
{
	float u, v, w;
	struct rvec as, bs, cs;
	UnitCell *cell;
	char line[1024];

	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read a-star\n");
		return NULL;
	}
	as.u = u*1e9;  as.v = v*1e9;  as.w = w*1e9;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "bstar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read b-star\n");
		return NULL;
	}
	bs.u = u*1e9;  bs.v = v*1e9;  bs.w = w*1e9;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "cstar = %f %f %f", &u, &v, &w) != 3 ) {
		ERROR("Couldn't read c-star\n");
		return NULL;
	}
	cs.u = u*1e9;  cs.v = v*1e9;  cs.w = w*1e9;
	cell = cell_new_from_axes(as, bs, cs);

	return cell;
}


static void write_reflections(struct image *image, FILE *ofh)
{
}


static void write_peaks(struct image *image, FILE *ofh)
{
	int i;

	fprintf(ofh, "Peaks from peak search in %s\n", image->filename);
	fprintf(ofh, "  x/px     y/px   (1/d)/nm^-1    Intensity\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec r;
		double q;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->x, f->y, NULL, 1.0/image->lambda);
		q = modulus(r.u, r.v, r.w);

		fprintf(ofh, "%8.3f %8.3f %8.3f    %12.3f\n",
		       f->x, f->y, q/1.0e9, f->intensity);

	}

	fprintf(ofh, "\n");
}



void write_chunk(FILE *ofh, struct image *image, int flags)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;

	fprintf(ofh, "----- Begin chunk -----\n");

	fprintf(ofh, "Image filename: %s\n", image->filename);

	cell_get_parameters(image->indexed_cell, &a, &b, &c, &al, &be, &ga);
	fprintf(ofh, "Cell parameters %7.5f %7.5f %7.5f nm,"
	             " %7.5f %7.5f %7.5f deg\n",
	             a*1.0e9, b*1.0e9, c*1.0e9,
	             rad2deg(al), rad2deg(be), rad2deg(ga));

	cell_get_reciprocal(image->indexed_cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	fprintf(ofh, "astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	        asx/1e9, asy/1e9, asz/1e9);
	fprintf(ofh, "bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	        bsx/1e9, bsy/1e9, bsz/1e9);
	fprintf(ofh, "cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       csx/1e9, csy/1e9, csz/1e9);

	if ( image->f0_available ) {
		fprintf(ofh, "I0 = %7.5f (arbitrary gas detector units)\n",
		       image->f0);
	} else {
		fprintf(ofh, "I0 = invalid\n");
	}

	fprintf(ofh, "photon_energy_eV = %f\n",
	        J_to_eV(ph_lambda_to_en(image->lambda)));

	write_peaks(image, ofh);
	write_reflections(image, ofh);

	fprintf(ofh, "----- End chunk -----\n\n");
}


int find_chunk(FILE *fh, UnitCell **cell, char **filename, double *ev)
{
	char line[1024];
	char *rval = NULL;
	int have_ev = 0;
	int have_cell = 0;
	int have_filename = 0;
	long start_of_chunk = 0;

	do {

		const long start_of_line = ftell(fh);

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;

		chomp(line);

		if ( strncmp(line, "Reflections from indexing", 25) == 0 ) {

			*filename = strdup(line+29);
			*cell = NULL;
			*ev = 0.0;
			have_cell = 0;
			have_ev = 0;
			have_filename = 1;
			start_of_chunk = ftell(fh);

		}

		if ( !have_filename ) continue;

		if ( strncmp(line, "astar = ", 8) == 0 ) {
			fseek(fh, start_of_line, 0);
			*cell = read_orientation_matrix(fh);
			have_cell = 1;
		}

		if ( strncmp(line, "photon_energy_eV = ", 19) == 0 ) {
			*ev = atof(line+19);
			have_ev = 1;
		}

		if ( strlen(line) == 0 ) {
			if ( have_filename && have_cell && have_ev ) {
				fseek(fh, start_of_chunk, 0);
				return 0;
			}
		}

	} while ( rval != NULL );

	return 1;
}
