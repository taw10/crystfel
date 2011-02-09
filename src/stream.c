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


static int find_cell(FILE *fh)
{
	int done = 0;
	int found = 0;

	do {

		long pos;
		char *rval;
		float u, v, w;
		char line[1024];

		pos = ftell(fh);

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) {
			STATUS("Read error in find_cell()\n");
			done = 1;
		}

		chomp(line);

		if ( strncmp(line, "Reflections from indexing", 25) == 0 ) {
			done = 1;
		}

		if ( sscanf(line, "astar = %f %f %f", &u, &v, &w) == 3 ) {
			fseek(fh, pos, SEEK_SET);
			done = 1;
			found = 1;
		}

	} while ( !done );

	return found;
}


static UnitCell *read_orientation_matrix(FILE *fh)
{
	float u, v, w;
	struct rvec as, bs, cs;
	UnitCell *cell;
	char line[1024];

	if ( find_cell(fh) == 0 ) {
		ERROR("Couldn't find orientation matrix.\n");
		return NULL;
	}

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


static UnitCell *read_orientation_matrix_rick(FILE *fh)
{
	float a, b, c;
	struct rvec as, bs, cs;
	UnitCell *cell;
	char line[1024];

	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "%f %f %f", &a, &b, &c) != 3 ) {
		ERROR("Couldn't read a-star\n");
		return NULL;
	}
	as.u = a*1e10;  bs.u = b*1e10;  cs.u = c*1e10;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "%f %f %f", &a, &b, &c) != 3 ) {
		ERROR("Couldn't read b-star\n");
		return NULL;
	}
	as.v = a*1e10;  bs.v = b*1e10;  cs.v = c*1e10;
	if ( fgets(line, 1023, fh) == NULL ) return NULL;
	if ( sscanf(line, "%f %f %f", &a, &b, &c) != 3 ) {
		ERROR("Couldn't read c-star\n");
		return NULL;
	}
	as.w = -a*1e10;  bs.w = -b*1e10;  cs.w = -c*1e10;
	cell = cell_new_from_axes(as, bs, cs);

	return cell;
}


int find_chunk(FILE *fh, UnitCell **cell, char **filename, double *ev)
{
	char line[1024];
	char *rval = NULL;

	do {

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;

		chomp(line);

		if ( strncmp(line, "photon_energy_eV = ", 19) == 0 ) {
			*ev = atof(line+19);
		}

		/* Look for the first line of a chunk */
		if ( (strncmp(line, "Reflections from indexing", 25) != 0)
		  && (strncmp(line, "## h5FilePath:", 14) != 0 ) ) {
			continue;
		}

		/* Read in "Tom Mode"? */
		if ( strncmp(line, "Reflections from indexing", 25) == 0 ) {

			*filename = strdup(line+29);
			*cell = read_orientation_matrix(fh);

		}

		/* Read in "Rick Mode"? */
		if ( strncmp(line, "## h5FilePath:", 14) == 0 ) {

			/* Filename is on next line */
			rval = fgets(line, 1023, fh);
			if ( rval == NULL ) continue;
			chomp(line);
			*filename = strdup(line);
			/* Look for the start of the orientation matrix */
			do {
				rval = fgets(line, 1023, fh);
				if ( rval == NULL ) continue;
			} while ( strncmp(line, "## A:", 5) != 0 );
			*cell = read_orientation_matrix_rick(fh);
		}

		if ( *cell == NULL ) {
			STATUS("Got filename but no cell for %s\n", *filename);
			continue;
		}

		return 0;

	} while ( rval != NULL );

	return 1;
}
