/*
 * beam-parameters.c
 *
 * Beam parameters
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010,2012,2014 Thomas White <taw@physics.org>
 *   2012           Chunhong Yoon
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

#include <stdio.h>
#include <stdlib.h>

#include "beam-parameters.h"
#include "utils.h"
#include "hdf5-file.h"


struct beam_params *get_beam_parameters(const char *filename)
{
	FILE *fh;
	struct beam_params *b;
	char *rval;
	int reject;

	fh = fopen(filename, "r");
	if ( fh == NULL ) return NULL;

	b = calloc(1, sizeof(struct beam_params));
	if ( b == NULL ) return NULL;

	b->fluence = -1.0;
	b->beam_radius = -1.0;
	b->photon_energy = -1.0;
	b->photon_energy_scale = 1.0;
	b->bandwidth = -1.0;
	b->divergence = -1.0;
	b->profile_radius = -1.0;

	do {

		int n1;
		char line[1024];
		int i;
		char **bits;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) break;
		chomp(line);

		n1 = assplode(line, " \t", &bits, ASSPLODE_NONE);
		if ( n1 < 3 ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		if ( bits[0][0] == ';' ) {
			for ( i=0; i<n1; i++ ) free(bits[i]);
			free(bits);
			continue;
		}

		if ( strcmp(bits[0], "beam/fluence") == 0 ) {
			b->fluence = atof(bits[2]);
		} else if ( strcmp(bits[0], "beam/radius") == 0 ) {
			b->beam_radius = atof(bits[2]);
		} else if ( strcmp(bits[0], "beam/photon_energy") == 0 ) {
			if ( strncmp(bits[2], "/", 1) == 0 ) {
				b->photon_energy = 0; // 0 means special case
				b->photon_energy_from = strdup(bits[2]);
			} else {
				b->photon_energy = atof(bits[2]);
			}
		} else if ( strcmp(bits[0], "beam/photon_energy_scale") == 0 ) {
			b->photon_energy_scale = atof(bits[2]);
		} else if ( strcmp(bits[0], "beam/bandwidth") == 0 ) {
			b->bandwidth = atof(bits[2]);
		} else if ( strcmp(bits[0], "beam/divergence") == 0 ) {
			b->divergence = atof(bits[2]);
		} else if ( strcmp(bits[0], "profile_radius") == 0 ) {
			b->profile_radius = atof(bits[2]);
		} else {
			ERROR("Unrecognised field '%s'\n", bits[0]);
		}

		for ( i=0; i<n1; i++ ) free(bits[i]);
		free(bits);

	} while ( rval != NULL );
	fclose(fh);

	reject = 0;

	if ( b->fluence < 0.0 ) {
		ERROR("Invalid or unspecified value for 'beam/fluence'.\n");
		reject = 1;
	}
	if ( b->beam_radius < 0.0 ) {
		ERROR("Invalid or unspecified value for 'beam/radius'.\n");
		reject = 1;
	}
	if ( b->photon_energy < 0.0 ) { // 0 is ok
		ERROR("Invalid or unspecified value for"
		      " 'beam/photon_energy'.\n");
		reject = 1;
	}
	if ( b->bandwidth < 0.0 ) {
		ERROR("Invalid or unspecified value for 'beam/bandwidth'.\n");
		reject = 1;
	}
	if ( b->divergence < 0.0 ) {
		ERROR("Invalid or unspecified value for 'beam/divergence'.\n");
		reject = 1;
	}
	if ( b->profile_radius < 0.0 ) {
		ERROR("WARNING: Invalid or unspecified value for "
		      "'profile_radius' in the beam parameters file.\n");
		ERROR("Using the default value of 0.001 nm^-1.\n");
		b->profile_radius = 0.001e9;
	}

	if ( reject ) {
		ERROR("Please fix the above problems with the beam"
		      " parameters file and try again.\n");
		return NULL;
	}

	return b;
}


void free_beam_parameters(struct beam_params *beam)
{
	free(beam->photon_energy_from);
	free(beam);
}


void fill_in_beam_parameters(struct beam_params *beam, struct hdfile *f)
{
	if ( beam->photon_energy_from != NULL ) {
		beam->photon_energy = get_value(f, beam->photon_energy_from);
		beam->photon_energy *= beam->photon_energy_scale;
	}
}
