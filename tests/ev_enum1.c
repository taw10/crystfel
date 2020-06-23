/*
 * ev_enum1.c
 *
 * Check that event enumeration works
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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
#include <stdarg.h>
#include <hdf5.h>

#include "../libcrystfel/src/image-hdf5.c"

int main(int argc, char *argv[])
{
	hid_t fh;
	char **event_ids;
	int n_event_ids;
	int i;

	fh = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fh < 0 ) {
		ERROR("Couldn't open file\n");
		return 1;
	}

	event_ids = expand_paths(fh,
	                         "/data/panelA/%/panel_data1t/%/array",
	                         &n_event_ids);

	if ( event_ids == NULL ) {
		STATUS("event_ids = NULL\n");
		return 1;
	}

	if ( n_event_ids != 4 ) {
		STATUS("Number of event IDs = %i\n", n_event_ids);
		return 1;
	}

	if ( strcmp(event_ids[0], "/ev_1/dataABCset") != 0 ) {
		STATUS("Wrong event id '%s'\n", event_ids[0]);
		return 1;
	}

	if ( strcmp(event_ids[1], "/ev_2/dataDEFset") != 0 ) {
		STATUS("Wrong event id '%s'\n", event_ids[1]);
		return 1;
	}

	if ( strcmp(event_ids[2], "/ev_3/dataGHIset") != 0 ) {
		STATUS("Wrong event id '%s'\n", event_ids[2]);
		return 1;
	}

	if ( strcmp(event_ids[3], "/ev_5/dataNOPset") != 0 ) {
		STATUS("Wrong event id '%s'\n", event_ids[3]);
		return 1;
	}

	for ( i=0; i<n_event_ids; i++ ) {
		free(event_ids[i]);
	}
	free(event_ids);

	H5Fclose(fh);

	return 0;
}
