/*
 * ev_enum2.c
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
#include <stdlib.h>
#include <string.h>

#include <image.h>

int main(int argc, char *argv[])
{
	char **event_ids;
	int n_event_ids;
	int i;
	DataTemplate *dtempl;

	dtempl = data_template_new_from_file(argv[2]);
	if ( dtempl == NULL ) {
		ERROR("Failed ot load data template\n");
		return 1;
	}

	event_ids = image_expand_frames(dtempl, argv[1], &n_event_ids);

	if ( event_ids == NULL ) {
		printf("event_ids = NULL\n");
		return 1;
	}

	for ( i=0; i<n_event_ids; i++ ) {
		char tmp[64];
		char c = i < 100 ? 'a' : 'b';
		int n = i < 100 ? i : (i-100);
		snprintf(tmp, 64, "%c//%i", c, n);
		if ( strcmp(tmp, event_ids[i]) != 0 ) {
			printf("Event ID %i is wrong '%s'\n",
			       i, event_ids[i]);
			return 1;
		}
		free(event_ids[i]);
	}
	free(event_ids);

	data_template_free(dtempl);

	return 0;
}
