/*
 * ev_enum3.c
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

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <image.h>

int main(int argc, char *argv[])
{
	char **event_ids;
	int n_event_ids;
	DataTemplate *dtempl;

	dtempl = data_template_new_from_file(argv[2]);
	if ( dtempl == NULL ) {
		ERROR("Failed to load data template\n");
		return 1;
	}

	event_ids = image_expand_frames(dtempl, argv[1], &n_event_ids);

	if ( n_event_ids != 1 ) {
		printf("n_event_ids = %i\n", n_event_ids);
		return 1;
	}

	if ( event_ids == NULL ) {
		printf("event_ids not NULL\n");
		return 1;
	}

	if ( strcmp(event_ids[0], "//") != 0 ) {
		printf("Event is not '//' ('%s')\n", event_ids[0]);
		return 1;
	}

	data_template_free(dtempl);

	return 0;
}
