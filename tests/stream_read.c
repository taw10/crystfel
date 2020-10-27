/*
 * stream_read.c
 *
 * Simple test of stream reading
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

#include "stream.h"
#include "image.h"

int main(int argc, char *argv[])
{
	Stream *st;
	int n;
	char *stream_filename = argv[1];

	st = stream_open_for_read(stream_filename);
	if ( st == NULL ) {
		fprintf(stderr, "Failed to open '%s'\n",
		        stream_filename);
		return 1;
	}

	n = 0;
	do {

		struct image *image = stream_read_chunk(st, 0);
		if ( image == NULL ) break;
		n++;

		image_free(image);

	} while ( 1 );

	printf("Got %i chunks\n", n);

	return (n != 70);
}
