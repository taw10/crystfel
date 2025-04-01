/*
 * header_read.c
 *
 * Check that HDF5 headers can be read correctly
 *
 * Copyright © 2025 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2025 Thomas White <taw@physics.org>
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
#include <unistd.h>

#include <stream.h>
#include <image.h>
#include <datatemplate.h>

int main(int argc, char *argv[])
{
	DataTemplate *dtempl;
	struct image *image;
	Stream *st;
	long long int v;

	dtempl = data_template_new_from_file(argv[1]);
	if ( dtempl == NULL ) {
		ERROR("Failed to load data template\n");
		return 1;
	}

	data_template_add_copy_header(dtempl, argv[4]);

	image = image_read(dtempl, argv[2], argv[3], 0, 0, NULL);
	if ( image == NULL ) {
		ERROR("Failed to load image\n");
		return 1;
	}

	st = stream_open_for_write(argv[5], dtempl);
	if ( st == NULL ) {
		ERROR("Failed to open stream for writing\n");
		return 1;
	}

	stream_write_geometry_file(st, argv[1]);

	if ( stream_write_chunk(st, image, 0) ) {
		ERROR("Failed to write stream chunk\n");
		return 1;
	}

	stream_close(st);
	image_free(image);

	st = stream_open_for_read(argv[5]);
	if ( st == NULL ) {
		ERROR("Failed to open stream for reading\n");
		return 1;
	}

	image = stream_read_chunk(st, 0);
	if ( image == NULL ) {
		ERROR("Failed to read stream chunk\n");
		return 1;
	}
	stream_close(st);

	image_read_header_int(image, argv[4], &v);
	if ( v != 1234567890123456789 ) {
		ERROR("Wrong value read (%lli)\n", v);
		return 1;
	}

	image_free(image);
	unlink("header_read.stream");
	return 0;
}
