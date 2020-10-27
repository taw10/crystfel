/*
 * wavelength_geom.c
 *
 * Check that wavelength reading works
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
#include <string.h>
#include <stdlib.h>

#include <image.h>

int main(int argc, char *argv[])
{
	DataTemplate *dtempl;
	struct image *image;
	const char *geom_filename;
	const char *image_filename;
	double expected_wavelength_m;

	image_filename = argv[1];
	geom_filename = argv[2];
	expected_wavelength_m = atof(argv[3]);

	dtempl = data_template_new_from_file(geom_filename);
	if ( dtempl == NULL ) {
		ERROR("Failed to load data template\n");
		return 1;
	}

	image = image_read(dtempl, image_filename, NULL, 0, 0);
	if ( image == NULL ) return 1;

	printf("wavelength = %e\n", image->lambda);
	printf("should be %e\n", expected_wavelength_m);

	if ( !within_tolerance(image->lambda,
	                       expected_wavelength_m,
	                       0.1) ) return 1;

	data_template_free(dtempl);

	return 0;
}
