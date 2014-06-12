/*
 * whirligig.c
 *
 * Find and combine rotation series
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2014 Thomas White <taw@physics.org>
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <pthread.h>
#include <gsl/gsl_errno.h>

#include <image.h>
#include <utils.h>
#include <symmetry.h>
#include <stream.h>
#include <geometry.h>
#include <peaks.h>
#include <thread-pool.h>
#include <reflist.h>
#include <reflist-utils.h>

#include "version.h"


static void show_help(const char *s)
{
	printf("Syntax: %s <input.stream> [options]\n\n", s);
	printf(
"Find and combine rotation series.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"      --no-polarisation      Disable polarisation correction.\n");
}


int main(int argc, char *argv[])
{
	int c;
	Stream *st;
	int n_images = 0;
	int n_crystals = 0;
	int polarisation = 1;

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"version",            0, NULL,                3 },

		{"no-polarisation",    0, &polarisation,       0},
		{"no-polarization",    0, &polarisation,       0},
		{"polarisation",       0, &polarisation,       1}, /* compat */
		{"polarization",       0, &polarisation,       1}, /* compat */

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "h",
	                        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	if ( argc != (optind+1) ) {
		ERROR("Please provide exactly one stream to process.\n");
		return 1;
	}

	st = open_stream_for_read(argv[optind++]);
	if ( st == NULL ) {
		ERROR("Failed to open input stream '%s'\n", argv[optind-1]);
		return 1;
	}

	do {

		int i;
		struct image cur;

		cur.div = NAN;
		cur.bw = NAN;
		cur.det = NULL;
		if ( read_chunk_2(st, &cur, STREAM_READ_REFLECTIONS
		                            | STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		if ( isnan(cur.div) || isnan(cur.bw) ) {
			ERROR("Chunk doesn't contain beam parameters.\n");
			return 1;
		}

		n_images++;

		for ( i=0; i<cur.n_crystals; i++ ) {

			Crystal *cr;
			RefList *cr_refl;

			cr = cur.crystals[i];

			/* This is the raw list of reflections */
			cr_refl = crystal_get_reflections(cr);

			if ( polarisation ) {
				polarisation_correction(cr_refl,
						        crystal_get_cell(cr),
						        &cur);
			}

			n_crystals++;

		}

	} while ( 1 );

	close_stream(st);

	return 0;
}
