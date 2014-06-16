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

#include <image.h>
#include <utils.h>
#include <stream.h>

#include "version.h"


static void process_series(struct image *images, signed int *ser, int len)
{
	STATUS("Found a rotation series of %i views\n", len);
}


static int gatinator(UnitCell *a, UnitCell *b)
{
	return 0;
}


static int try_all(struct image *a, struct image *b, int *c1, int *c2)
{
	int i, j;

	for ( i=0; i<a->n_crystals; i++ ) {
		for ( j=0; j<b->n_crystals; j++ ) {
			if ( gatinator(crystal_get_cell(a->crystals[i]),
			               crystal_get_cell(b->crystals[j])) )
			{
				*c1 = i;
				*c2 = j;
				return 1;
			}
		}
	}

	return 0;
}


static void dump(struct image *win, signed int *series, int window_len, int pos)
{
	int i;

	for ( i=0; i<pos; i++ ) {
		free_all_crystals(&win[i]);
	}

	memmove(win, &win[pos], (window_len-pos)*sizeof(struct image *));
	memmove(series, &series[pos], (window_len-pos)*sizeof(signed int));
}


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
	int polarisation = 1;
	int pos = 0;
	struct image *win;
	signed int *series;
	int window_len = 64;

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

	win = calloc(window_len, sizeof(struct image));
	series = calloc(window_len, sizeof(int));
	if ( (win == NULL) || (series == NULL) ) {
		ERROR("Failed to allocate series buffers\n");
		return 1;
	}

	pos = 0;
	do {

		struct image *cur;
		int c1, c2;

		cur = &win[pos];

		cur->div = NAN;
		cur->bw = NAN;
		cur->det = NULL;
		if ( read_chunk_2(st, cur, STREAM_READ_REFLECTIONS
		                            | STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		if ( isnan(cur->div) || isnan(cur->bw) ) {
			ERROR("Chunk doesn't contain beam parameters.\n");
			return 1;
		}

		/* Need at least two images to compare */
		if ( pos == 0 ) {
			series[0] = -1;
			pos++;
			continue;
		}

		if ( try_all(&win[pos-1], cur, &c1, &c2) ) {
			series[pos-1] = c1;
			series[pos] = c2;
			printf("-");
		} else {
			series[pos] = -1;
			printf(".");
		}
		fflush(stdout);

		if ( series[0] == -1 ) {
			dump(win, series, window_len, 1);
			pos--;
		}

		if ( (series[pos] == -1) && (series[pos-1] == -1) ) {
			/* Series ready to process */
			process_series(win, series, pos-2);
			dump(win, series, window_len, pos);
			pos = 0;
		}

		pos++;
		if ( pos == window_len ) {
			window_len *= 2;
			win = realloc(win, window_len*sizeof(struct image));
			series = realloc(series, window_len*sizeof(signed int));
			if ( (win == NULL) || (series == NULL) ) {
				ERROR("Failed to expand series buffers\n");
				return 1;
			}
		}

	} while ( 1 );
	printf("\n");

	close_stream(st);

	free(win);
	free(series);

	return 0;
}
