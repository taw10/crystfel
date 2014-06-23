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
#include "cell-utils.h"
#include "integer_matrix.h"


static RefList *transform_reflections(RefList *in, IntegerMatrix *m)
{
}


static void process_series(struct image *images, signed int *ser,
                           IntegerMatrix **mat, int len)
{
	int i;
	RefList *list;

	printf("\n");
	STATUS("Found a rotation series of %i views\n", len);

	for ( i=0; i<len; i++ ) {
		Crystal *cr = images[i].crystals[ser[i]];
		RefList *p = transform_reflections(crystal_get_reflections(cr),
		                                   mat[i]);
		reflist_free(p);
	}
}


static double moduli_check(double ax, double ay, double az,
                           double bx, double by, double bz)
{
	double ma = modulus(ax, ay, az);
	double mb = modulus(bx, by, bz);
	return fabs(ma-mb)/ma;
}


static int cells_are_similar(UnitCell *cell1, UnitCell *cell2)
{
	double asx1, asy1, asz1, bsx1, bsy1, bsz1, csx1, csy1, csz1;
	double asx2, asy2, asz2, bsx2, bsy2, bsz2, csx2, csy2, csz2;
	UnitCell *pcell1, *pcell2;
	const double atl = deg2rad(5.0);
	const double ltl = 0.1;

	/* Compare primitive cells, not centered */
	pcell1 = uncenter_cell(cell1, NULL);
	pcell2 = uncenter_cell(cell2, NULL);

	cell_get_reciprocal(pcell1, &asx1, &asy1, &asz1,
	                            &bsx1, &bsy1, &bsz1,
	                            &csx1, &csy1, &csz1);

	cell_get_reciprocal(pcell2, &asx2, &asy2, &asz2,
	                            &bsx2, &bsy2, &bsz2,
	                            &csx2, &csy2, &csz2);


	cell_free(pcell1);
	cell_free(pcell2);

	if ( angle_between(asx1, asy1, asz1, asx2, asy2, asz2) > atl ) return 0;
	if ( angle_between(bsx1, bsy1, bsz1, bsx2, bsy2, bsz2) > atl ) return 0;
	if ( angle_between(csx1, csy1, csz1, csx2, csy2, csz2) > atl ) return 0;

	if ( moduli_check(asx1, asy1, asz1, asx2, asy2, asz2) > ltl ) return 0;
	if ( moduli_check(bsx1, bsy1, bsz1, bsx2, bsy2, bsz2) > ltl ) return 0;
	if ( moduli_check(csx1, csy1, csz1, csx2, csy2, csz2) > ltl ) return 0;

	return 1;
}


static int gatinator(UnitCell *a, UnitCell *b, IntegerMatrix **pmb)
{
	IntegerMatrix *m;
	int i[9];

	m = intmat_new(3, 3);

	for ( i[0]=-1; i[0]<=+1; i[0]++ ) {
	for ( i[1]=-1; i[1]<=+1; i[1]++ ) {
	for ( i[2]=-1; i[2]<=+1; i[2]++ ) {
	for ( i[3]=-1; i[3]<=+1; i[3]++ ) {
	for ( i[4]=-1; i[4]<=+1; i[4]++ ) {
	for ( i[5]=-1; i[5]<=+1; i[5]++ ) {
	for ( i[6]=-1; i[6]<=+1; i[6]++ ) {
	for ( i[7]=-1; i[7]<=+1; i[7]++ ) {
	for ( i[8]=-1; i[8]<=+1; i[8]++ ) {

		UnitCellTransformation *tfn;
		UnitCell *nc;
		int j, k;
		int l = 0;

		for ( j=0; j<3; j++ )
			for ( k=0; k<3; k++ )
				intmat_set(m, j, k, i[l++]);

		if ( intmat_det(m) != +1 ) continue;

		tfn = tfn_from_intmat(m);
		nc = cell_transform(b, tfn);

		if ( cells_are_similar(a, nc) ) {
			*pmb = m;
			tfn_free(tfn);
			cell_free(nc);
			return 1;
		}

		tfn_free(tfn);
		cell_free(nc);

	}
	}
	}
	}
	}
	}
	}
	}
	}

	intmat_free(m);
	return 0;
}


/* Try all combinations of crystals from the two patterns, in the hope of
 * finding the start of a rotation series */
static int try_all(struct image *a, struct image *b, int *c1, int *c2,
                   IntegerMatrix **m2)
{
	int i, j;

	for ( i=0; i<a->n_crystals; i++ ) {
		for ( j=0; j<b->n_crystals; j++ ) {
			if ( gatinator(crystal_get_cell(a->crystals[i]),
			               crystal_get_cell(b->crystals[j]), m2) )
			{
				*c1 = i;
				*c2 = j;
				return 1;
			}
		}
	}

	return 0;
}


/* Try to continue the rotation series from crystal c1 in image a, using any
 * crystal from image b */
static int try_cont(struct image *a, struct image *b, int c1, int *c2,
                   IntegerMatrix **m2)
{
	int j;

	for ( j=0; j<b->n_crystals; j++ ) {
		if ( gatinator(crystal_get_cell(a->crystals[c1]),
		               crystal_get_cell(b->crystals[j]), m2) )
		{
			*c2 = j;
			return 1;
		}
	}

	return 0;
}

static void dump(struct image *win, signed int *ser, IntegerMatrix **mat,
                 int window_len, int pos)
{
	int i;

	for ( i=0; i<pos; i++ ) {
		free_all_crystals(&win[i]);
		intmat_free(mat[i]);
		free(win[i].filename);
	}

	memmove(win, &win[pos], (window_len-pos)*sizeof(struct image *));
	memmove(ser, &ser[pos], (window_len-pos)*sizeof(signed int));
	memmove(mat, &mat[pos], (window_len-pos)*sizeof(IntegerMatrix *));
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
	signed int *ser;
	IntegerMatrix **mat;
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
	ser = calloc(window_len, sizeof(int));
	mat = calloc(window_len, sizeof(IntegerMatrix *));
	if ( (win == NULL) || (ser == NULL) || (mat == NULL) ) {
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
			ser[0] = -1;
			pos++;
			continue;
		}

		if ( ser[pos-1] == -1 ) {
			if ( try_all(&win[pos-1], cur, &c1, &c2, &mat[pos]) ) {
				ser[pos-1] = c1;
				ser[pos] = c2;
			} else {
				ser[pos] = -1;
				mat[pos] = NULL;
			}
		} else {
			if ( try_cont(&win[pos-1], cur, ser[pos-1], &c2,
				     &mat[pos]) )
			{
				ser[pos] = c2;
			} else {
				ser[pos] = -1;
				mat[pos] = NULL;
			}
		}

		if ( ser[pos-1] != -1 ) {
			printf("*");
		} else {
			printf("-");
		}
		fflush(stdout);

		if ( ser[0] == -1 ) {
			dump(win, ser, mat, window_len, 1);
			pos--;
		}

		if ( (pos > 0) && (ser[pos] == -1) && (ser[pos-1] == -1) ) {
			/* Series ready to process */
			process_series(win, ser, mat, pos-1);
			dump(win, ser, mat, window_len, pos);
			pos = 0;
		}

		pos++;
		if ( pos == window_len ) {
			window_len *= 2;
			win = realloc(win, window_len*sizeof(struct image));
			ser = realloc(ser, window_len*sizeof(signed int));
			mat = realloc(mat, window_len*sizeof(IntegerMatrix *));
			if ( (win == NULL) || (ser == NULL) || (mat == NULL) ) {
				ERROR("Failed to expand series buffers\n");
				return 1;
			}
		}

	} while ( 1 );
	printf("\n");

	close_stream(st);

	/* Final series to process? */
	if ( pos > 0 ) {
		process_series(win, ser, mat, pos);
		dump(win, ser, mat, window_len, 1);
	}
	free(win);
	free(ser);
	free(mat);

	return 0;
}
