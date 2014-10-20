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
#include "reflist.h"
#include "reflist-utils.h"


static void do_op(const IntegerMatrix *op,
                  signed int h, signed int k, signed int l,
                  signed int *he, signed int *ke, signed int *le)
{
	signed int v[3];
	signed int *ans;

	v[0] = h;  v[1] = k;  v[2] = l;

	ans = intmat_intvec_mult(op, v);
	assert(ans != NULL);

	*he = ans[0];  *ke = ans[1];  *le = ans[2];
	free(ans);
}


static RefList *transform_reflections(RefList *in, IntegerMatrix *m)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *out;

	if ( m == NULL ) return copy_reflist(in);

	out = reflist_new();

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l, he, ke, le;
		Reflection *n;
		get_indices(refl, &h, &k, &l);
		do_op(m, h, k, l, &he, &ke, &le);
		n = add_refl(out, he, ke, le);
		copy_data(n, refl);
	}

	return out;
}


static int find_common_reflections(RefList *list1, RefList *list2)
{
	Reflection *refl1;
	RefListIterator *iter;
	int ncom = 0;

	for ( refl1 = first_refl(list1, &iter);
	      refl1 != NULL;
	      refl1 = next_refl(refl1, iter) )
	{
		signed int h, k, l;
		Reflection *refl2;
		get_indices(refl1, &h, &k, &l);
		refl2 = find_refl(list2, h, k, l);
		if ( refl2 == NULL ) continue;
		ncom++;
	}

	return ncom;
}


static void process_series(struct image *images, signed int *ser,
                           IntegerMatrix **mat, int len)
{
	int i;
	RefList **p;

	printf("\n");
	STATUS("Found a rotation series of %i views\n", len);

	p = calloc(len, sizeof(RefList *));
	if ( p == NULL ) return;

	for ( i=0; i<len; i++ ) {
		Crystal *cr = images[i].crystals[ser[i]];
		p[i] = transform_reflections(crystal_get_reflections(cr),
		                             mat[i]);
	}

	for ( i=1; i<len; i++ ) {
		STATUS("%i -> %i: %i common reflections\n",
		       i-1, i, find_common_reflections(p[i-1], p[i]));
	}

	for ( i=0; i<len; i++ ) {
		reflist_free(p[i]);
		ser[i] = -1;
	}
	free(p);
}


static void check_for_series(struct image *win, signed int *ser,
                             IntegerMatrix **m, int ws, int is_last_frame)
{
	int i;
	int ser_len = 0;
	int ser_start = 0;

	for ( i=0; i<ws; i++ ) {
		if ( (win[i].serial != 0) && (ser[i] == -1) ) {
			if ( ser_len > 2 ) {
				process_series(win+ser_start, ser+ser_start,
				               m+ser_start, ser_len);
			}
		}
		if ( (win[i].serial != 0) && (ser[i] != -1) ) {
			ser_len++;
		} else {
			ser_start = i;
			ser_len = 0;
		}
		//STATUS("%3i: serial %i, series %i, matrix %p, len %i\n",
		//       i, win[i].serial, ser[i], m[i], ser_len);
	}

	if ( is_last_frame && (ser_len > 2) ) {
		process_series(win+ser_start, ser+ser_start,
		               m+ser_start, ser_len);
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


static int try_all(struct image *win, signed int *ser, IntegerMatrix **m,
                   int ws, signed int p1, signed int p2)
{
	int i, j;

	for ( i=0; i<win[p1].n_crystals; i++ ) {
	for ( j=0; j<win[p2].n_crystals; j++ ) {

		if ( gatinator(crystal_get_cell(win[p1].crystals[i]),
		               crystal_get_cell(win[p2].crystals[j]), &m[p2]) )
		{
			ser[p1] = i;
			ser[p2] = j;
			m[p1] = intmat_identity(3);
			return 1;
		}

	}
	}

	return 0;
}


/* Try to fit p1 in with p2 */
static int try_join(struct image *win, signed int *ser, IntegerMatrix **m,
                    int ws, signed int p1, signed int p2)
{
	int j;
	UnitCell *ref;
	UnitCellTransformation *tfn;

	if ( ser[p2] == -1 ) {
		return try_all(win, ser, m, ws, p1, p2);
	}

	tfn = tfn_from_intmat(m[p2]);
	ref = cell_transform(crystal_get_cell(win[p2].crystals[ser[p2]]), tfn);
	tfn_free(tfn);

	for ( j=0; j<win[p1].n_crystals; j++ ) {
		if ( gatinator(ref, crystal_get_cell(win[p1].crystals[j]),
		               &m[p1]) ) {
			ser[p1] = j;
			cell_free(ref);
			return 1;
		}
	}

	cell_free(ref);

	return 0;
}


static void try_connect(struct image *win, signed int *ser, IntegerMatrix **m,
                        int ws, signed int pos)
{
	/* Try to connect to the left */
	if ( (pos > 0) && (win[pos-1].serial != 0) ) {
		try_join(win, ser, m, ws, pos, pos-1);
	}

	/* Try to connect to the right */
	if ( (pos < ws-1) && (win[pos+1].serial != 0) ) {
		try_join(win, ser, m, ws, pos, pos+1);
	}
}


static int series_fills_window(signed int *ser, int ws)
{
	int i;

	for ( i=0; i<ws; i++ ) {
		if ( ser[i] == -1 ) return 0;
	}

	return 1;
}


static int add_to_window(struct image *cur, struct image **pwin,
                         signed int **pser, IntegerMatrix ***pmat,
                         int *pws)
{
	int sf, pos;
	struct image *win = *pwin;
	signed int *ser = *pser;
	IntegerMatrix **mat = *pmat;
	int ws = *pws;

	if ( cur->serial > win[ws-1].serial ) {

		int i;

		if ( series_fills_window(ser, ws) ) {

			ws++;
			win = realloc(win, ws*sizeof(struct image));
			ser = realloc(ser, ws*sizeof(signed int));
			mat = realloc(mat, ws*sizeof(IntegerMatrix *));
			if ( (win == NULL) || (ser == NULL) || (mat == NULL) ) {
				ERROR("Failed to expand series buffers\n");
				exit(1);
			}
			*pws = ws;
			*pwin = win;
			*pser = ser;
			*pmat = mat;

			pos = ws-1;

		} else {

			sf = cur->serial - win[ws-1].serial;

			memmove(win, win+sf, (ws-sf)*sizeof(struct image));
			memmove(ser, ser+sf, (ws-sf)*sizeof(signed int));
			memmove(mat, mat+sf, (ws-sf)*sizeof(IntegerMatrix *));

			for ( i=0; i<sf; i++ ) {
				win[ws-sf+i].serial = 0;
				ser[ws-sf+i] = -1;
			}

			pos = ws - 1;

		}

	} else {

		pos = ws-(win[ws-1].serial - cur->serial)-1;

	}

	win[pos] = *cur;

	return pos;
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


static void display_progress(int n_images)
{
	if ( !isatty(STDERR_FILENO) ) return;
	if ( tcgetpgrp(STDERR_FILENO) != getpgrp() ) return;

	pthread_mutex_lock(&stderr_lock);
	fprintf(stderr, "\r%i images processed.", n_images);
	pthread_mutex_unlock(&stderr_lock);

	fflush(stdout);
}


int main(int argc, char *argv[])
{
	int c;
	Stream *st;
	int polarisation = 1;
	struct image *win;
	signed int *ser;
	IntegerMatrix **mat;
	int default_window_size = 16;
	int ws, i;
	int n_images = 0;

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

	/* Allocate initial window */
	ws = default_window_size;
	win = calloc(ws, sizeof(struct image));
	ser = calloc(ws, sizeof(signed int));
	mat = calloc(ws, sizeof(IntegerMatrix *));
	if ( (win == NULL) || (ser == NULL) || (mat == NULL) ) {
		ERROR("Failed to allocate series buffers\n");
		return 1;
	}

	for ( i=0; i<ws; i++ ) {
		win[i].serial = 0;
		ser[i] = -1;
	}

	do {

		struct image cur;
		int pos;

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

		if ( cur.serial < 1 ) {
			ERROR("Serial numbers must be greater than zero.\n");
			return 1;
		}

		pos = add_to_window(&cur, &win, &ser, &mat, &ws);
		try_connect(win, ser, mat, ws, pos);
		check_for_series(win, ser, mat, ws, 0);

		display_progress(n_images++);

	} while ( 1 );
	display_progress(n_images);
	printf("\n");

	close_stream(st);

	check_for_series(win, ser, mat, ws, 1);

	free(win);
	free(ser);
	free(mat);

	return 0;
}
