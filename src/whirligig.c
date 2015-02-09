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

/* Maximum number of series which can overlap at once */
#define MAX_SER 8

struct window
{
	struct image   *img;
	int             ws;
	int             add_ptr;   /* First empty slot (for adding frames) */
	int             join_ptr;  /* First unjoined slot */

	int            *ser[MAX_SER];
	IntegerMatrix **mat[MAX_SER];
};


struct series_stats
{
	int n_series;            /* Number of series */
	int in_series;           /* Number of frames with at least one series */
	int max_series_length;   /* Length of longest series */
	int total_series_steps;  /* For calculating mean series length */
	int late_frames;         /* Number of frames which arrived too late */
	int missed_frames;       /* Number of frames which scrolled out of the
	                          * window before they could be analysed */
};


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
                           IntegerMatrix **mat, int len, const char *outdir,
                           struct series_stats *ss)
{
	int i;
	RefList **p;
	char filename[256];
	FILE *fh;
	int snum = ss->n_series;

	printf("\n");
	STATUS("Found a rotation series of %i views\n", len);

	ss->n_series++;
	if ( len > ss->max_series_length ) ss->max_series_length = len;
	ss->total_series_steps += len;

	snprintf(filename, 256, "%s/series-%i.log", outdir, snum);
	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Failed to open log file '%s'\n", filename);
		goto out;
	}

	p = calloc(len, sizeof(RefList *));
	if ( p == NULL ) return;

	fprintf(fh, "%i frames in series\n\n", len);
	fprintf(fh, "   # Serial Filename   EventID   Crystal\n");
	for ( i=0; i<len; i++ ) {
		Crystal *cr = images[i].crystals[ser[i]];
		fprintf(fh, "%4i %5i %s %s %i\n", i, images[i].serial,
		                             images[i].filename,
		                             get_event_string(images[i].event),
		                             ser[i]);
		p[i] = transform_reflections(crystal_get_reflections(cr),
		                             mat[i]);
	}

	for ( i=1; i<len; i++ ) {
		STATUS("%i -> %i: %i common reflections\n",
		       i-1, i, find_common_reflections(p[i-1], p[i]));
	}

	for ( i=0; i<len; i++ ) {
		reflist_free(p[i]);
	}
	free(p);
	fclose(fh);

out:
	for ( i=0; i<len; i++ ) {
		ser[i] = -1;
		intmat_free(mat[i]);
	}
}


static void count_series_frames(int **ser, int ser_start, int ser_len,
                                struct series_stats *ss)
{
	int i;

	for ( i=0; i<ser_len; i++ ) {

		int j;
		int clean = 1;

		for ( j=0; j<MAX_SER; j++ ) {
			if ( ser[j][ser_start+i] != -1 ) clean = 0;
		}

		if ( clean ) {
			ss->in_series++;
		}

	}
}


static void find_ser(struct window *win, int sn, int is_last_frame,
                     struct series_stats *ss, const char *outdir)
{
	int i;
	int ser_len = 0;
	int ser_start = 0;
	int in_series = 0;

	assert(win->join_ptr <= win->ws);

	for ( i=0; i<win->join_ptr; i++ ) {

		if ( in_series && win->ser[sn][i] == -1 ) {

			process_series(win->img+ser_start,
			               win->ser[sn]+ser_start,
			               win->mat[sn]+ser_start,
			               ser_len, outdir, ss);

			count_series_frames(win->ser, ser_start, ser_len, ss);

			in_series = 0;

		}

		if ( win->ser[sn][i] != -1 ) {
			if ( in_series ) {
				ser_len++;
			} else {
				ser_start = i;
				ser_len = 1;
				in_series = 1;
			}
		}

	}

	if ( is_last_frame && (ser_len > 1) ) {
		process_series(win->img+ser_start,
		               win->ser[sn]+ser_start,
		               win->mat[sn]+ser_start,
		               ser_len, outdir, ss);
		count_series_frames(win->ser, ser_start, ser_len, ss);
	}
}


static void find_and_process_series(struct window *win, int is_last_frame,
                                    struct series_stats *ss, const char *outdir)
{
	int i;

	for ( i=0; i<MAX_SER; i++ ) {
		find_ser(win, i, is_last_frame, ss, outdir);
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


static int crystal_used(struct window *win, int pos, int cn)
{
	int i;

	for ( i=0; i<MAX_SER; i++ ) {
		if ( win->ser[i][pos] == cn ) return 1;
	}

	return 0;
}


static IntegerMatrix *try_all(struct window *win, int n1, int n2,
                              int *c1, int *c2)
{
	int i, j;
	IntegerMatrix *m;
	struct image *i1;
	struct image *i2;

	assert(n1 >= 0);
	assert(n2 >= 0);
	assert(n1 < win->ws);
	assert(n2 < win->ws);

	i1 = &win->img[n1];
	i2 = &win->img[n2];

	for ( i=0; i<i1->n_crystals; i++ ) {
	for ( j=0; j<i2->n_crystals; j++ ) {

		if ( gatinator(crystal_get_cell(i1->crystals[i]),
		               crystal_get_cell(i2->crystals[j]), &m) )
		{
			if ( !crystal_used(win, n1, i)
			  && !crystal_used(win, n2, j) )
			{
				*c1 = i;
				*c2 = j;
				return m;
			}
		}

	}
	}

	return NULL;
}


/* Return a series number which can be used at the current join_ptr */
static int find_available_series(struct window *win)
{
	int i;

	for ( i=0; i<MAX_SER; i++ ) {

		/* Series must not be in use at the moment */
		if ( win->ser[i][win->join_ptr] != -1 ) continue;

		/* Series must not have been in use recently */
		if ( win->join_ptr > 0 ) {
			if ( win->ser[i][win->join_ptr-1] != -1 ) continue;
		}

		if ( win->join_ptr > 1 ) {
			if ( win->ser[i][win->join_ptr-2] != -1 ) continue;
		}


		return i;

	}

	ERROR("Too many overlapping series!\n");
	abort();
}


/* Try to fit p1 in with p2 */
static int try_join(struct window *win, int sn)
{
	int j;
	Crystal *cr;
	UnitCell *ref;
	UnitCellTransformation *tfn;
	const int sp = win->join_ptr - 1;

	/* Get the appropriately transformed cell from the last crystal in this
	 * series */
	tfn = tfn_from_intmat(win->mat[sn][sp]);
	cr = win->img[sp].crystals[win->ser[sn][sp]];
	ref = cell_transform(crystal_get_cell(cr), tfn);
	tfn_free(tfn);

	for ( j=0; j<win->img[win->join_ptr].n_crystals; j++ ) {
		Crystal *cr2;
		cr2 = win->img[win->join_ptr].crystals[j];
		if ( gatinator(ref, crystal_get_cell(cr2),
		               &win->mat[sn][win->join_ptr]) ) {
			win->ser[sn][win->join_ptr] = j;
			cell_free(ref);
			return 1;
		}
	}

	cell_free(ref);

	return 0;
}


static void connect_series(struct window *win)
{
	while ( win->join_ptr < win->ws ) {

		int i;

		if ( win->join_ptr == 0 ) {
			win->join_ptr++;
			continue;
		}

		/* Stop if we found a missing frame */
		if ( win->img[win->join_ptr].serial == 0 ) break;

		/* Try to join this frame to each of the active series */
		if ( win->join_ptr > 1 ) {
			for ( i=0; i<MAX_SER; i++ ) {
				if ( win->ser[i][win->join_ptr-1] != -1 ) {
					try_join(win, i);
				}
			}
		}

		/* Try to nucleate a new series here */
		if ( (win->join_ptr > 0)
		  && (win->img[win->join_ptr-1].serial != 0) )
		{
			IntegerMatrix *m;
			int c1, c2;
			m = try_all(win, win->join_ptr-1, win->join_ptr,
			            &c1, &c2);
			if ( m != NULL ) {
				int sn = find_available_series(win);
				win->ser[sn][win->join_ptr-1] = c1;
				win->mat[sn][win->join_ptr-1] = intmat_identity(3);
				win->ser[sn][win->join_ptr] = c2;
				win->mat[sn][win->join_ptr] = m;
			}
		}

		win->join_ptr++;

	};
}


static int series_fills_window(struct window *win)
{
	int i;
	int cont[MAX_SER];

	for ( i=0; i<MAX_SER; i++ ) cont[i] = 1;

	for ( i=0; i<win->ws; i++ ) {
		int j;
		for ( j=0; j<MAX_SER; j++ ) {
			if ( (win->img[i].serial != 0)
			  && (win->ser[j][i] == -1) ) {
				cont[j] = 0;
			}
		}
	}

	for ( i=0; i<MAX_SER; i++ ) {
		if ( cont[i] ) return 1;
	}

	return 0;
}


static void add_to_window(struct image *cur, struct window *win,
                          struct series_stats *ss)
{
	int pos;

	pos = cur->serial - win->img[win->add_ptr-1].serial;
	pos += win->add_ptr - 1;

	if ( pos < 0 ) {
		/* Frame arrived too late */
		ss->late_frames++;
		return;
	}

	if ( pos >= win->ws ) {

		int sf, i;

		sf = (pos - win->ws) + 1;

		if ( series_fills_window(win) ) {

			win->ws += sf;
			win->img = realloc(win->img,
			                   win->ws*sizeof(struct image));
			if ( win->img == NULL ) {
				ERROR("Failed to expand series buffers\n");
				exit(1);
			}
			for ( i=0; i<MAX_SER; i++ ) {
				win->ser[i] = realloc(win->ser[i],
				               win->ws*sizeof(signed int));
				win->mat[i] = realloc(win->mat[i],
				               win->ws*sizeof(IntegerMatrix *));
				if ( (win->ser[i] == NULL)
				  || (win->mat[i] == NULL) )
				{
					ERROR("Failed to expand buffers\n");
					exit(1);
				}
			}

		} else {

			pos -= sf;
			if ( sf > win->join_ptr ) {

				int i;

				for ( i=0; i<sf-win->join_ptr; i++ ) {
					if ( win->img[i].serial != 0 ) {
						ss->missed_frames++;
					}
				}

				win->join_ptr = 0;

			} else {
				win->join_ptr -= sf;
			}

			if ( sf > win->ws ) {
				sf = win->ws;
			}

			for ( i=0; i<sf; i++ ) {
				if ( win->img[i].serial != 0 ) {
					free_all_crystals(&win->img[i]);
				}
			}

			memmove(win->img, win->img+sf,
			        (win->ws-sf)*sizeof(struct image));

			for ( i=0; i<MAX_SER; i++ ) {
				memmove(win->ser[i], win->ser[i]+sf,
					(win->ws-sf)*sizeof(signed int));
				memmove(win->mat[i], win->mat[i]+sf,
					(win->ws-sf)*sizeof(IntegerMatrix *));
			}
		}

		for ( i=0; i<sf; i++ ) {
			int j;
			win->img[win->ws-sf+i].serial = 0;
			for ( j=0; j<MAX_SER; j++ ) {
				win->ser[j][win->ws-sf+i] = -1;
				win->mat[j][win->ws-sf+i] = NULL;
			}
		}

	}

	win->img[pos] = *cur;
	if ( pos >= win->add_ptr ) win->add_ptr = pos+1;
}


static void show_help(const char *s)
{
	printf("Syntax: %s <input.stream> [options]\n\n", s);
	printf(
"Find and combine rotation series.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"      --version              Print CrystFEL version number and exit.\n"
"\n"
"      --window-size=n        History size for finding connected crystals.\n"
"      --output-dir=folder    Put output files in <folder>.\n");
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
	struct window win;
	int i;
	char *rval;
	struct series_stats ss;
	int n_images = 0;

	/* Defaults */
	int default_window_size = 16;
	char *outdir = ".";
	int verbose = 0;

	/* Long options */
	const struct option longopts[] = {

		{"help",               0, NULL,               'h'},
		{"verbose",            0, NULL,               'v'},

		{"version",            0, NULL,                3 },
		{"window-size",        1, NULL,                4 },
		{"output-dir",         1, NULL,                5 },

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

			case 'v' :
			verbose = 1;
			break;

			case 3 :
			printf("CrystFEL: " CRYSTFEL_VERSIONSTRING "\n");
			printf(CRYSTFEL_BOILERPLATE"\n");
			return 0;

			case 4 :
			errno = 0;
			default_window_size = strtol(optarg, &rval, 10);
			if ( (*rval != '\0') || (default_window_size < 2) ) {
				ERROR("Invalid value for --window-size.\n");
				return 1;
			}
			break;

			case 5 :
			outdir = strdup(optarg);
			break;

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
	win.ws = default_window_size;
	win.img = calloc(win.ws, sizeof(struct image));
	if ( (win.img == NULL) ) {
		ERROR("Failed to allocate series buffers\n");
		return 1;
	}
	for ( i=0; i<win.ws; i++ ) {
		win.img[i].serial = 0;
	}
	for ( i=0; i<MAX_SER; i++ ) {

		int j;

		win.ser[i] = calloc(win.ws, sizeof(signed int));
		win.mat[i] = calloc(win.ws, sizeof(IntegerMatrix *));
		if ( (win.ser[i] == NULL) || (win.mat[i] == NULL) ) {
			ERROR("Failed to allocate series buffers\n");
			return 1;
		}
		for ( j=0; j<win.ws; j++ ) {
			win.ser[i][j] = -1;
			win.mat[i][j] = NULL;
		}

	}

	ss.n_series = 0;
	ss.in_series = 0;
	ss.max_series_length = 0;
	ss.total_series_steps = 0;
	ss.late_frames = 0;
	ss.missed_frames = 0;

	win.add_ptr = 1;  /* Horrendous bodge, but it works */
	win.join_ptr = 0;
	do {

		struct image cur;

		cur.div = NAN;
		cur.bw = NAN;
		cur.det = NULL;
		if ( read_chunk_2(st, &cur, STREAM_READ_REFLECTIONS
		                            | STREAM_READ_UNITCELL) != 0 ) {
			break;
		}

		if ( verbose ) printf("\n\nIncoming serial %i\n", cur.serial);

		if ( isnan(cur.div) || isnan(cur.bw) ) {
			ERROR("Chunk doesn't contain beam parameters.\n");
			return 1;
		}

		if ( cur.serial < 1 ) {
			ERROR("Serial numbers must be greater than zero.\n");
			return 1;
		}

		add_to_window(&cur, &win, &ss);
		connect_series(&win);

		if ( verbose ) {
			for ( i=0; i<win.ws; i++ ) {
				int j;
				printf("%3i %4i %c %c", i, win.img[i].serial,
				       (i==win.add_ptr)?'*':' ',
				       (i==win.join_ptr)?'<':' ');
				for ( j=0; j<MAX_SER; j++ ) {
					printf(" %3i", win.ser[j][i]);
				}
				printf(" %s\n", win.img[i].filename);
			}
			printf("(%3i) %c %c\n", i,
				       (win.ws==win.add_ptr)?'*':' ',
				       (win.ws==win.join_ptr)?'<':' ');

		}

		find_and_process_series(&win, 0, &ss, outdir);

		display_progress(n_images++);

	} while ( 1 );
	display_progress(n_images);
	printf("\n");

	close_stream(st);

	find_and_process_series(&win, 1, &ss, outdir);

	STATUS("-----------------------------------------------------\n");
	STATUS("            Number of frames processed: %i\n", n_images);
	STATUS("              Frames arriving too late: %i", ss.late_frames);
	if ( ss.late_frames > 0 ) {
		STATUS(" (consider increasing the window size)");
	}
	STATUS("\n");
	STATUS(" Frames leaving window before analysis: %i", ss.missed_frames);
	if ( ss.missed_frames > 0 ) {
		STATUS(" (consider increasing the window size)");
	}
	STATUS("\n");
	STATUS("             Number of rotation series: %i\n", ss.n_series);
	STATUS("                 Average series length: %-6.2f frames\n",
	       (double)ss.total_series_steps/ss.n_series);
	STATUS("              Length of longest series: %-6i frames\n",
	       ss.max_series_length);
	STATUS("            Number of frames in series: %i\n", ss.in_series);
	STATUS("          Fraction of frames in series: %-6.2f %%\n",
	       (double)ss.in_series*100.0 / n_images);

	return 0;
}
