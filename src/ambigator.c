/*
 * ambigator.c
 *
 * Resolve indexing ambiguities
 *
 * Copyright © 2014 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2014 Wolfgang Brehm
 *
 * Authors:
 *   2014 Thomas White <taw@physics.org>
 *   2014 Wolfgang Brehm <wolfgang.brehm@gmail.com>
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
#include <symmetry.h>
#include <stream.h>
#include <reflist.h>
#include <reflist-utils.h>
#include <cell.h>
#include <cell-utils.h>


static void show_help(const char *s)
{
	printf("Syntax: %s [options] input.stream\n\n", s);
	printf(
"Resolve indexing ambiguities.\n"
"\n"
"  -h, --help                 Display this help message.\n"
"\n"
"  -o, --output=<filename>    Output stream.\n"
"  -y, --symmetry=<sym>       Actual (\"target\") symmetry.\n"
"  -w <sym>                   Apparent (\"source\" or \"twinned\") symmetry.\n"
"  -n, --iterations=<n>       Iterate <n> times.\n"
"      --highres=<n>          High resolution cutoff in A.\n"
"      --lowres=<n>           Low resolution cutoff in A.\n"
);
}

#define SERIAL(h, k, l) ((((h)+512)<<20) + (((k)+512)<<10) + ((l)+512))
#define GET_H(serial) ((((serial) & 0x3ff00000)>>20)-512)
#define GET_K(serial) ((((serial) & 0x000ffc00)>>10)-512)
#define GET_L(serial) (((serial) & 0x000003ff)-512)

struct flist
{
	int n;
	unsigned int *s;
	float *i;
};


static struct flist *asymm_and_merge(RefList *in, const SymOpList *sym,
                                     UnitCell *cell, double rmin, double rmax)
{
	Reflection *refl;
	RefListIterator *iter;
	RefList *asym;
	struct flist *f;
	int n;

	asym = reflist_new();
	if ( asym == NULL ) return NULL;

	for ( refl = first_refl(in, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		signed int ha, ka, la;
		Reflection *cr;
		double res;

		get_indices(refl, &h, &k, &l);

		res = 2.0*resolution(cell, h, k, l);
		if ( res < rmin ) continue;
		if ( res > rmax ) continue;

		get_asymm(sym, h, k, l, &ha, &ka, &la);

		cr = find_refl(asym, ha, ka, la);
		if ( cr == NULL ) {
			cr = add_refl(asym, ha, ka, la);
			assert(cr != NULL);
			copy_data(cr, refl);
		} else {
			const double i = get_intensity(cr);
			const int r = get_redundancy(cr);
			set_intensity(cr, (r*i + get_intensity(refl))/(r+1));
			set_redundancy(cr, r+1);
		}
	}

	f = malloc(sizeof(struct flist));
	if ( f == NULL ) {
		ERROR("Failed to allocate flist\n");
		return NULL;
	}

	n = num_reflections(asym);
	f->s = malloc(n*sizeof(unsigned int));
	f->i = malloc(n*sizeof(float));
	if ( (f->s == NULL) || (f->i == NULL) ) {
		ERROR("Failed to allocate flist\n");
		return NULL;
	}

	f->n = 0;
	for ( refl = first_refl(asym, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);
		f->s[f->n] = SERIAL(h, k, l);
		f->i[f->n] = get_intensity(refl);
		f->n++;
	}
	reflist_free(asym);
	assert(f->n == n);

	return f;
}


static float corr(struct flist *a, struct flist *b, int *pn)
{
	float s_xy = 0.0;
	float s_x = 0.0;
	float s_y = 0.0;
	float s_x2 = 0.0;
	float s_y2 = 0.0;
	int n = 0;
	float t1, t2;
	int ap = 0;
	int bp = 0;
	int done = 0;

	while ( 1 ) {

		while ( a->s[ap] > b->s[bp] ) {
			if ( ++bp == b->n ) {
				done = 1;
				break;
			}
		}
		if ( done ) break;

		while ( a->s[ap] < b->s[bp] ) {
			if ( ++ap == a->n ) {
				done = 1;
				break;
			}
		}
		if ( done ) break;

		if ( a->s[ap] == b->s[bp] ) {

			float aint, bint;

			aint = a->i[ap];
			bint = b->i[bp];

			s_xy += aint*bint;
			s_x += aint;
			s_y += bint;
			s_x2 += aint*aint;
			s_y2 += bint*bint;
			n++;

			if ( ++ap == a->n ) break;
			if ( ++bp == b->n ) break;

		}

	}
	*pn = n;

	t1 = s_x2 - s_x*s_x / n;
	t2 = s_y2 - s_y*s_y / n;

	if ( (t1 <= 0.0) || (t2 <= 0.0) ) return 0.0;

	return (s_xy - s_x*s_y/n) / sqrt(t1*t2);
}


static void detwin(struct flist **crystals, int n_crystals, SymOpList *amb,
                   int *assignments)
{
	int i;
	int nch = 0;
	float mf = 0.0;
	int nmf = 0;
	int ndud = 0;

	for ( i=0; i<n_crystals; i++ ) {

		int j;
		float f = 0.0;
		float g = 0.0;;
		int p = 0;
		int q = 0;

		progress_bar(i, n_crystals-1, "Calculating");

		for ( j=0; j<n_crystals; j++ ) {

			float cc;
			int n;

			if ( i == j ) continue;

			cc = corr(crystals[i], crystals[j], &n);

			if ( n < 3 ) continue;

			if ( assignments[i] == assignments[j] ) {
				f += cc;
				p++;
			} else {
				g += cc;
				q++;
			}

		}

		if ( (p==0) || (q==0) ) {
			ndud++;
			continue;
		}

		f /= p;
		g /= q;

		mf += f;
		nmf++;

		if ( f < g ) {
			assignments[i] = 1 - assignments[i];
			nch++;
		}

	}

	if ( ndud > 0 ) {
		STATUS("Warning: %i crystals had no correlation\n", ndud);
	}

	STATUS("Mean f = %f, changed %i assignments this time.\n", mf/nmf, nch);
}


int main(int argc, char *argv[])
{
	int c;
	const char *infile;
	char *outfile = NULL;
	char *s_sym_str = NULL;
	SymOpList *s_sym;
	char *e_sym_str = NULL;
	SymOpList *e_sym;
	SymOpList *amb;
	int n_iter = 1;
	int n_crystals, n_chunks, max_crystals;
	int n_dif;
	struct flist **crystals;
	Stream *st;
	int i;
	int *assignments;
	int *orig_assignments;
	gsl_rng *rng;
	float highres, lowres;
	double rmin = 0.0;  /* m^-1 */
	double rmax = INFINITY;  /* m^-1 */

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},

		{"highres",            1, NULL,                2},
		{"lowres",             1, NULL,                3},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:y:n:w:",
	                        longopts, NULL)) != -1)
	{

		switch (c) {

			case 'h' :
			show_help(argv[0]);
			return 0;

			case 'o' :
			outfile = strdup(optarg);
			break;

			case 'y' :
			s_sym_str = strdup(optarg);
			break;

			case 'w' :
			e_sym_str = strdup(optarg);
			break;

			case 'n' :
			n_iter = atoi(optarg);
			break;

			case 2 :
			if ( sscanf(optarg, "%e", &highres) != 1 ) {
				ERROR("Invalid value for --highres\n");
				return 1;
			}
			rmax = 1.0 / (highres/1e10);
			break;

			case 3 :
			if ( sscanf(optarg, "%e", &lowres) != 1 ) {
				ERROR("Invalid value for --lowres\n");
				return 1;
			}
			rmin = 1.0 / (lowres/1e10);
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
		ERROR("Please provide exactly one stream filename.\n");
		return 1;
	}

	infile = argv[optind++];
	st = open_stream_for_read(infile);
	if ( st == NULL ) {
		ERROR("Failed to open input stream '%s'\n", infile);
		return 1;
	}

	/* Sanitise output filename */
	if ( outfile == NULL ) {
		outfile = strdup("partialator.hkl");
	}

	if ( s_sym_str == NULL ) {
		ERROR("You must specify the input symmetry (with -y)\n");
		return 1;
	}
	s_sym = get_pointgroup(s_sym_str);
	free(s_sym_str);

	if ( e_sym_str == NULL ) {
		e_sym = NULL;
		amb = NULL;
	} else {
		e_sym = get_pointgroup(e_sym_str);
		free(e_sym_str);
		if ( e_sym == NULL ) return 1;
		amb = get_ambiguities(s_sym, e_sym);
		if ( amb == NULL ) return 1;
		STATUS("Ambiguity operations:\n");
		describe_symmetry(amb);
	}

	crystals = NULL;
	n_crystals = 0;
	max_crystals = 0;
	n_chunks = 0;
	do {

		struct image cur;
		int i;

		cur.det = NULL;

		if ( read_chunk(st, &cur) != 0 ) {
			break;
		}

		image_feature_list_free(cur.features);

		for ( i=0; i<cur.n_crystals; i++ ) {

			Crystal *cr;
			RefList *list;
			UnitCell *cell;

			cr = cur.crystals[i];
			cell = crystal_get_cell(cr);

			if ( n_crystals == max_crystals ) {

				struct flist **crystals_new;
				size_t nsz;

				nsz = (max_crystals+1024)*sizeof(struct flist *);
				crystals_new = realloc(crystals, nsz);
				if ( crystals_new == NULL ) {
					fprintf(stderr, "Failed to allocate "
					        "memory for crystals.\n");
					break;
				}

				max_crystals += 1024;
				crystals = crystals_new;

			}

			list = crystal_get_reflections(cr);
			crystals[n_crystals] = asymm_and_merge(list, s_sym,
			                                       cell,
			                                       rmin, rmax);
			cell_free(cell);
			n_crystals++;

			reflist_free(list);

		}

		fprintf(stderr, "Loaded %i crystals from %i chunks\r",
		        n_crystals, ++n_chunks);

	} while ( 1 );
	fprintf(stderr, "\n");

	close_stream(st);

	assignments = malloc(n_crystals*sizeof(int));
	if ( assignments == NULL ) {
		ERROR("Couldn't allocate memory for assignments.\n");
		return 1;
	}

	orig_assignments = malloc(n_crystals*sizeof(int));
	if ( orig_assignments == NULL ) {
		ERROR("Couldn't allocate memory for original assignments.\n");
		return 1;
	}

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	for ( i=0; i<n_crystals; i++ ) {
		assignments[i] = (random_flat(rng, 1.0) > 0.5);
		orig_assignments[i] = assignments[i];
	}

	for ( i=0; i<n_iter; i++ ) {
		detwin(crystals, n_crystals, amb, assignments);
	}

	n_dif = 0;
	for ( i=0; i<n_crystals; i++ ) {
		if ( orig_assignments[i] != assignments[i] ) n_dif++;
	}
	STATUS("%i assignments are different from their starting values.\n",
	       n_dif);

	free(assignments);
	free(crystals);
	gsl_rng_free(rng);

	return 0;
}
