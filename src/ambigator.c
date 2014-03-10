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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>

#include <image.h>
#include <utils.h>
#include <symmetry.h>
#include <stream.h>
#include <reflist.h>
#include <reflist-utils.h>
#include <cell.h>
#include <cell-utils.h>
#include <thread-pool.h>


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
"      --end-assignments=<fn> Save end assignments to file <fn>.\n"
"      --fg-graph=<fn>        Save f and g correlation values to file <fn>.\n"
"      --ncorr=<n>            Use <n> correlations per crystal.  Default 1000\n"
"      --stop-after=<n>       Use at most the first <n> crystals.\n"
"  -j <n>                     Use <n> threads for CC calculation.\n"
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

	unsigned int *s_reidx;
	float *i_reidx;
};


static struct flist *asymm_and_merge(RefList *in, const SymOpList *sym,
                                     UnitCell *cell, double rmin, double rmax,
                                     SymOpList *amb)
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

		if ( amb != NULL ) {

			signed int hr, kr, lr;
			signed int hra, kra, lra;

			get_equiv(amb, NULL, 1, ha, ka, la, &hr, &kr, &lr);
			get_asymm(sym, hr, kr, lr, &hra, &kra, &lra);

			/* Skip twin-proof reflections */
			if ( (ha==hra) && (ka==kra) && (la==lra) ) {
				//STATUS("%i %i %i is twin proof\n", h, k, l);
				continue;
			}

		}

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
	f->s_reidx = malloc(n*sizeof(unsigned int));
	f->i = malloc(n*sizeof(float));
	f->i_reidx = malloc(n*sizeof(float));
	if ( (f->s == NULL) || (f->i == NULL)
	   || (f->s_reidx == NULL) || (f->i_reidx == NULL) ) {
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
	assert(f->n == n);

	if ( amb != NULL ) {

		RefList *reidx = reflist_new();
		if ( reidx == NULL ) return NULL;

		for ( refl = first_refl(asym, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			signed int hr, kr, lr;
			signed int hra, kra, lra;
			Reflection *cr;

			get_indices(refl, &h, &k, &l);
			get_equiv(amb, NULL, 1, h, k, l, &hr, &kr, &lr);
			get_asymm(sym, hr, kr, lr, &hra, &kra, &lra);

			cr = add_refl(reidx, hra, kra, lra);
			copy_data(cr, refl);
		}

		n = 0;
		for ( refl = first_refl(reidx, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			signed int h, k, l;
			get_indices(refl, &h, &k, &l);
			f->s_reidx[n] = SERIAL(h, k, l);
			f->i_reidx[n++] = get_intensity(refl);
		}

		reflist_free(reidx);
	}

	reflist_free(asym);

	return f;
}


static float corr(struct flist *a, struct flist *b, int *pn, int a_reidx)
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
	unsigned int *sa;
	float *ia;

	if ( a_reidx ) {
		sa = a->s_reidx;
		ia = a->i_reidx;
	} else {
		sa = a->s;
		ia = a->i;
	}

	if ( (a->n == 0) || (b->n == 0) ) {
		*pn = 0;
		return 0.0;
	}

	while ( 1 ) {

		while ( sa[ap] > b->s[bp] ) {
			if ( ++bp == b->n ) {
				done = 1;
				break;
			}
		}
		if ( done ) break;

		while ( sa[ap] < b->s[bp] ) {
			if ( ++ap == a->n ) {
				done = 1;
				break;
			}
		}
		if ( done ) break;

		if ( sa[ap] == b->s[bp] ) {

			float aint, bint;

			aint = ia[ap];
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


struct cc_list
{
	signed int *ind;
	float *cc;

	signed int *ind_reidx;
	float *cc_reidx;
};


struct queue_args
{
	int n_started;
	int n_finished;
	int n_to_do;
	int mean_nac;
	int nmean_nac;
	gsl_permutation *p;

	struct cc_list *ccs;
	struct flist **crystals;
	int n_crystals;
	int ncorr;
	SymOpList *amb;
	gsl_rng **rngs;
};


struct cc_job
{
	struct cc_list *ccs;
	int i;
	int mean_nac;
	int nmean_nac;

	struct flist **crystals;
	int n_crystals;
	int ncorr;
	SymOpList *amb;
	gsl_rng **rngs;
};


static void *get_task(void *vp)
{
	struct queue_args *qargs = vp;
	struct cc_job *job;

	if ( qargs->n_started == qargs->n_to_do ) return NULL;

	job = malloc(sizeof(struct cc_job));
	if ( job == NULL ) return NULL;

	job->ccs = qargs->ccs;
	job->i = qargs->n_started++;

	job->crystals = qargs->crystals;
	job->n_crystals = qargs->n_crystals;
	job->ncorr = qargs->ncorr;
	job->amb = qargs->amb;
	job->rngs = qargs->rngs;

	return job;
}


static void final(void *qp, void *wp)
{
	struct queue_args *qargs = qp;
	struct cc_job *job = wp;

	qargs->mean_nac += job->mean_nac;
	qargs->nmean_nac += job->nmean_nac;

	free(job);

	qargs->n_finished++;
	progress_bar(qargs->n_finished, qargs->n_to_do, "Calculating CCs");
}


static void work(void *wp, int cookie)
{
	struct cc_job *job = wp;
	int i = job->i;
	int k, l;
	struct cc_list *ccs = job->ccs;
	struct flist **crystals = job->crystals;
	int n_crystals = job->n_crystals;
	int ncorr = job->ncorr;
	SymOpList *amb = job->amb;
	int mean_nac = 0;
	int nmean_nac = 0;
	gsl_permutation *p;

	p = gsl_permutation_alloc(n_crystals);
	gsl_permutation_init(p);
	gsl_ran_shuffle(job->rngs[cookie], p->data, n_crystals, sizeof(size_t));

	ccs[i].ind = malloc(ncorr*sizeof(int));
	ccs[i].cc = malloc(ncorr*sizeof(float));
	ccs[i].ind_reidx = calloc(ncorr, sizeof(int));
	ccs[i].cc_reidx = calloc(ncorr, sizeof(float));
	if ( (ccs[i].ind==NULL) || (ccs[i].cc==NULL) ||
	     (ccs[i].ind_reidx==NULL) ||  (ccs[i].cc_reidx==NULL) ) {
		return;
	}

	k = 0;
	for ( l=0; l<n_crystals; l++ ) {

		int n;
		int j;
		float cc;

		j = gsl_permutation_get(p, l);
		if ( i == j ) continue;

		cc = corr(crystals[i], crystals[j], &n, 0);

		if ( n < 4 ) continue;

		ccs[i].ind[k] = j+1;
		ccs[i].cc[k] = cc;
		k++;

		if ( k == ncorr-1 ) break;

	}
	ccs[i].ind[k] = 0;
	mean_nac += k;
	nmean_nac++;

	if ( amb != NULL ) {

		k = 0;
		for ( l=0; l<n_crystals; l++ ) {

			int n;
			int j;
			float cc;

			j = gsl_permutation_get(p, l);
			if ( i == j ) continue;

			cc = corr(crystals[i], crystals[j], &n, 1);

			if ( n < 4 ) continue;

			ccs[i].ind_reidx[k] = j+1;
			ccs[i].cc_reidx[k] = cc;
			k++;

			if ( k == ncorr-1 ) break;

		}
		ccs[i].ind_reidx[k] = 0;
		mean_nac += k;
		nmean_nac++;

	}

	gsl_permutation_free(p);

	job->mean_nac = mean_nac;
	job->nmean_nac = nmean_nac;
}


static gsl_rng **setup_random(gsl_rng *rng, int n)
{
	gsl_rng **rngs;
	int i;

	rngs = malloc(n * sizeof(gsl_rng *));
	if ( rngs == NULL ) return NULL;

	for ( i=0; i<n; i++ ) {
		rngs[i] = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(rngs[i], gsl_rng_get(rng));
	}

	return rngs;
}


static struct cc_list *calc_ccs(struct flist **crystals, int n_crystals,
                                int ncorr, SymOpList *amb, gsl_rng *rng,
                                float *pmean_nac, int nthreads)
{
	struct cc_list *ccs;
	struct queue_args qargs;
	int i;

	assert(n_crystals >= ncorr);
	ncorr++;  /* Extra value at end for sentinel */

	qargs.rngs = setup_random(rng, nthreads);

	ccs = malloc(n_crystals*sizeof(struct cc_list));
	if ( ccs == NULL ) return NULL;

	qargs.n_started = 0;
	qargs.n_finished = 0;
	qargs.n_to_do = n_crystals;
	qargs.ccs = ccs;
	qargs.mean_nac = 0;
	qargs.nmean_nac = 0;

	qargs.crystals = crystals;
	qargs.n_crystals = n_crystals;
	qargs.ncorr = ncorr;
	qargs.amb = amb;

	run_threads(nthreads, work, get_task, final, &qargs, n_crystals,
	            0, 0, 0);

	for ( i=0; i<nthreads; i++ ) {
		gsl_rng_free(qargs.rngs[i]);
	}
	gsl_permutation_free(qargs.p);

	*pmean_nac = (float)qargs.mean_nac/qargs.nmean_nac;

	return ccs;
}


static void detwin(struct cc_list *ccs, int n_crystals, int *assignments,
                   FILE *fh, struct flist **crystals)
{
	int i;
	int nch = 0;
	float mf = 0.0;
	int nmf = 0;
	int ndud = 0;

	for ( i=0; i<n_crystals; i++ ) {

		int k;
		float f = 0.0;
		float g = 0.0;;
		int p = 0;
		int q = 0;

		for ( k=0; (ccs[i].ind[k] != 0); k++ ) {

			int j = ccs[i].ind[k]-1;
			float cc = ccs[i].cc[k];

			if ( assignments[i] == assignments[j] ) {
				f += cc;
				p++;
			} else {
				g += cc;
				q++;
			}

		}

		for ( k=0; (ccs[i].ind_reidx[k] != 0); k++ ) {

			int j = ccs[i].ind_reidx[k]-1;
			float cc = ccs[i].cc_reidx[k];

			if ( assignments[i] == assignments[j] ) {
				g += cc;
				q++;
			} else {
				f += cc;
				p++;
			}

		}

		if ( (p==0) || (q==0) ) {
			ndud++;
			continue;
		}

		f /= p;
		g /= q;

		if ( fh != NULL ) fprintf(fh, "%5.3f %5.3f\n", f, g);

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
	char *end_ass_fn = NULL;
	char *fg_graph_fn = NULL;
	char *s_sym_str = NULL;
	SymOpList *s_sym;
	char *w_sym_str = NULL;
	SymOpList *w_sym;
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
	FILE *fgfh = NULL;
	struct cc_list *ccs;
	int ncorr;
	int ncorr_set = 0;
	int stop_after = 0;
	float mean_nac;
	int n_threads = 1;

	/* Long options */
	const struct option longopts[] = {
		{"help",               0, NULL,               'h'},
		{"output",             1, NULL,               'o'},
		{"symmetry",           1, NULL,               'y'},
		{"iterations",         1, NULL,               'n'},

		{"highres",            1, NULL,                2},
		{"lowres",             1, NULL,                3},
		{"end-assignments",    1, NULL,                4},
		{"fg-graph",           1, NULL,                5},
		{"ncorr",              1, NULL,                6},
		{"stop-after",         1, NULL,                7},

		{0, 0, NULL, 0}
	};

	/* Short options */
	while ((c = getopt_long(argc, argv, "ho:y:n:w:j:",
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
			w_sym_str = strdup(optarg);
			break;

			case 'n' :
			n_iter = atoi(optarg);
			break;

			case 'j' :
			if ( sscanf(optarg, "%i", &n_threads) != 1 ) {
				ERROR("Invalid value for -j\n");
				return 1;
			}
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

			case 4 :
			end_ass_fn = strdup(optarg);
			break;

			case 5 :
			fg_graph_fn = strdup(optarg);
			break;

			case 6 :
			if ( sscanf(optarg, "%i", &ncorr) != 1 ) {
				ERROR("Invalid value for --ncorr\n");
				return 1;
			} else {
				ncorr_set = 1;
			}
			break;

			case 7 :
			if ( sscanf(optarg, "%i", &stop_after) != 1 ) {
				ERROR("Invalid value for --stop-after\n");
				return 1;
			}
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

	if ( w_sym_str == NULL ) {
		w_sym = NULL;
		amb = NULL;
	} else {
		w_sym = get_pointgroup(w_sym_str);
		free(w_sym_str);
		if ( w_sym == NULL ) return 1;
		amb = get_ambiguities(w_sym, s_sym);
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
			                                       rmin, rmax,
			                                       amb);
			cell_free(cell);
			n_crystals++;
			reflist_free(list);

			if ( stop_after && (n_crystals == stop_after) ) break;

		}

		fprintf(stderr, "Loaded %i crystals from %i chunks\r",
		        n_crystals, ++n_chunks);

		if ( stop_after && (n_crystals == stop_after) ) break;

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

	if ( fg_graph_fn != NULL ) {
		fgfh = fopen(fg_graph_fn, "w");
		if ( fgfh == NULL ) {
			ERROR("Failed to open '%s'\n", fg_graph_fn);
		}
	}

	if ( !ncorr_set || (ncorr > n_crystals) ) {
		ncorr = n_crystals;
	}

	ccs = calc_ccs(crystals, n_crystals, ncorr, amb, rng, &mean_nac,
	               n_threads);
	if ( ccs == NULL ) {
		ERROR("Failed to allocate CCs\n");
		return 1;
	}
	STATUS("Mean number of correlations per crystal: %.1f\n", mean_nac);

	/* FIXME: Free crystals */

	for ( i=0; i<n_iter; i++ ) {
		detwin(ccs, n_crystals, assignments, fgfh, crystals);
	}

	if ( fgfh != NULL ) {
		fclose(fgfh);
	}

	if ( end_ass_fn != NULL ) {
		FILE *fh = fopen(end_ass_fn, "w");
		if ( fh == NULL ) {
			ERROR("Failed to open '%s'\n", end_ass_fn);
		} else {
			for ( i=0; i<n_crystals; i++ ) {
				fprintf(fh, "%i\n", assignments[i]);
			}
		}
		fclose(fh);
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
