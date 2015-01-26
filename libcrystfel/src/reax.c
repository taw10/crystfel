/*
 * reax.c
 *
 * A new auto-indexer
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2014 Takanori Nakane
 *
 * Authors:
 *   2011-2014 Thomas White <taw@physics.org>
 *   2014 Takanori Nakane <nakane.t@gmail.com>
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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <fftw3.h>
#include <fenv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell.h"
#include "cell-utils.h"
#include "index.h"
#include "reax.h"


/* Minimum number of standard deviations above the mean a peak must be in the
 * 1D FT to qualify as a candidate vector */
#define MIN_SIGMAS (3.0)


/* Maximum number of times the angular tolerance that vectors are permitted to
 * be together before they get merged by squash_vectors() */
#define INC_TOL_MULTIPLIER (3.0)


/* Maximum number of candidate vectors to find (we will take the best ones) */
#define MAX_CANDIDATES (1024)


/* Choose the best solution from this many candidate cells */
#define MAX_REAX_CELL_CANDIDATES (32)

struct dvec
{
	double x;
	double y;
	double z;
	double th;
	double ph;
};


struct reax_candidate
{
	struct dvec v;   /* This is the vector for the candidate */
	double fom;
};


struct reax_search_v
{
	unsigned int smin;
	unsigned int smax;  /* Search for vector in this range */

	struct reax_candidate *cand;  /* Candidate vectors go here */
	int n_cand;                   /* There are this many candidates */
	int max_warned;
};


struct reax_search
{
	struct reax_search_v *search;  /* Search for these vectors */
	int n_search;                  /* There are this many vectors to find */
	double pmax;                   /* The maximum feature resolution */
};


struct reax_private
{
	IndexingMethod indm;
	struct dvec *directions;
	int n_dir;
	double angular_inc;
	UnitCell *cell;

	double *fft_in;
	fftw_complex *fft_out;
	fftw_plan plan;
	int nel;

	fftw_complex *r_fft_in;
	fftw_complex *r_fft_out;
	fftw_plan r_plan;
	int ch;
	int cw;
};


static void fill_and_transform(struct dvec *dir, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        const struct rigid_group *rg, struct detector *det)
{
	int n, i;

	for ( i=0; i<nel; i++ ) {
		fft_in[i] = 0.0;
	}

	n = image_feature_count(flist);
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		double val;
		int idx;

		f = image_get_feature(flist, i);
		if ( f == NULL ) continue;

		if ( rg != NULL ) {

			struct panel *p;

			p = find_panel(det, f->fs, f->ss);
			assert(p != NULL);

			if ( !panel_is_in_rigid_group(rg, p) ) continue;

		}

		val = f->rx*dir->x + f->ry*dir->y + f->rz*dir->z;

		idx = nel/2 + nel*val/(2.0*pmax);
		fft_in[idx]++;

	}

	fftw_execute_dft_r2c(plan, fft_in, fft_out);
}


static void add_candidate(struct reax_search_v *s, struct reax_candidate *c)
{
	int idx;

	if ( s->n_cand == MAX_CANDIDATES ) {
		if ( !s->max_warned ) {
			ERROR("WARNING: Too many candidates.\n");
			s->max_warned = 1;
		}
		return;
	}

	idx = s->n_cand++;

	s->cand[idx].v = c->v;
	s->cand[idx].fom = c->fom;
}


static double check_dir(struct dvec *dir, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        struct reax_search *s,
                        const struct rigid_group *rg, struct detector *det)
{
	int i;
	double tot;

	fill_and_transform(dir, flist, nel, pmax, fft_in, fft_out,
	                   plan, rg, det);

	tot = 0.0;
	for ( i=0; i<s->n_search; i++ ) {

		double tot = 0.0;
		double peak = 0.0;
		double peak_mod = 0.0;
		double mean;
		double sd = 0.0;
		int j;
		int n = 0;

		for ( j=0; j<nel/2+1; j++ ) {

			double re, im, am;

			re = fft_out[j][0];
			im = fft_out[j][1];
			am = sqrt(re*re + im*im);

			tot += am;
			n++;

			if ( ( j >= s->search[i].smin )
			  && ( j <= s->search[i].smax ) ) {
				if ( am > peak ) {
					peak = am;
					peak_mod = (double)j/(2.0*pmax);
				}
			}

		}
		mean = tot/(double)n;

		for ( j=0; j<nel/2+1; j++ ) {

			double re, im, am;

			re = fft_out[j][0];
			im = fft_out[j][1];
			am = sqrt(re*re + im*im);

			sd += pow(am - mean, 2.0);

		}
		sd = sqrt(sd/(double)n);

		/* If sufficiently strong, add to list of candidates */
		if ( peak > mean+MIN_SIGMAS*sd ) {

			struct reax_candidate c;

			c.v.x = dir->x * peak_mod;
			c.v.y = dir->y * peak_mod;
			c.v.z = dir->z * peak_mod;
			c.fom = peak;

			add_candidate(&s->search[i],  &c);

		}

	}

	return tot;
}


/* Refine a direct space vector.  From Clegg (1984)
 * with added iteration because more reflections might get included as the
 * refinement proceeds. */
static double iterate_refine_vector(double *x, double *y, double *z,
                                    ImageFeatureList *flist)
{
	int fi, n, err;
	gsl_matrix *C;
	gsl_vector *A;
	gsl_vector *t;
	gsl_matrix *s_vec;
	gsl_vector *s_val;
	double dtmax;

	A = gsl_vector_calloc(3);
	C = gsl_matrix_calloc(3, 3);

	n = image_feature_count(flist);
	fesetround(1);
	for ( fi=0; fi<n; fi++ ) {

		struct imagefeature *f;
		double val;
		double kn, kno;
		double xv[3];
		int i, j;

		f = image_get_feature(flist, fi);
		if ( f == NULL ) continue;

		kno = f->rx*(*x) + f->ry*(*y) + f->rz*(*z);  /* Sorry ... */
		kn = nearbyint(kno);
		if ( fabs(kn - kno) > 0.3 ) continue;

		xv[0] = f->rx;  xv[1] = f->ry;  xv[2] = f->rz;

		for ( i=0; i<3; i++ ) {

			val = gsl_vector_get(A, i);
			gsl_vector_set(A, i, val+xv[i]*kn);

			for ( j=0; j<3; j++ ) {
				val = gsl_matrix_get(C, i, j);
				gsl_matrix_set(C, i, j, val+xv[i]*xv[j]);
			}

		}

	}

	s_val = gsl_vector_calloc(3);
	s_vec = gsl_matrix_calloc(3, 3);
	err = gsl_linalg_SV_decomp_jacobi(C, s_vec, s_val);
	if ( err ) {
		ERROR("SVD failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		gsl_matrix_free(C);
		gsl_vector_free(A);
		return 0.0;
	}

	t = gsl_vector_calloc(3);
	err = gsl_linalg_SV_solve(C, s_vec, s_val, A, t);
	if ( err ) {
		ERROR("Matrix solution failed: %s\n", gsl_strerror(err));
		gsl_matrix_free(s_vec);
		gsl_vector_free(s_val);
		gsl_matrix_free(C);
		gsl_vector_free(A);
		gsl_vector_free(t);
		return 0.0;
	}

	gsl_matrix_free(s_vec);
	gsl_vector_free(s_val);

	dtmax  = fabs(*x - gsl_vector_get(t, 0));
	dtmax += fabs(*y - gsl_vector_get(t, 1));
	dtmax += fabs(*z - gsl_vector_get(t, 2));

	*x = gsl_vector_get(t, 0);
	*y = gsl_vector_get(t, 1);
	*z = gsl_vector_get(t, 2);

	gsl_matrix_free(C);
	gsl_vector_free(A);
	gsl_vector_free(t);

	return dtmax;
}


static void refine_vector(ImageFeatureList *flist, struct dvec *dir)
{
	int i;
	double sm;

	i = 0;
	do {

		sm  = iterate_refine_vector(&dir->x, &dir->y, &dir->z, flist);
		i++;

	} while ( (sm > 0.001e-9) && (i<10) );
}


static void squash_vectors(struct reax_search *s, double tol)
{
	int i;

	for ( i=0; i<s->n_search; i++ ) {

		struct reax_search_v *sv;
		struct reax_candidate *new;
		int j, k;
		int n_invalid = 0;
		int n_copied;

		sv = &s->search[i];

		for ( j=0; j<sv->n_cand; j++ ) {
		for ( k=0; k<sv->n_cand; k++ ) {

			struct reax_candidate *v1, *v2;

			if ( j == k ) continue;

			v1 = &sv->cand[j];
			v2 = &sv->cand[k];

			if ( angle_between(v1->v.x, v1->v.y, v1->v.z,
			                   v2->v.x, v2->v.y, v2->v.z) < tol )
			{
				if ( !isnan(v1->fom) && !isnan(v2->fom ) ) {
					if ( v1->fom > v2->fom ) {
						v2->fom = NAN;
					} else {
						v1->fom = NAN;
					}
					n_invalid++;
				}
			}

		}
		}

		new = calloc(sv->n_cand - n_invalid,
		             sizeof(struct reax_candidate));
		if ( new == NULL ) {
			ERROR("Failed to allocate memory for squashed"
			      " candidate list.\n");
			return;
		}

		n_copied = 0;
		for ( j=0; j<sv->n_cand; j++ ) {
			if ( !isnan(sv->cand[j].fom) ) {

				new[n_copied] = sv->cand[j];
				n_copied++;

			}
		}
		assert(sv->n_cand - n_invalid == n_copied);

		free(sv->cand);
		//STATUS("Search vector %i:", i);
		//STATUS(" squashed %i candidates down to %i\n",
		//       sv->n_cand, n_copied);
		sv->n_cand = n_copied;
		sv->cand = new;

	}
}


static UNUSED void show_vectors(struct reax_search *s, const char *pre)
{
	int i;

	/* For each direction being searched for */
	for ( i=0; i<s->n_search; i++ ) {

		int j;

		for ( j=0; j<s->search[i].n_cand; j++ ) {
			STATUS("%s: %i/%i: %+6.2f %+6.2f %+6.2f nm %.2f\n", pre,
			       i, j, s->search[i].cand[j].v.x*1e9,
			       s->search[i].cand[j].v.y*1e9,
			       s->search[i].cand[j].v.z*1e9,
			       s->search[i].cand[j].fom);
		}

	}
}


static void find_candidates(struct reax_private *p,
                            ImageFeatureList *flist, double pmax,
                            double *fft_in, fftw_complex *fft_out,
                            struct reax_search *s,
                            const char *rg, struct detector *det)
{
	int i;

	for ( i=0; i<s->n_search; i++ ) {
		s->search[i].cand = calloc(MAX_CANDIDATES,
		                           sizeof(struct reax_candidate));
		s->search[i].n_cand = 0;
	}

	for ( i=0; i<p->n_dir; i++ ) {
		check_dir(&p->directions[i], flist,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    s, NULL, NULL);
	}

	squash_vectors(s, INC_TOL_MULTIPLIER*p->angular_inc);

	//show_vectors(s, "BEFORE");

	for ( i=0; i<s->n_search; i++ ) {

		struct reax_search_v *sv;
		int j;

		sv = &s->search[i];

		for ( j=0; j<sv->n_cand; j++ ) {
			refine_vector(flist, &sv->cand[j].v);
		}

	}

	//show_vectors(s, "FINAL");
}


/* Set up search parameters to look for all three cell axes */
static struct reax_search *search_all_axes(UnitCell *cell, double pmax)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double mod_a, mod_b, mod_c;
	double amin, amax;
	double bmin, bmax;
	double cmin, cmax;
	unsigned int smin, smax;
	const double ltol = 10.0; /* Direct space axis length tolerance in % */
	struct reax_search *s;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	mod_a = modulus(ax, ay, az);
	amin = mod_a * (1.0-ltol/100.0);
	amax = mod_a * (1.0+ltol/100.0);

	mod_b = modulus(bx, by, bz);
	bmin = mod_b * (1.0-ltol/100.0);
	bmax = mod_b * (1.0+ltol/100.0);

	mod_c = modulus(cx, cy, cz);
	cmin = mod_c * (1.0-ltol/100.0);
	cmax = mod_c * (1.0+ltol/100.0);

	s = malloc(3*sizeof(*s));
	s->pmax = pmax;
	s->n_search = 3;
	s->search = malloc(3*sizeof(struct reax_search_v));
	smin = 2.0*pmax * amin;  smax = 2.0*pmax * amax;
	s->search[0].smin = smin;  s->search[0].smax = smax;
	s->search[0].max_warned = 0;
	smin = 2.0*pmax * bmin;  smax = 2.0*pmax * bmax;
	s->search[1].smin = smin;  s->search[1].smax = smax;
	s->search[1].max_warned = 0;
	smin = 2.0*pmax * cmin;  smax = 2.0*pmax * cmax;
	s->search[2].smin = smin;  s->search[2].smax = smax;
	s->search[2].max_warned = 0;

	return s;
}


static double get_model_phase(double x, double y, double z, ImageFeatureList *f,
                              int nel, double pmax, double *fft_in,
                              fftw_complex *fft_out, fftw_plan plan,
                              int smin, int smax, const struct rigid_group *rg,
                              struct detector *det)
{
	struct dvec dir;
	int s, i;
	double max;
	double re, im;

	dir.x = x;  dir.y = y;  dir.z = z;

	fill_and_transform(&dir, f, nel, pmax, fft_in, fft_out, plan, rg, det);

	s = -1;
	max = 0.0;
	for ( i=smin; i<=smax; i++ ) {

		double re, im, m;

		re = fft_out[i][0];
		im = fft_out[i][1];
		m = sqrt(re*re + im*im);
		if ( m > max ) {
			max = m;
			s = i;
		}

	}

	re = fft_out[s][0];
	im = fft_out[s][1];

	return atan2(im, re);
}


static void refine_rigid_group(struct image *image, UnitCell *cell,
                               const struct rigid_group *rg, double pmax,
                               double *fft_in, fftw_complex *fft_out,
                               fftw_plan plan, int smin, int smax,
                               struct detector *det, struct reax_private *pr)
{
	double ax, ay, az, ma;
	double bx, by, bz, mb;
	double cx, cy, cz, mc;
	double pha, phb, phc;
	int i, j;
	fftw_complex *r_fft_in;
	fftw_complex *r_fft_out;
	double m2m;
	signed int aix, aiy;
	signed int bix, biy;
	signed int cix, ciy;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	ma = modulus(ax, ay, az);
	mb = modulus(bx, by, bz);
	mc = modulus(cx, cy, cz);

	pha = get_model_phase(ax/ma, ay/ma, az/ma, image->features,
	                      pr->nel, pmax, fft_in, fft_out, plan,
	                      smin, smax, rg, det);
	phb = get_model_phase(bx/mb, by/mb, bz/mb, image->features,
	                      pr->nel, pmax, fft_in, fft_out, plan,
	                      smin, smax, rg, det);
	phc = get_model_phase(cx/mc, cy/mc, cz/mc, image->features,
	                      pr->nel, pmax, fft_in, fft_out, plan,
	                      smin, smax, rg, det);

	r_fft_in = fftw_malloc(pr->cw*pr->ch*sizeof(fftw_complex));
	r_fft_out = fftw_malloc(pr->cw*pr->ch*sizeof(fftw_complex));
	for ( i=0; i<pr->cw; i++ ) {
	for ( j=0; j<pr->ch; j++ ) {
		r_fft_in[i+pr->cw*j][0] = 0.0;
		r_fft_in[i+pr->cw*j][1] = 0.0;
	}
	}

	ma = modulus(ax, ay, 0.0);
	mb = modulus(bx, by, 0.0);
	mc = modulus(cx, cy, 0.0);
	m2m = ma;
	if ( mb > m2m ) m2m = mb;
	if ( mc > m2m ) m2m = mc;

	aix = (pr->cw/2)*ax/m2m;  aiy = (pr->ch/2)*ay/m2m;
	bix = (pr->cw/2)*bx/m2m;  biy = (pr->ch/2)*by/m2m;
	cix = (pr->cw/2)*cx/m2m;  ciy = (pr->ch/2)*cy/m2m;

	if ( aix < 0 ) aix += pr->cw/2;
	if ( bix < 0 ) bix += pr->cw/2;
	if ( cix < 0 ) cix += pr->cw/2;

	if ( aiy < 0 ) aiy += pr->ch/2;
	if ( biy < 0 ) biy += pr->ch/2;
	if ( ciy < 0 ) ciy += pr->ch/2;

	r_fft_in[aix + pr->cw*aiy][0] = cos(pha);
	r_fft_in[aix + pr->cw*aiy][1] = sin(pha);
	r_fft_in[pr->cw-aix + pr->cw*(pr->ch-aiy)][0] = cos(pha);
	r_fft_in[pr->cw-aix + pr->cw*(pr->ch-aiy)][1] = -sin(pha);

	r_fft_in[bix + pr->cw*biy][0] = cos(phb);
	r_fft_in[bix + pr->cw*biy][1] = sin(phb);
	r_fft_in[pr->cw-bix + pr->cw*(pr->ch-biy)][0] = cos(phb);
	r_fft_in[pr->cw-bix + pr->cw*(pr->ch-biy)][1] = -sin(phb);

	r_fft_in[cix + pr->cw*ciy][0] = cos(phc);
	r_fft_in[cix + pr->cw*ciy][1] = sin(phc);
	r_fft_in[pr->cw-cix + pr->cw*(pr->ch-ciy)][0] = cos(phc);
	r_fft_in[pr->cw-cix + pr->cw*(pr->ch-ciy)][1] = -sin(phc);

	const int tidx = 1;
	r_fft_in[tidx][0] = 1.0;
	r_fft_in[tidx][1] = 0.0;

//	STATUS("%i %i\n", aix, aiy);
//	STATUS("%i %i\n", bix, biy);
//	STATUS("%i %i\n", cix, ciy);

	fftw_execute_dft(pr->r_plan, r_fft_in, r_fft_out);

//	max = 0.0;
//	FILE *fh = fopen("centering.dat", "w");
//	for ( i=0; i<pr->cw; i++ ) {
//	for ( j=0; j<pr->ch; j++ ) {
//
//		double re, im, am, ph;
//
//		re = r_fft_out[i + pr->cw*j][0];
//		im = r_fft_out[i + pr->cw*j][1];
//		am = sqrt(re*re + im*im);
//		ph = atan2(im, re);
//
//		if ( am > max ) {
//			max = am;
//			max_i = i;
//			max_j = j;
//		}
//
//		fprintf(fh, "%f ", am);
//
//	}
//	fprintf(fh, "\n");
//	}
//	STATUS("Max at %i, %i\n", max_i, max_j);
//	fclose(fh);
//	exit(1);

//	STATUS("Offsets for '%s': %.2f, %.2f pixels\n", rg, dx, dy);
}


static UNUSED void refine_all_rigid_groups(struct image *image, UnitCell *cell,
                                           double pmax,
                                           double *fft_in,
                                           fftw_complex *fft_out,
                                           fftw_plan plan, int smin, int smax,
                                           struct detector *det,
                                           struct reax_private *p)
{
	int i;

	for ( i=0; i<image->det->n_rigid_groups; i++ ) {
		refine_rigid_group(image, cell, image->det->rigid_groups[i],
		                   pmax, fft_in, fft_out, plan, smin, smax,
		                   det, p);
	}
}


static double max_feature_resolution(ImageFeatureList *flist)
{
	double pmax;
	int i, n;

	pmax = 0.0;
	n = image_feature_count(flist);
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		double val;

		f = image_get_feature(flist, i);
		if ( f == NULL ) continue;

		val = modulus(f->rx, f->ry, f->rz);
		if ( val > pmax ) pmax = val;

	}

	return pmax;
}


static int right_handed_vec(struct rvec a, struct rvec b, struct rvec c)
{
	struct rvec aCb;
	double aCb_dot_c;

	/* "a" cross "b" */
	aCb.u = a.v*b.w - a.w*b.v;
	aCb.v = - (a.u*b.w - a.w*b.u);
	aCb.w = a.u*b.v - a.v*b.u;

	/* "a cross b" dot "c" */
	aCb_dot_c = aCb.u*c.u + aCb.v*c.v + aCb.w*c.w;

	if ( aCb_dot_c > 0.0 ) return 1;
	return 0;
}


struct cell_candidate
{
	UnitCell *cell;
	double fom;
};


struct cell_candidate_list
{
	struct cell_candidate *cand;
	int n_cand;
};


static int check_twinning(UnitCell *c1, UnitCell *c2, int verbose)
{
	int i;
	int n_dup;
	const int n_trials = 40;

	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;

	gsl_rng *rng;

	/* This is a rubbish RNG, but it serves for this purpose: nothing more
	 * than "I couldn't be bothered to think of n_trials sets of random
	 * indices". */
	rng = gsl_rng_alloc(gsl_rng_taus2);

	cell_get_reciprocal(c1, &asx, &asy, &asz,
	                        &bsx, &bsy, &bsz,
	                        &csx, &csy, &csz);
	cell_get_cartesian(c2, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	n_dup = 0;
	for ( i=0; i<n_trials; i++ ) {

		signed int h, k, l;
		double h2, k2, l2;
		double rx, ry, rz;
		double dev;
		signed int h2i, k2i, l2i;

		h = gsl_rng_uniform_int(rng, 10);
		k = gsl_rng_uniform_int(rng, 10);
		l = gsl_rng_uniform_int(rng, 10);

		/* Position of this (randomly selected)
		 * reciprocal lattice point */
		rx = h*asx + k*bsx + l*csx;
		ry = h*asy + k*bsy + l*csy;
		rz = h*asz + k*bsz + l*csz;

		/* Indices of this point in the basis of the other cell */
		h2 = rx*ax + ry*ay + rz*az;
		k2 = rx*bx + ry*by + rz*bz;
		l2 = rx*cx + ry*cy + rz*cz;

		h2i = lrint(h2);
		k2i = lrint(k2);
		l2i = lrint(l2);

		dev = pow(h2i-h2, 2.0) + pow(k2i-k2, 2.0) + pow(l2i-l2, 2.0);

		if ( verbose ) {
			STATUS("%3i %3i %3i -> %5.2f %5.2f %5.2f -> "
			       "%3i %3i %3i -> %5.2f\n", h, k, l,
			       h2, k2, l2, h2i, k2i, l2i, dev);
		}

		if ( dev < 0.1 ) {
			n_dup++;
		}

	}

	if ( verbose ) {
		STATUS("%i duplicates.\n", n_dup);
	}

	gsl_rng_free(rng);

	if ( n_dup > 10 ) return 1;
	return 0;
}


/* Return true if "cnew" accounts for more than 25% of the peaks predicted by
 * any of the "ncells" cells in "cells". */
static int twinned(UnitCell *cnew, struct cell_candidate_list *cl)
{
	int i;

	for ( i=0; i<cl->n_cand; i++ ) {
		if ( check_twinning(cnew, cl->cand[i].cell, 0) ) return 1;
	}

	return 0;
}


static int check_vector_combination(struct dvec *vi, struct dvec *vj,
                                    struct dvec *vk, UnitCell *cell)
{
	double ang;
	double a, b, c, al, be, ga;
	const double angtol = deg2rad(5.0);

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	ang = angle_between(vi->x, vi->y, vi->z, vj->x, vj->y, vj->z);
	if ( fabs(ang-ga) > angtol ) return 0;

	ang = angle_between(vi->x, vi->y, vi->z, vk->x, vk->y, vk->z);
	if ( fabs(ang-be) > angtol ) return 0;

	ang = angle_between(vj->x, vj->y, vj->z, vk->x, vk->y, vk->z);
	if ( fabs(ang-al) > angtol ) return 0;

	return 1;
}


static void add_cell_candidate(struct cell_candidate_list *cl, UnitCell *cnew,
                               double fom)
{
	struct cell_candidate cshift;
	int i, cpos;

	cpos = cl->n_cand;
	for ( i=0; i<cl->n_cand; i++ ) {
		if ( fom > cl->cand[i].fom ) {
			cpos = i;
			break;
		}
	}

	cshift.cell = cnew;
	cshift.fom = fom;

	for ( i=cpos; i<cl->n_cand; i++ ) {

		struct cell_candidate cshift2;
		cshift2 = cl->cand[i];
		cl->cand[i] = cshift;
		cshift = cshift2;

	}

	if ( cl->n_cand >= MAX_REAX_CELL_CANDIDATES ) {
		/* "cshift" just fell off the end of the list */
	} else {
		cl->cand[cl->n_cand++] = cshift;
	}
}


static int check_cell(struct reax_private *dp, struct image *image,
                      UnitCell *cell)
{
	UnitCell *out;
	Crystal *cr;

	out = cell_new_from_cell(cell);

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 0;
	}

	crystal_set_cell(cr, out);

	if ( dp->indm & INDEXING_CHECK_PEAKS ) {
		if ( !peak_sanity_check(image, &cr, 1) ) {
			crystal_free(cr);  /* Frees the cell as well */
			return 0;
		}
	}

	image_add_crystal(image, cr);

	return 1;
}


static int assemble_cells_from_candidates(struct image *image,
                                           struct reax_search *s,
                                           UnitCell *cell,
                                           struct reax_private *p)
{
	int i, j, k;
	signed int ti, tj, tk;
	int rval;
	struct cell_candidate_list cl;

	cl.cand = calloc(MAX_REAX_CELL_CANDIDATES,
	                 sizeof(struct cell_candidate));
	if ( cl.cand == NULL ) {
		ERROR("Failed to allocate cell candidate list.\n");
		return 0;
	}
	cl.n_cand = 0;

	/* Find candidates for axes 0 and 1 which have the right angle */
	for ( i=0; i<s->search[0].n_cand; i++ ) {
	for ( j=0; j<s->search[1].n_cand; j++ ) {
	for ( k=0; k<s->search[2].n_cand; k++ ) {
	for ( ti=-1; ti<=1; ti+=2 ) {
	for ( tj=-1; tj<=1; tj+=2 ) {
	for ( tk=-1; tk<=1; tk+=2 ) {

		struct dvec vi, vj, vk;
		struct rvec ai, bi, ci;
		UnitCell *cnew;
		double fom;

		vi = s->search[0].cand[i].v;
		vj = s->search[1].cand[j].v;
		vk = s->search[2].cand[k].v;

		vi.x *= ti;  vi.y *= ti;  vi.z *= ti;
		vj.x *= tj;  vj.y *= tj;  vj.z *= tj;
		vk.x *= tk;  vk.y *= tk;  vk.z *= tk;

		if ( !check_vector_combination(&vi, &vj, &vk, cell) ) continue;

		ai.u = vi.x;  ai.v = vi.y;  ai.w = vi.z;
		bi.u = vj.x;  bi.v = vj.y;  bi.w = vj.z;
		ci.u = vk.x;  ci.v = vk.y;  ci.w = vk.z;

		if ( !right_handed_vec(ai, bi, ci) ) continue;

		/* We have three vectors with the right angles */
		cnew = cell_new_from_direct_axes(ai, bi, ci);

		if ( twinned(cnew, &cl) ) {
			cell_free(cnew);
			continue;
		}

		fom = 1.0;
		add_cell_candidate(&cl, cnew, fom);

	}
	}
	}
	}
	}
	}

	for ( i=0; i<cl.n_cand; i++ ) {
		double a, b, c, al, be, ga;
		double aA, bA, cA, alA, beA, gaA;
		int w = 0;
//		STATUS("%i: %f\n", i, cl.cand[i].fom);
		cell_get_parameters(cl.cand[i].cell, &a, &b, &c, &al, &be, &ga);
		cell_get_parameters(cl.cand[i].cell, &aA, &bA, &cA,
		                                     &alA, &beA, &gaA);
		if ( fabs(a - aA) > aA/10.0 ) w = 1;
		if ( fabs(b - bA) > bA/10.0 ) w = 1;
		if ( fabs(c - cA) > cA/10.0 ) w = 1;
		if ( fabs(al - alA) > deg2rad(5.0) ) w = 1;
		if ( fabs(be - beA) > deg2rad(5.0) ) w = 1;
		if ( fabs(ga - gaA) > deg2rad(5.0) ) w = 1;
		if ( w ) {
			STATUS("This cell is a long way from that sought:\n");
			cell_print(cl.cand[i].cell);
		}
	}

	rval = 0;
	for ( i=0; i<cl.n_cand; i++ ) {
		if ( check_cell(p,  image, cl.cand[i].cell) ) {
			rval = 1;
			break;
		}
	}

	free(cl.cand);
	return rval;
}


int reax_index(IndexingPrivate *pp, struct image *image)
{
	struct reax_private *p;
	double *fft_in;
	fftw_complex *fft_out;
	double pmax;
	struct reax_search *s;
	int i;
	int rval = 0;

	p = (struct reax_private *)pp;

	fft_in = fftw_malloc(p->nel*sizeof(double));
	fft_out = fftw_malloc((p->nel/2 + 1)*sizeof(fftw_complex));

	pmax = max_feature_resolution(image->features);

	/* Sanity check */
	if ( pmax < 1e4 ) {
		fftw_free(fft_in);
		fftw_free(fft_out);
		return 0;
	}

	s = search_all_axes(p->cell, pmax);
	find_candidates(p, image->features, pmax, fft_in, fft_out, s,
	                NULL, image->det);

//	refine_all_rigid_groups(image, image->candidate_cells[0], pmax,
//	                        fft_in, fft_out, p->plan, smin, smax,
//	                        image->det, p);

	rval = assemble_cells_from_candidates(image, s, p->cell, p);

	for ( i=0; i<s->n_search; i++ ) {
		free(s->search[i].cand);
	}
	free(s->search);
	free(s);
	fftw_free(fft_in);
	fftw_free(fft_out);
	return rval;
}


IndexingPrivate *reax_prepare(IndexingMethod *indm, UnitCell *cell,
                              struct detector *det, float *ltl)
{
	struct reax_private *p;
	int samp;
	double th;

	if ( cell == NULL ) {
		ERROR("ReAx needs a unit cell.\n");
		return NULL;
	}

	p = calloc(1, sizeof(*p));
	if ( p == NULL ) return NULL;

	p->cell = cell;

	/* Flags that ReAx knows about */
	*indm &= INDEXING_METHOD_MASK | INDEXING_CHECK_PEAKS;

	/* Flags that ReAx requires */
	*indm |= INDEXING_USE_LATTICE_TYPE;
	*indm |= INDEXING_USE_CELL_PARAMETERS;

	p->angular_inc = deg2rad(1.0);

	/* Reserve memory, over-estimating the number of directions */
	samp = 2.0*M_PI / p->angular_inc;
	p->directions = malloc(samp*samp*sizeof(struct dvec));
	if ( p == NULL) {
		free(p);
		return NULL;
	}
	STATUS("Allocated space for %i directions\n", samp*samp);

	/* Generate vectors for 1D Fourier transforms */
	fesetround(1);  /* Round to nearest */
	p->n_dir = 0;
	for ( th=0.0; th<M_PI_2; th+=p->angular_inc ) {

		double ph, phstep, n_phstep;

		n_phstep = 2.0*M_PI*sin(th)/p->angular_inc;
		n_phstep = nearbyint(n_phstep);
		phstep = 2.0*M_PI/n_phstep;

		for ( ph=0.0; ph<2.0*M_PI; ph+=phstep ) {

			struct dvec *dir;

			assert(p->n_dir<samp*samp);

			dir = &p->directions[p->n_dir++];

			dir->x = cos(ph) * sin(th);
			dir->y = sin(ph) * sin(th);
			dir->z = cos(th);
			dir->th = th;
			dir->ph = ph;

		}

	}
	STATUS("Generated %i directions (angular increment %.3f deg)\n",
	       p->n_dir, rad2deg(p->angular_inc));

	p->nel = 1024;

	/* These arrays are not actually used */
	p->fft_in = fftw_malloc(p->nel*sizeof(double));
	p->fft_out = fftw_malloc((p->nel/2 + 1)*sizeof(fftw_complex));

	p->plan = fftw_plan_dft_r2c_1d(p->nel, p->fft_in, p->fft_out,
	                               FFTW_MEASURE);

	p->cw = 128; p->ch = 128;

	/* Also not used */
	p->r_fft_in = fftw_malloc(p->cw*p->ch*sizeof(fftw_complex));
	p->r_fft_out = fftw_malloc(p->cw*p->ch*sizeof(fftw_complex));

	p->r_plan = fftw_plan_dft_2d(p->cw, p->ch, p->r_fft_in, p->r_fft_out,
	                             1, FFTW_MEASURE);

	p->indm = *indm;

	return (IndexingPrivate *)p;
}


void reax_cleanup(IndexingPrivate *pp)
{
	struct reax_private *p;

	p = (struct reax_private *)pp;

	free(p->directions);

	fftw_destroy_plan(p->plan);
	fftw_free(p->fft_in);
	fftw_free(p->fft_out);

	fftw_destroy_plan(p->r_plan);
	fftw_free(p->r_fft_in);
	fftw_free(p->r_fft_out);

	free(p);
}
