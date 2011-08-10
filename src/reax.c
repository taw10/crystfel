/*
 * reax.c
 *
 * A new auto-indexer
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
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

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "cell.h"
#include "index.h"
#include "index-priv.h"


struct dvec
{
	double x;
	double y;
	double z;
	double th;
	double ph;
};


struct reax_private
{
	IndexingPrivate base;
	struct dvec *directions;
	int n_dir;
	double angular_inc;
	double *fft_in;
	fftw_complex *fft_out;
	fftw_plan plan;
	int nel;
};


static double check_dir(struct dvec *dir, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        int smin, int smax)
{
	int n, i;
	double tot;

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

		val = f->rx*dir->x + f->ry*dir->y + f->rz*dir->z;

		idx = nel/2 + nel*val/(2.0*pmax);
		fft_in[idx]++;

	}

	fftw_execute_dft_r2c(plan, fft_in, fft_out);

	tot = 0.0;
	for ( i=smin; i<=smax; i++ ) {
		double re, im;
		re = fft_out[i][0];
		im = fft_out[i][1];
		tot += sqrt(re*re + im*im);
	}

	return tot;
}


#define idx_to_m(a) ( 1.0/((2.0*pmax/(double)(a))) )

static void walk_graph(double *x, double *y, double *z, int smin, int smax,
                       fftw_complex *fft_out, int nel, double pmax,
                       double modv_exp)
{
	int i, s;
	double max, modv;

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
	assert(s>0);

	modv = 2.0*pmax / (double)s;
	*x *= modv;  *y *= modv;  *z *= modv;
}


static void fine_search(struct reax_private *p, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        int smin, int smax, double th_cen, double ph_cen,
                        double *x, double *y, double *z,
                        double modv_exp)
{
	double fom = 0.0;
	double th, ph;
	double inc;
	struct dvec dir;

	inc = p->angular_inc / 100.0;

	for ( th=th_cen-p->angular_inc; th<=th_cen+p->angular_inc; th+=inc ) {
	for ( ph=ph_cen-p->angular_inc; ph<=ph_cen+p->angular_inc; ph+=inc ) {

		double new_fom;

		dir.x = cos(ph) * sin(th);
		dir.y = sin(ph) * sin(th);
		dir.z = cos(th);

		new_fom = check_dir(&dir, flist, nel, pmax, fft_in, fft_out,
		                    plan, smin, smax);

		if ( new_fom > fom ) {
			fom = new_fom;
			*x = dir.x;  *y = dir.y;  *z = dir.z;
		}

	}
	}

	dir.x = *x;  dir.y = *y;  dir.z = *z;
	check_dir(&dir, flist, nel, pmax, fft_in, fft_out, plan, smin, smax);
	walk_graph(x, y, z, smin, smax, fft_out, nel, pmax, modv_exp);
}


void reax_index(IndexingPrivate *pp, struct image *image, UnitCell *cell)
{
	int i;
	struct reax_private *p;
	double fom;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double mod_as, mod_bs, mod_cs;
	double als, bes, gas;
	double th, ph;
	double *fft_in;
	fftw_complex *fft_out;
	int smin, smax;
	double astmin, astmax;
	double bstmin, bstmax;
	double cstmin, cstmax;
	double pmax;
	int n;
	const double ltol = 5.0;             /* Direct space axis length
	                                      * tolerance in percent */
	const double angtol = deg2rad(1.5);  /* Reciprocal angle tolerance
	                                      * in radians */

	assert(pp->indm == INDEXING_REAX);
	p = (struct reax_private *)pp;

	fft_in = fftw_malloc(p->nel*sizeof(double));
	fft_out = fftw_malloc((p->nel/2 + 1)*sizeof(fftw_complex));

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	mod_as = modulus(asx, asy, asz);
	astmin = mod_as * (1.0-ltol/100.0);
	astmax = mod_as * (1.0+ltol/100.0);

	mod_bs = modulus(bsx, bsy, bsz);
	bstmin = mod_bs * (1.0-ltol/100.0);
	bstmax = mod_bs * (1.0+ltol/100.0);

	mod_cs = modulus(csx, csy, csz);
	cstmin = mod_cs * (1.0-ltol/100.0);
	cstmax = mod_cs * (1.0+ltol/100.0);

	als = angle_between(bsx, bsy, bsz, csx, csy, csz);
	bes = angle_between(asx, asy, asz, csx, csy, csz);
	gas = angle_between(asx, asy, asz, bsx, bsy, bsz);

	pmax = 0.0;
	n = image_feature_count(image->features);
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		double val;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		val = modulus(f->rx, f->ry, f->rz);
		if ( val > pmax ) pmax = val;

	}

	smin = 2.0*pmax / astmax;
	smax = 2.0*pmax / astmin;

	/* Search for a* */
	fom = 0.0;  th = 0.0;  ph = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom;

		new_fom = check_dir(&p->directions[i], image->features,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    smin, smax);
		if ( new_fom > fom ) {
			fom = new_fom;
			th = p->directions[i].th;
			ph = p->directions[i].ph;
		}

	}
	fine_search(p, image->features, p->nel, pmax, fft_in, fft_out,
	            p->plan, smin, smax, th, ph, &asx, &asy, &asz, mod_as);

	/* Search for b* */
	smin = 2.0*pmax / bstmax;
	smax = 2.0*pmax / bstmin;
	fom = 0.0;  th = 0.0;  ph = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom, ang;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, asx, asy, asz);
		if ( fabs(ang-gas) > angtol ) continue;

		new_fom = check_dir(&p->directions[i], image->features,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    smin, smax);
		if ( new_fom > fom ) {
			fom = new_fom;
			th = p->directions[i].th;
			ph = p->directions[i].ph;
		}

	}
	fine_search(p, image->features, p->nel, pmax, fft_in, fft_out,
	            p->plan, smin, smax, th, ph, &bsx, &bsy, &bsz, mod_bs);

	/* Search for c* */
	smin = 2.0*pmax / cstmax;
	smax = 2.0*pmax / cstmin;
	fom = 0.0;  th = 0.0;  ph = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom, ang;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, asx, asy, asz);
		if ( fabs(ang-bes) > angtol ) continue;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, bsx, bsy, bsz);
		if ( fabs(ang-als) > angtol ) continue;

		new_fom = check_dir(&p->directions[i], image->features,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    smin, smax);
		if ( new_fom > fom ) {
			fom = new_fom;
			th = p->directions[i].th;
			ph = p->directions[i].ph;
		}

	}
	fine_search(p, image->features, p->nel, pmax, fft_in, fft_out,
	            p->plan, smin, smax, th, ph, &csx, &csy, &csz, mod_cs);

	image->indexed_cell = cell_new();
	cell_set_reciprocal(image->indexed_cell, asx, asy, asz,
	                    bsx, bsy, bsz, csx, csy, csz);

	fftw_free(fft_in);
	fftw_free(fft_out);
}


IndexingPrivate *reax_prepare()
{
	struct reax_private *p;
	int samp;
	double th;

	p = calloc(1, sizeof(*p));
	if ( p == NULL ) return NULL;

	p->base.indm = INDEXING_REAX;

	p->angular_inc = deg2rad(1.7);  /* From Steller (1997) */

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


	return (IndexingPrivate *)p;
}


void reax_cleanup(IndexingPrivate *pp)
{
	struct reax_private *p;

	assert(pp->indm == INDEXING_REAX);
	p = (struct reax_private *)pp;

	fftw_destroy_plan(p->plan);
	fftw_free(p->fft_in);
	fftw_free(p->fft_out);


	free(p);
}
