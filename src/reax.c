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
};


struct reax_private
{
	IndexingPrivate base;
	struct dvec *directions;
	int n_dir;
	double angular_inc;
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

	fftw_execute(plan);

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
	int i, s, mult;
	double max, modv;

	FILE *fh = fopen("fft.dat", "w");
	for ( i=0; i<nel/2+1; i++ ) {
		double re, im, m;
		re = fft_out[i][0];
		im = fft_out[i][1];
		m = sqrt(re*re + im*im);
		fprintf(fh, "%i %f\n", i, m);
	}
	fclose(fh);

	//*x *= modv_exp;  *y *= modv_exp;  *z *= modv_exp;
	//return;

	mult = 1;
	do {

		double new_s;

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

		//STATUS("Estimated axis length:%.5f nm\n",
		//       (idx_to_m(s)/mult)*1e9);

		new_s = (double)s*(mult+1)/mult;
		smin = new_s - 1;
		smax = new_s + 1;
		mult++;

	} while ( mult<5 );

	modv = 2.0*pmax / (double)s;
	modv *= mult-1;
	*x *= modv;  *y *= modv;  *z *= modv;

	//exit(1);
}


static void fine_search(struct reax_private *p, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        int smin, int smax,
                        double dx, double dy, double dz,
                        double *x, double *y, double *z,
                        double modv_exp)
{
	const int fine_samp = 100;
	double fom = 0.0;
	signed int ui, vi;
	struct dvec dir;

	for ( ui=-fine_samp; ui<fine_samp; ui++ ) {
	for ( vi=-fine_samp; vi<fine_samp; vi++ ) {

		double u, v;
		double new_fom;
		double tx, ty, tz;

		u = (double)ui/fine_samp;
		v = (double)vi/fine_samp;

		u *= p->angular_inc/fine_samp;
		v *= p->angular_inc/fine_samp;

		tx = dx*cos(u) - dy*sin(u);
		ty = dx*sin(u) + dy*cos(u);
		tz = dz;
		dir.x = tx;
		dir.y = tz*sin(v) + ty*cos(v);
		dir.z = tz*cos(v) - ty*sin(v);

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
	double dx, dy, dz;
	int nel, n;
	double pmax;
	double *fft_in;
	fftw_complex *fft_out;
	fftw_plan plan;
	int smin, smax;
	double astmin, astmax;
	double bstmin, bstmax;
	double cstmin, cstmax;
	const double cellmax = 50.0e-9;  /* 50 nm max cell size */
	const double ltol = 5.0;             /* Direct space axis length
	                                      * tolerance in percent */
	const double angtol = deg2rad(1.5);  /* Reciprocal angle tolerance
	                                      * in radians */

	assert(pp->indm == INDEXING_REAX);
	p = (struct reax_private *)pp;

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

	n = image_feature_count(image->features);
	pmax = 0.0;
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		double val;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		val = modulus(f->rx, f->ry, f->rz);

		if ( val > pmax ) pmax = val;

	}
	nel = 2.0*pmax*5.0*cellmax;

	smin = 2.0*pmax / astmax;
	smax = 2.0*pmax / astmin;

	/* Take the smallest power of 2 greater than "nel"
	 * to make the FFT go fast */
	nel = pow(2.0, floor((log(nel)/log(2.0))) + 1.0);

	fft_in = fftw_malloc(nel*sizeof(double));
	fft_out = fftw_malloc((nel/2 + 1)*sizeof(fftw_complex));

	plan = fftw_plan_dft_r2c_1d(nel, fft_in, fft_out, FFTW_ESTIMATE);

	/* Search for a* */
	fom = 0.0;  dx = 0.0;  dy = 0.0;  dz = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom;

		new_fom = check_dir(&p->directions[i], image->features,
		                    nel, pmax, fft_in, fft_out, plan,
		                    smin, smax);
		if ( new_fom > fom ) {
			fom = new_fom;
			dx = p->directions[i].x;
			dy = p->directions[i].y;
			dz = p->directions[i].z;
		}

	}
	fine_search(p, image->features, nel, pmax, fft_in, fft_out, plan,
	            smin, smax, dx, dy, dz, &asx, &asy, &asz, mod_as);

	/* Search for b* */
	smin = 2.0*pmax / bstmax;
	smax = 2.0*pmax / bstmin;
	fom = 0.0;  dx = 0.0;  dy = 0.0;  dz = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom, ang;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, asx, asy, asz);
		if ( fabs(ang-gas) > angtol ) continue;

		new_fom = check_dir(&p->directions[i], image->features,
		                    nel, pmax, fft_in, fft_out, plan,
		                    smin, smax);
		if ( new_fom > fom ) {
			fom = new_fom;
			dx = p->directions[i].x;
			dy = p->directions[i].y;
			dz = p->directions[i].z;
		}

	}
	fine_search(p, image->features, nel, pmax, fft_in, fft_out, plan,
	            smin, smax, dx, dy, dz, &bsx, &bsy, &bsz, mod_bs);

	/* Search for c* */
	smin = 2.0*pmax / cstmax;
	smax = 2.0*pmax / cstmin;
	fom = 0.0;  dx = 0.0;  dy = 0.0;  dz = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom, ang;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, asx, asy, asz);
		if ( fabs(ang-bes) > angtol ) continue;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, bsx, bsy, bsz);
		if ( fabs(ang-als) > angtol ) continue;

		new_fom = check_dir(&p->directions[i], image->features,
		                    nel, pmax, fft_in, fft_out, plan,
		                    smin, smax);
		if ( new_fom > fom ) {
			fom = new_fom;
			dx = p->directions[i].x;
			dy = p->directions[i].y;
			dz = p->directions[i].z;
		}

	}
	fine_search(p, image->features, nel, pmax, fft_in, fft_out, plan,
	            smin, smax, dx, dy, dz, &csx, &csy, &csz, mod_cs);

	fftw_destroy_plan(plan);
	fftw_free(fft_in);
	fftw_free(fft_out);

	image->indexed_cell = cell_new();
	cell_set_reciprocal(image->indexed_cell, asx, asy, asz,
	                    bsx, bsy, bsz, csx, csy, csz);
}


IndexingPrivate *reax_prepare()
{
	struct reax_private *priv;
	int ui, vi;
	int samp;

	priv = calloc(1, sizeof(*priv));
	if ( priv == NULL ) return NULL;

	priv->base.indm = INDEXING_REAX;

	/* Decide on sampling interval */
	priv->angular_inc = 0.03;  /* From Steller (1997) */
	samp = (2.0 * M_PI) / priv->angular_inc;

	priv->n_dir = samp*samp;
	priv->directions = malloc(priv->n_dir*sizeof(struct dvec));
	if ( priv == NULL) {
		free(priv);
		return NULL;
	}

	/* Generate vectors for 1D Fourier transforms */
	for ( ui=0; ui<samp; ui++ ) {
	for ( vi=0; vi<samp; vi++ ) {

		double u, v;
		double th, ph;
		struct dvec *dir;

		u = (double)ui/samp;
		v = (double)vi/samp;

		/* Uniform sampling of a hemisphere */
		th = 2.0 * M_PI * u;  /* "Longitude" */
		ph = acos(v);         /* "Latitude" */

		dir = &priv->directions[ui + vi*samp];

		dir->x = cos(th) * sin(ph);
		dir->y = sin(th) * sin(ph);
		dir->z = cos(ph);

	}
	}

	return (IndexingPrivate *)priv;
}
