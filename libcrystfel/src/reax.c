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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

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

	fftw_complex *r_fft_in;
	fftw_complex *r_fft_out;
	fftw_plan r_plan;
	int ch;
	int cw;
};


static double check_dir(struct dvec *dir, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        int smin, int smax,
                        const char *rg, struct detector *det)
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

		if ( rg != NULL ) {

			struct panel *p;

			p = find_panel(det, f->fs, f->ss);
			assert(p != NULL);

			if ( p->rigid_group != rg ) continue;

		}

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


/* Refine a direct space vector.  From Clegg (1984) */
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
		if ( kn - kno > 0.3 ) continue;

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

	return dtmax;
}


static void refine_cell(struct image *image, UnitCell *cell,
                        ImageFeatureList *flist)
{
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	int i;
	double sm;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	i = 0;
	do {

		sm  = iterate_refine_vector(&ax, &ay, &az, flist);
		sm += iterate_refine_vector(&bx, &by, &bz, flist);
		sm += iterate_refine_vector(&cx, &cy, &cz, flist);
		i++;

	} while ( (sm > 0.001e-9) && (i<10) );

	cell_set_cartesian(cell, ax, ay, az, bx, by, bz, cx, cy, cz);

	if ( i == 10 ) {
		cell_free(image->indexed_cell);
		image->indexed_cell = NULL;
	}
}


static void fine_search(struct reax_private *p, ImageFeatureList *flist,
                        int nel, double pmax, double *fft_in,
                        fftw_complex *fft_out, fftw_plan plan,
                        int smin, int smax, double th_cen, double ph_cen,
                        double *x, double *y, double *z)
{
	double fom = 0.0;
	double th, ph;
	double inc;
	struct dvec dir;
	int i, s;
	double max, modv;

	inc = p->angular_inc / 100.0;

	for ( th=th_cen-p->angular_inc; th<=th_cen+p->angular_inc; th+=inc ) {
	for ( ph=ph_cen-p->angular_inc; ph<=ph_cen+p->angular_inc; ph+=inc ) {

		double new_fom;

		dir.x = cos(ph) * sin(th);
		dir.y = sin(ph) * sin(th);
		dir.z = cos(th);

		new_fom = check_dir(&dir, flist, nel, pmax, fft_in, fft_out,
		                    plan, smin, smax, NULL, NULL);

		if ( new_fom > fom ) {
			fom = new_fom;
			*x = dir.x;  *y = dir.y;  *z = dir.z;
		}

	}
	}

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

	modv = (double)s / (2.0*pmax);
	*x *= modv;  *y *= modv;  *z *= modv;
}


static double get_model_phase(double x, double y, double z, ImageFeatureList *f,
                              int nel, double pmax, double *fft_in,
                              fftw_complex *fft_out, fftw_plan plan,
                              int smin, int smax, const char *rg,
                              struct detector *det)
{
	struct dvec dir;
	int s, i;
	double max;
	double re, im;

	dir.x = x;  dir.y = y;  dir.z = z;

	check_dir(&dir, f, nel, pmax, fft_in,fft_out, plan, smin, smax,
	          rg, det);

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
                               const char *rg, int nel, double pmax,
                               double *fft_in, fftw_complex *fft_out,
                               fftw_plan plan, int smin, int smax,
                               struct detector *det, struct reax_private *pr)
{
	double ax, ay, az, ma;
	double bx, by, bz, mb;
	double cx, cy, cz, mc;
	double pha, phb, phc;
	struct panel *p;
	int i, j;
	fftw_complex *r_fft_in;
	fftw_complex *r_fft_out;
	double m2m;
	signed int aix, aiy;
	signed int bix, biy;
	signed int cix, ciy;
	double max;
	int max_i, max_j;

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	ma = modulus(ax, ay, az);
	mb = modulus(bx, by, bz);
	mc = modulus(cx, cy, cz);

	pha = get_model_phase(ax/ma, ay/ma, az/ma, image->features, nel, pmax,
	                      fft_in, fft_out, plan, smin, smax, rg, det);
	phb = get_model_phase(bx/mb, by/mb, bz/mb, image->features, nel, pmax,
	                      fft_in, fft_out, plan, smin, smax, rg, det);
	phc = get_model_phase(cx/mc, cy/mc, cz/mc, image->features, nel, pmax,
	                      fft_in, fft_out, plan, smin, smax, rg, det);

	for ( i=0; i<det->n_panels; i++ ) {
		if ( det->panels[i].rigid_group == rg ) {
			p = &det->panels[i];
			break;
		}
	}

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


static void refine_all_rigid_groups(struct image *image, UnitCell *cell,
                                    int nel, double pmax,
                                    double *fft_in, fftw_complex *fft_out,
                                    fftw_plan plan, int smin, int smax,
                                    struct detector *det,
                                    struct reax_private *p)
{
	int i;

	for ( i=0; i<image->det->num_rigid_groups; i++ ) {
		refine_rigid_group(image, cell, image->det->rigid_groups[i],
		                   nel, pmax, fft_in, fft_out, plan, smin, smax,
		                   det, p);
	}
}


void reax_index(IndexingPrivate *pp, struct image *image, UnitCell *cell)
{
	int i;
	struct reax_private *p;
	double fom;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double mod_a, mod_b, mod_c;
	double al, be, ga;
	double th, ph;
	double *fft_in;
	fftw_complex *fft_out;
	int smin, smax;
	double amin, amax;
	double bmin, bmax;
	double cmin, cmax;
	double pmax;
	int n;
	const double ltol = 5.0;             /* Direct space axis length
	                                      * tolerance in percent */
	const double angtol = deg2rad(1.5);  /* Direct space angle tolerance
	                                      * in radians */

	assert(pp->indm == INDEXING_REAX);
	p = (struct reax_private *)pp;

	fft_in = fftw_malloc(p->nel*sizeof(double));
	fft_out = fftw_malloc((p->nel/2 + 1)*sizeof(fftw_complex));

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

	al = angle_between(bx, by, bz, cx, cy, cz);
	be = angle_between(ax, ay, az, cx, cy, cz);
	ga = angle_between(ax, ay, az, bx, by, bz);

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

	/* Sanity check */
	if ( pmax < 1e4 ) return;

	/* Search for a */
	smin = 2.0*pmax * amin;
	smax = 2.0*pmax * amax;
	fom = 0.0;  th = 0.0;  ph = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom;

		new_fom = check_dir(&p->directions[i], image->features,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    smin, smax, NULL, NULL);
		if ( new_fom > fom ) {
			fom = new_fom;
			th = p->directions[i].th;
			ph = p->directions[i].ph;
		}

	}
	fine_search(p, image->features, p->nel, pmax, fft_in, fft_out,
	            p->plan, smin, smax, th, ph, &ax, &ay, &az);

	/* Search for b */
	smin = 2.0*pmax * bmin;
	smax = 2.0*pmax * bmax;
	fom = 0.0;  th = 0.0;  ph = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom, ang;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, ax, ay, az);
		if ( fabs(ang-ga) > angtol ) continue;

		new_fom = check_dir(&p->directions[i], image->features,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    smin, smax, NULL, NULL);
		if ( new_fom > fom ) {
			fom = new_fom;
			th = p->directions[i].th;
			ph = p->directions[i].ph;
		}

	}
	fine_search(p, image->features, p->nel, pmax, fft_in, fft_out,
	            p->plan, smin, smax, th, ph, &bx, &by, &bz);

	/* Search for c */
	smin = 2.0*pmax * cmin;
	smax = 2.0*pmax * cmax;
	fom = 0.0;  th = 0.0;  ph = 0.0;
	for ( i=0; i<p->n_dir; i++ ) {

		double new_fom, ang;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, ax, ay, az);
		if ( fabs(ang-be) > angtol ) continue;

		ang = angle_between(p->directions[i].x, p->directions[i].y,
		                    p->directions[i].z, bx, by, bz);
		if ( fabs(ang-al) > angtol ) continue;

		new_fom = check_dir(&p->directions[i], image->features,
		                    p->nel, pmax, fft_in, fft_out, p->plan,
		                    smin, smax, NULL, NULL);
		if ( new_fom > fom ) {
			fom = new_fom;
			th = p->directions[i].th;
			ph = p->directions[i].ph;
		}

	}
	fine_search(p, image->features, p->nel, pmax, fft_in, fft_out,
	            p->plan, smin, smax, th, ph, &cx, &cy, &cz);

	image->candidate_cells[0] = cell_new();
	cell_set_cartesian(image->candidate_cells[0],
	                   ax, ay, az, bx, by, bz, cx, cy, cz);

	refine_all_rigid_groups(image, image->candidate_cells[0], p->nel, pmax,
	                        fft_in, fft_out, p->plan, smin, smax,
	                        image->det, p);
	map_all_peaks(image);
	refine_cell(image, image->candidate_cells[0], image->features);

	fftw_free(fft_in);
	fftw_free(fft_out);

	image->ncells = 1;
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

	p->cw = 128; p->ch = 128;

	/* Also not used */
	p->r_fft_in = fftw_malloc(p->cw*p->ch*sizeof(fftw_complex));
	p->r_fft_out = fftw_malloc(p->cw*p->ch*sizeof(fftw_complex));

	p->r_plan = fftw_plan_dft_2d(p->cw, p->ch, p->r_fft_in, p->r_fft_out,
	                             1, FFTW_MEASURE);

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

	fftw_destroy_plan(p->r_plan);
	fftw_free(p->r_fft_in);
	fftw_free(p->r_fft_out);

	free(p);
}
