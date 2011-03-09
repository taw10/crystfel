/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
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
#include <string.h>
#include <assert.h>
#include <gsl/gsl_statistics_int.h>
#include <pthread.h>
#include <fenv.h>

#include "image.h"
#include "utils.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "diffraction.h"


/* How close a peak must be to an indexed position to be considered "close"
 * for the purposes of double hit detection and sanity checking. */
#define PEAK_CLOSE (30.0)

/* How close a peak must be to an indexed position to be considered "close"
 * for the purposes of integration. */
#define PEAK_REALLY_CLOSE (10.0)

/* Degree of polarisation of X-ray beam */
#define POL (1.0)

/* Window size for Zaefferer peak detection */
#define PEAK_WINDOW_SIZE (10)


static int cull_peaks_in_panel(struct image *image, struct panel *p)
{
	int i, n;
	int nelim = 0;

	n = image_feature_count(image->features);

	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		int j, ncol;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		if ( f->x < p->min_fs ) continue;
		if ( f->x > p->max_fs ) continue;
		if ( f->y < p->min_ss ) continue;
		if ( f->y > p->max_ss ) continue;

		/* How many peaks are in the same column? */
		ncol = 0;
		for ( j=0; j<n; j++ ) {

			struct imagefeature *g;

			if ( i==j ) continue;

			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;

			if ( p->badrow == 'x' ) {
				if ( fabs(f->y - g->y) < 2.0 ) ncol++;
			} else if ( p->badrow == 'y' ) {
				if ( fabs(f->x - g->x) < 2.0 ) ncol++;
			} /* else do nothing */

		}

		/* More than three? */
		if ( ncol <= 3 ) continue;

		/* Yes?  Delete them all... */
		nelim = 0;
		for ( j=0; j<n; j++ ) {
			struct imagefeature *g;
			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;
			if ( p->badrow == 'x' ) {
				if ( fabs(f->y - g->y) < 2.0 ) {
					image_remove_feature(image->features,
					                     j);
					nelim++;
				}
			} else if ( p->badrow == 'y' ) {
				if ( fabs(f->x - g->x) < 2.0 ) {
					image_remove_feature(image->features,
					                     j);
					nelim++;
				}
			} else {
				ERROR("Invalid badrow direction.\n");
				abort();
			}

		}

	}

	return nelim;
}


/* Post-processing of the peak list to remove noise */
static int cull_peaks(struct image *image)
{
	int nelim = 0;
	struct panel *p;
	int i;

	for ( i=0; i<image->det->n_panels; i++ ) {
		p = &image->det->panels[i];
		if ( p->badrow != '-' ) {
			nelim += cull_peaks_in_panel(image, p);
		}
	}

	return nelim;
}


/* Returns non-zero if peak has been vetoed.
 * i.e. don't use result if return value is not zero. */
int integrate_peak(struct image *image, int xp, int yp,
                   double *xc, double *yc, double *intensity,
                   double *pbg, double *pmax,
                   int do_polar, int centroid)
{
	signed int x, y;
	int lim, out_lim;
	double total = 0.0;
	double xct = 0.0;
	double yct = 0.0;
	double noise = 0.0;
	int noise_counts = 0;
	double max = 0.0;
	struct panel *p = NULL;

	p = find_panel(image->det, xp, yp);
	if ( p == NULL ) return 1;
	if ( p->no_index ) return 1;

	lim = p->peak_sep/4.0;
	out_lim = lim + 1;

	for ( x=-out_lim; x<+out_lim; x++ ) {
	for ( y=-out_lim; y<+out_lim; y++ ) {

		double val;
		double tt = 0.0;
		double phi, pa, pb, pol;
		uint16_t flags;
		struct panel *p2;

		/* Outer mask radius */
		if ( x*x + y*y > out_lim ) continue;

		if ( ((x+xp)>=image->width) || ((x+xp)<0) ) continue;
		if ( ((y+yp)>=image->height) || ((y+yp)<0) ) continue;

		/* Strayed off one panel? */
		p2 = find_panel(image->det, x+xp, y+yp);
		if ( p2 != p ) return 1;

		/* Veto this peak if we tried to integrate in a bad region */
		if ( image->flags != NULL ) {
			flags = image->flags[(x+xp)+image->width*(y+yp)];
			if ( !(flags & 0x01) ) return 1;
		}

		val = image->data[(x+xp)+image->width*(y+yp)];

		/* Inner mask */
		if ( x*x + y*y > lim ) {
			/* Estimate noise from this region */
			noise += fabs(val);
			noise_counts++;
			continue;
		}

		if ( val > max ) max = val;

		if ( do_polar ) {

			tt = get_tt(image, x+xp, y+yp);

			phi = atan2(y+yp, x+xp);
			pa = pow(sin(phi)*sin(tt), 2.0);
			pb = pow(cos(tt), 2.0);
			pol = 1.0 - 2.0*POL*(1-pa) + POL*(1.0+pb);

			val /= pol;

		}

		total += val;

		xct += val*(xp+x);
		yct += val*(yp+y);

	}
	}

	/* The centroid is excitingly undefined if there is no intensity */
	if ( centroid && (total != 0) ) {
		*xc = (double)xct / total;
		*yc = (double)yct / total;
		*intensity = total;
	} else {
		*xc = (double)xp;
		*yc = (double)yp;
		*intensity = total;
	}

	if ( in_bad_region(image->det, *xc, *yc) ) {
		return 1;
	}

	if ( pbg != NULL ) {
		*pbg = (noise / noise_counts);
	}
	if ( pmax != NULL ) {
		*pmax = max;
	}

	return 0;
}


static void search_peaks_in_panel(struct image *image, float threshold,
                                  float min_gradient, struct panel *p)
{
	int fs, ss, stride;
	float *data;
	double d;
	int idx;
	double f_fs = 0.0;
	double f_ss = 0.0;
	double intensity = 0.0;
	int nrej_dis = 0;
	int nrej_hot = 0;
	int nrej_pro = 0;
	int nrej_fra = 0;
	int nrej_bad = 0;
	int nacc = 0;
	int ncull;

	data = image->data;
	stride = image->width;

	for ( fs = p->min_fs+1; fs <= p->max_fs-1; fs++ ) {
	for ( ss = p->min_ss+1; ss <= p->max_ss-1; ss++ ) {

		double dx1, dx2, dy1, dy2;
		double dxs, dys;
		double grad;
		int mask_fs, mask_ss;
		int s_fs, s_ss;
		double max;
		unsigned int did_something;
		int r;

		/* Overall threshold */
		if ( data[fs+stride*ss] < threshold ) continue;

		/* Get gradients */
		dx1 = data[fs+stride*ss] - data[(fs+1)+stride*ss];
		dx2 = data[(fs-1)+stride*ss] - data[fs+stride*ss];
		dy1 = data[fs+stride*ss] - data[(fs+1)+stride*(ss+1)];
		dy2 = data[fs+stride*(ss-1)] - data[fs+stride*ss];

		/* Average gradient measurements from both sides */
		dxs = ((dx1*dx1) + (dx2*dx2)) / 2;
		dys = ((dy1*dy1) + (dy2*dy2)) / 2;

		/* Calculate overall gradient */
		grad = dxs + dys;

		if ( grad < min_gradient ) continue;

		mask_fs = fs;
		mask_ss = ss;

		do {

			max = data[mask_fs+stride*mask_ss];
			did_something = 0;

			for ( s_ss=biggest(mask_ss-PEAK_WINDOW_SIZE/2,
			                   p->min_ss);
			      s_ss<=smallest(mask_ss+PEAK_WINDOW_SIZE/2,
			                     p->max_ss);
			      s_ss++ ) {
			for ( s_fs=biggest(mask_fs-PEAK_WINDOW_SIZE/2,
			                   p->min_fs);
			      s_fs<=smallest(mask_fs+PEAK_WINDOW_SIZE/2,
			                     p->max_fs);
			      s_fs++ ) {

				if ( data[s_fs+stride*s_ss] > max ) {
					max = data[s_fs+stride*s_ss];
					mask_fs = s_fs;
					mask_ss = s_ss;
					did_something = 1;
				}

			}
			}

			/* Abort if drifted too far from the foot point */
			if ( distance(mask_fs, mask_ss, fs, ss)
			   > p->peak_sep ) {
				break;
			}

		} while ( did_something );

		/* Too far from foot point? */
		if ( distance(mask_fs, mask_ss, fs, ss) > p->peak_sep ) {
			nrej_dis++;
			continue;
		}

		/* Should be enforced by bounds used above.  Muppet check. */
		assert(mask_fs <= p->max_fs);
		assert(mask_ss <= p->max_ss);
		assert(mask_fs >= p->min_fs);
		assert(mask_ss >= p->min_ss);

		/* Centroid peak and get better coordinates.
		 * Don't bother doing polarisation/SA correction, because the
		 * intensity of this peak is only an estimate at this stage. */
		r = integrate_peak(image, mask_fs, mask_ss,
		                   &f_fs, &f_ss, &intensity, NULL, NULL, 0, 1);
		if ( r ) {
			/* Bad region - don't detect peak */
			nrej_bad++;
			continue;
		}

		/* It is possible for the centroid to fall outside the image */
		if ( (f_fs < p->min_fs) || (f_fs > p->max_fs)
		  || (f_ss < p->min_ss) || (f_ss > p->max_ss) ) {
			nrej_fra++;
			continue;
		}

		/* Check for a nearby feature */
		image_feature_closest(image->features, f_fs, f_ss, &d, &idx);
		if ( d < p->peak_sep ) {
			nrej_pro++;
			continue;
		}

		/* Add using "better" coordinates */
		image_add_feature(image->features, f_fs, f_ss, image, intensity,
		                  NULL);
		nacc++;

	}
	}

	if ( image->det != NULL ) {
		ncull = cull_peaks(image);
		nacc -= ncull;
	} else {
		STATUS("Not culling peaks because I don't have a "
		       "detector geometry file.\n");
		ncull = 0;
	}

	STATUS("%i accepted, %i box, %i hot, %i proximity, %i outside panel, "
	       "%i in bad regions, %i badrow culled.\n",
	       nacc, nrej_dis, nrej_hot, nrej_pro, nrej_fra, nrej_bad, ncull);
}


void search_peaks(struct image *image, float threshold, float min_gradient)
{
	int i;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();


	for ( i=0; i<image->det->n_panels; i++ ) {

		struct panel *p = &image->det->panels[i];

		if ( p->no_index ) continue;
		search_peaks_in_panel(image, threshold, min_gradient, p);

	}
}


void dump_peaks(struct image *image, FILE *ofh, pthread_mutex_t *mutex)
{
	int i;

	/* Get exclusive access to the output stream if necessary */
	if ( mutex != NULL ) pthread_mutex_lock(mutex);

	fprintf(ofh, "Peaks from peak search in %s\n", image->filename);
	fprintf(ofh, "  x/px     y/px   (1/d)/nm^-1    Intensity\n");

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec r;
		double q;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;

		r = get_q(image, f->x, f->y, NULL, 1.0/image->lambda);
		q = modulus(r.u, r.v, r.w);

		fprintf(ofh, "%8.3f %8.3f %8.3f    %12.3f\n",
		       f->x, f->y, q/1.0e9, f->intensity);

	}

	fprintf(ofh, "\n");

	if ( mutex != NULL ) pthread_mutex_unlock(mutex);
}


RefList *find_projected_peaks(struct image *image, UnitCell *cell,
                              int circular_domain, double domain_r)
{
	int fs, ss;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	RefList *reflections;
	double alen, blen, clen;
	int n_reflections = 0;

	reflections = reflist_new();

	/* "Borrow" direction values to get reciprocal lengths */
	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	alen = modulus(ax, ay, az);
	blen = modulus(bx, by, bz);
	clen = modulus(cx, cy, cz);

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	fesetround(1);  /* Round towards nearest */
	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		double hd, kd, ld;  /* Indices with decimal places */
		double dh, dk, dl;  /* Distances in h,k,l directions */
		signed int h, k, l;
		struct rvec q;
		double dist;
		Reflection *refl;
		double cur_dist;

		q = get_q(image, fs, ss, NULL, 1.0/image->lambda);

		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;

		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		dh = hd - h;
		dk = kd - k;
		dl = ld - l;

		if ( circular_domain ) {
			/* Circular integration domain */
			dist = sqrt(pow(dh*alen, 2.0) + pow(dk*blen, 2.0)
			                              + pow(dl*clen, 2.0));
			if ( dist > domain_r ) continue;
		} else {
			/* "Crystallographic" integration domain */
			dist = sqrt(pow(dh, 2.0) + pow(dk, 2.0) + pow(dl, 2.0));
			if ( dist > domain_r ) continue;
		}

		refl = find_refl(reflections, h, k, l);
		if ( refl != NULL ) {
			cur_dist = get_excitation_error(refl);
			if ( dist < cur_dist ) {
				set_detector_pos(refl, dist, fs, ss);
			}
		} else {
			Reflection *new;
			new = add_refl(reflections, h, k, l);
			set_detector_pos(new, dist, fs, ss);
			n_reflections++;
		}

	}
	}

	optimise_reflist(reflections);

	STATUS("Found %i reflections\n", n_reflections);
	return reflections;
}


int peak_sanity_check(struct image *image, UnitCell *cell,
                      int circular_domain, double domain_r)
{
	int i;
	int n_feat = 0;
	int n_sane = 0;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double aslen, bslen, cslen;

	/* "Borrow" direction values to get reciprocal lengths */
	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	aslen = modulus(ax, ay, az);
	bslen = modulus(bx, by, bz);
	cslen = modulus(cx, cy, cz);

	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	fesetround(1);  /* Round towards nearest */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		double dist;
		struct rvec q;
		struct imagefeature *f;
		double hd, kd, ld;
		signed int h, k, l;
		double dh, dk, dl;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		n_feat++;

		/* Get closest hkl */
		q = get_q(image, f->x, f->y, NULL, 1.0/image->lambda);

		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;

		h = lrint(hd);  k = lrint(kd);  l = lrint(ld);

		dh = hd - h;  dk = kd - k;  dl = ld - l;

		if ( circular_domain ) {

			/* Circular integration domain */
			dist = sqrt(pow(dh*aslen, 2.0) + pow(dk*bslen, 2.0)
			                              + pow(dl*cslen, 2.0));
			if ( dist <= domain_r ) n_sane++;

		} else {

			/* "Crystallographic" integration domain */
			dist = sqrt(pow(dh, 2.0) + pow(dk, 2.0) + pow(dl, 2.0));
			if ( dist <= domain_r ) n_sane++;
		}



	}

	STATUS("Sanity factor: %f / %f = %f\n", (float)n_sane, (float)n_feat,
                                                (float)n_sane / (float)n_feat);
	if ( (float)n_sane / (float)n_feat < 0.1 ) return 0;

	return 1;
}


static void output_header(FILE *ofh, UnitCell *cell, struct image *image)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	double a, b, c, al, be, ga;

	fprintf(ofh, "Reflections from indexing in %s\n", image->filename);

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);
	fprintf(ofh, "Cell parameters %7.5f %7.5f %7.5f nm, %7.5f %7.5f %7.5f deg\n",
	       a*1.0e9, b*1.0e9, c*1.0e9,
	       rad2deg(al), rad2deg(be), rad2deg(ga));

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	fprintf(ofh, "astar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       asx/1e9, asy/1e9, asz/1e9);
	fprintf(ofh, "bstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       bsx/1e9, bsy/1e9, bsz/1e9);
	fprintf(ofh, "cstar = %+9.7f %+9.7f %+9.7f nm^-1\n",
	       csx/1e9, csy/1e9, csz/1e9);

	if ( image->f0_available ) {
		fprintf(ofh, "f0 = %7.5f (arbitrary gas detector units)\n",
		       image->f0);
	} else {
		fprintf(ofh, "f0 = invalid\n");
	}

	fprintf(ofh, "photon_energy_eV = %f\n",
	        J_to_eV(ph_lambda_to_en(image->lambda)));

}


void output_intensities(struct image *image, UnitCell *cell,
                        RefList *reflections, pthread_mutex_t *mutex, int polar,
                        int use_closer, FILE *ofh)
{
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;
	Reflection *refl;
	RefListIterator *iter;

	/* Get exclusive access to the output stream if necessary */
	if ( mutex != NULL ) pthread_mutex_lock(mutex);

	output_header(ofh, cell, image);

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double x, y, intensity;
		double d;
		int idx;
		double bg, max;
		struct panel *p;
		double px, py;
		signed int h, k, l;

		get_detector_pos(refl, &px, &py);
		p = find_panel(image->det, px, py);
		if ( p == NULL ) continue;
		if ( p->no_index ) continue;

		/* Wait.. is there a really close feature which was detected? */
		if ( use_closer ) {

			struct imagefeature *f;

			if ( image->features != NULL ) {
				f = image_feature_closest(image->features,
					                  px, py, &d, &idx);
			} else {
				f = NULL;
			}
			if ( (f != NULL) && (d < PEAK_REALLY_CLOSE) ) {

				int r;

				/* f->intensity was measured on the filtered
				 * pattern, so instead re-integrate using old
				 * coordinates. This will produce further
				 * revised coordinates. */
				r = integrate_peak(image, f->x, f->y, &x, &y,
					           &intensity, &bg, &max,
					           polar, 1);
				if ( r ) {
					/* The original peak (which also went
					 * through integrate_peak(), but with
					 * the mangled image data) would have
					 * been rejected if it was in a bad
					 * region.  Integration of the same
					 * peak included a bad region this time.
					 */
					continue;
				}
				intensity = f->intensity;

			} else {

				int r;

				r = integrate_peak(image, px, py, &x, &y,
				                   &intensity, &bg, &max,
					           polar, 1);
				if ( r ) {
					/* Plain old ordinary peak veto */
					continue;
				}

			}

		} else {

			int r;

			r = integrate_peak(image, px, py, &x, &y,
			                   &intensity, &bg, &max, polar, 0);
			if ( r ) {
				/* Plain old ordinary peak veto */
				continue;
			}

		}

		/* Write h,k,l, integrated intensity and centroid coordinates */
		get_indices(refl, &h, &k, &l);
		fprintf(ofh, "%3i %3i %3i %6f (at %5.2f,%5.2f) max=%6f bg=%6f\n",
		        h, k, l, intensity, x, y, max, bg);

	}

	/* Blank line at end */
	fprintf(ofh, "\n");

	if ( mutex != NULL ) pthread_mutex_unlock(mutex);
}


void output_pixels(struct image *image, UnitCell *cell,
                   pthread_mutex_t *mutex, int do_polar,
                   FILE *ofh, int circular_domain, double domain_r)
{
	int i;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	int fs, ss;
	double aslen, bslen, cslen;
	double *intensities;
	double *xmom;
	double *ymom;
	ReflItemList *obs;

	/* Get exclusive access to the output stream if necessary */
	if ( mutex != NULL ) pthread_mutex_lock(mutex);

	output_header(ofh, cell, image);

	obs = new_items();
	intensities = new_list_intensity();
	xmom = new_list_intensity();
	ymom = new_list_intensity();

	/* "Borrow" direction values to get reciprocal lengths */
	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);
	aslen = modulus(ax, ay, az);
	bslen = modulus(bx, by, bz);
	cslen = modulus(cx, cy, cz);

	cell_get_cartesian(cell, &ax, &ay, &az,
	                         &bx, &by, &bz,
	                         &cx, &cy, &cz);
	/* For each pixel */
	fesetround(1);  /* Round towards nearest */
	for ( fs=0; fs<image->width; fs++ ) {
	for ( ss=0; ss<image->height; ss++ ) {

		double hd, kd, ld;  /* Indices with decimal places */
		double dh, dk, dl;  /* Distances in h,k,l directions */
		signed int h, k, l;
		struct rvec q;
		double dist;
		struct panel *p;
		double twotheta;

		p = find_panel(image->det, fs, ss);
		if ( p == NULL ) continue;
		if ( p->no_index ) continue;

		q = get_q(image, fs, ss, &twotheta, 1.0/image->lambda);

		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;

		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		dh = hd - h;  dk = kd - k;  dl = ld - l;

		if ( circular_domain ) {

			/* Circular integration domain */
			dist = sqrt(pow(dh*aslen, 2.0) + pow(dk*bslen, 2.0)
			                              + pow(dl*cslen, 2.0));

		} else {

			/* "Crystallographic" integration domain */
			dist = sqrt(pow(dh, 2.0) + pow(dk, 2.0) + pow(dl, 2.0));
		}

		if ( dist < domain_r ) {

			double val;
			double pix_area, Lsq, proj_area, dsq, sa;
			double phi, pa, pb, pol;
			double xs, ys, rx, ry;

			/* Veto if we want to integrate a bad region */
			if ( image->flags != NULL ) {
				int flags;
				flags = image->flags[fs+image->width*ss];
				if ( !(flags & 0x01) ) continue;
			}

			val = image->data[fs+image->width*ss];

			/* Area of one pixel */
			pix_area = pow(1.0/p->res, 2.0);
			Lsq = pow(p->clen, 2.0);

			/* Area of pixel as seen from crystal  */
			proj_area = pix_area * cos(twotheta);

			/* Calculate distance from crystal to pixel */
			xs = (fs-p->min_fs)*p->fsx + (ss-p->min_ss)*p->ssx;
			ys = (fs-p->min_fs)*p->fsy + (ss-p->min_ss)*p->ssy;
			rx = (xs + p->cnx) / p->res;
			ry = (ys + p->cny) / p->res;
			dsq = sqrt(pow(rx, 2.0) + pow(ry, 2.0));

			/* Projected area of pixel / distance squared */
			sa = 1.0e7 * proj_area / (dsq + Lsq);

			/* Solid angle correction is needed in this case */
			val /= sa;

			if ( do_polar ) {

				phi = atan2(ry, rx);
				pa = pow(sin(phi)*sin(twotheta), 2.0);
				pb = pow(cos(twotheta), 2.0);
				pol = 1.0 - 2.0*POL*(1-pa) + POL*(1.0+pb);

				val /= pol;

			}

			/* Add value to sum */
			integrate_intensity(intensities, h, k, l, val);

			integrate_intensity(xmom, h, k, l, val*fs);
			integrate_intensity(ymom, h, k, l, val*ss);

			if ( !find_item(obs, h, k, l) ) {
				add_item(obs, h, k, l);
			}

		}

	}
	}

	for ( i=0; i<num_items(obs); i++ ) {

		struct refl_item *it;
		double intensity, xmomv, ymomv;
		double xp, yp;

		it = get_item(obs, i);
		intensity = lookup_intensity(intensities, it->h, it->k, it->l);
		xmomv = lookup_intensity(xmom, it->h, it->k, it->l);
		ymomv = lookup_intensity(ymom, it->h, it->k, it->l);

		xp = xmomv / (double)intensity;
		yp = ymomv / (double)intensity;

		fprintf(ofh, "%3i %3i %3i %6f (at %5.2f,%5.2f)\n",
		       it->h, it->k, it->l, intensity, xp, yp);

	}

	fprintf(ofh, "No peak statistics, because output_pixels() was used.\n");
	/* Blank line at end */
	fprintf(ofh, "\n");

	free(xmom);
	free(ymom);
	free(intensities);
	delete_items(obs);

	if ( mutex != NULL ) pthread_mutex_unlock(mutex);
}
