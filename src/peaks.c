/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *	    2011 Andrew Martin
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

		if ( f->fs < p->min_fs ) continue;
		if ( f->fs > p->max_fs ) continue;
		if ( f->ss < p->min_ss ) continue;
		if ( f->ss > p->max_ss ) continue;

		/* How many peaks are in the same column? */
		ncol = 0;
		for ( j=0; j<n; j++ ) {

			struct imagefeature *g;

			if ( i==j ) continue;

			g = image_get_feature(image->features, j);
			if ( g == NULL ) continue;

			if ( p->badrow == 'f' ) {
				if ( fabs(f->ss - g->ss) < 2.0 ) ncol++;
			} else if ( p->badrow == 's' ) {
				if ( fabs(f->fs - g->fs) < 2.0 ) ncol++;
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
			if ( p->badrow == 'f' ) {
				if ( fabs(f->ss - g->ss) < 2.0 ) {
					image_remove_feature(image->features,
					                     j);
					nelim++;
				}
			} else if ( p->badrow == 's' ) {
				if ( fabs(f->fs - g->ss) < 2.0 ) {
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
int integrate_peak(struct image *image, int cfs, int css,
                   double *pfs, double *pss, double *intensity,
                   double *pbg, double *pmax, double *sigma,
                   int do_polar, int centroid, int bgsub)
{
	signed int fs, ss;
	double lim, out_lim;
	double lim_sq, out_lim_sq;
	double total = 0.0;
	double fsct = 0.0;
	double ssct = 0.0;
	double noise = 0.0;
	int noise_counts = 0;
	double max = 0.0;
	struct panel *p = NULL;
        int pixel_counts = 0;
        double noise_mean = 0.0;
	double noise_meansq = 0.0;

	p = find_panel(image->det, cfs, css);
	if ( p == NULL ) return 1;
	if ( p->no_index ) return 1;

	lim = p->integr_radius;
	out_lim = 2.0 + lim;
	lim_sq = pow(lim, 2.0);
	out_lim_sq = pow(out_lim, 2.0);

	for ( fs=-out_lim; fs<+out_lim; fs++ ) {
	for ( ss=-out_lim; ss<+out_lim; ss++ ) {

		double val;
		double tt = 0.0;
		double phi, pa, pb, pol;
		uint16_t flags;
		struct panel *p2;
		int idx;

		/* Outer mask radius */
		if ( fs*fs + ss*ss > out_lim_sq ) continue;

		if ( ((fs+cfs)>=image->width) || ((fs+cfs)<0) ) continue;
		if ( ((ss+css)>=image->height) || ((ss+css)<0) ) continue;

		/* Strayed off one panel? */
		p2 = find_panel(image->det, fs+cfs, ss+css);
		if ( p2 != p ) return 1;

		idx = fs+cfs+image->width*(ss+css);

		/* Veto this peak if we tried to integrate in a bad region */
		if ( image->flags != NULL ) {

			flags = image->flags[idx];

			/* It must have all the "good" bits to be valid */
			if ( !((flags & image->det->mask_good)
			                   == image->det->mask_good) ) return 1;

			/* If it has any of the "bad" bits, reject */
			if ( flags & image->det->mask_bad ) return 1;

		}

		val = image->data[idx];

		if ( do_polar ) {

			tt = get_tt(image, fs+cfs, ss+css);

			phi = atan2(ss+css, fs+cfs);
			pa = pow(sin(phi)*sin(tt), 2.0);
			pb = pow(cos(tt), 2.0);
			pol = 1.0 - 2.0*POL*(1-pa) + POL*(1.0+pb);

			val /= pol;

		}

		if ( val > max ) max = val;

		/* If outside inner mask, estimate noise from this region */
		if ( fs*fs + ss*ss > lim_sq ) {

			/* Noise */
			noise += val;
			noise_counts++;
			noise_meansq += pow(val, 2.0);

		} else {

			/* Peak */
			pixel_counts++;
			total += val;
			fsct += val*(cfs+fs);
			ssct += val*(css+ss);
		}

	}
	}

        noise_mean = noise / noise_counts;

	/* The centroid is excitingly undefined if there is no intensity */
	if ( centroid && (total != 0) ) {
		*pfs = ((double)fsct / total) + 0.5;
		*pss = ((double)ssct / total) + 0.5;
	} else {
		*pfs = (double)cfs + 0.5;
		*pss = (double)css + 0.5;
	}
	if ( bgsub ) {
		*intensity = total - pixel_counts*noise_mean;
	} else {
		*intensity = total;
	}

	if ( in_bad_region(image->det, *pfs, *pss) ) return 1;

	if ( sigma != NULL ) {
		/* First term is standard deviation of background per pixel
		 * sqrt(pixel_counts) - increase of error for integrated value
		 * sqrt(2) - increase of error for background subtraction  */
		*sigma = sqrt(noise_meansq/noise_counts-(noise_mean*noise_mean))
		          * sqrt(2.0*pixel_counts);
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
			if ( distance(mask_fs, mask_ss, fs, ss) >
			     p->peak_sep/2.0 )
			{
				break;
			}

		} while ( did_something );

		/* Too far from foot point? */
		if ( distance(mask_fs, mask_ss, fs, ss) > p->peak_sep/2.0 ) {
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
		                   &f_fs, &f_ss, &intensity,
		                   NULL, NULL, NULL, 0, 1, 0);
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
		if ( d < p->peak_sep/2.0 ) {
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

//	STATUS("%i accepted, %i box, %i proximity, %i outside panel, "
//	       "%i in bad regions, %i badrow culled.\n",
//	       nacc, nrej_dis, nrej_pro, nrej_fra, nrej_bad, ncull);
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
		q = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);

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

	if ( (float)n_sane / (float)n_feat < 0.1 ) return 0;

	return 1;
}


/* Integrate the list of predicted reflections in "image" */
void integrate_reflections(struct image *image, int polar, int use_closer,
                           int bgsub)
{
	Reflection *refl;
	RefListIterator *iter;

	for ( refl = first_refl(image->reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		double fs, ss, intensity;
		double d;
		int idx;
		double bg, max;
		double sigma;
		struct panel *p;
		double pfs, pss;
		int r;

		get_detector_pos(refl, &pfs, &pss);
		p = find_panel(image->det, pfs, pss);
		if ( p == NULL ) continue;
		if ( p->no_index ) continue;

		/* Is there a really close feature which was detected? */
		if ( use_closer ) {

			struct imagefeature *f;

			if ( image->features != NULL ) {
				f = image_feature_closest(image->features,
					                  pfs, pss, &d, &idx);
			} else {
				f = NULL;
			}
			if ( (f != NULL) && (d < PEAK_REALLY_CLOSE) ) {

				pfs = f->fs;
				pss = f->ss;

			}
		}

		r = integrate_peak(image, pfs, pss, &fs, &ss,
		                   &intensity, &bg, &max, &sigma, polar, 0,
		                   bgsub);

		/* Record intensity and set redundancy to 1 on success */
		if ( r == 0 ) {
			set_int(refl, intensity);
			set_esd_intensity(refl, sigma);
			set_redundancy(refl, 1);
		}

	}
}
