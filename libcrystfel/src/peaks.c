/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
 * Copyright © 2012 Andrew Martin <andrew.martin@desy.de>
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
#include <string.h>
#include <assert.h>
#include <gsl/gsl_statistics_int.h>
#include <pthread.h>
#include <fenv.h>

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "reflist-utils.h"
#include "beam-parameters.h"


/* How close a peak must be to an indexed position to be considered "close"
 * for the purposes of double hit detection and sanity checking. */
#define PEAK_CLOSE (30.0)

/* How close a peak must be to an indexed position to be considered "close"
 * for the purposes of integration. */
#define PEAK_REALLY_CLOSE (10.0)

/* Degree of polarisation of X-ray beam */
#define POL (1.0)

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
                   int centroid, int bgsub)
{
	signed int fs, ss;
	double lim, out_lim, mid_lim;
	double lim_sq, out_lim_sq, mid_lim_sq;
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
	struct beam_params *beam;
	double aduph;

	beam = image->beam;
	if ( beam != NULL ) {
		aduph = image->beam->adu_per_photon;
	} else {
		aduph = 1.0;
	}

	p = find_panel(image->det, cfs, css);
	if ( p == NULL ) return 1;
	if ( p->no_index ) return 1;

	lim = p->integr_radius;
	mid_lim = 3.0 + lim;
	out_lim = 6.0 + lim;
	lim_sq = pow(lim, 2.0);
	mid_lim_sq = pow(mid_lim, 2.0);
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

		if ( val > max ) max = val;

		/* If outside inner mask, estimate noise from this region */
		if ( fs*fs + ss*ss > mid_lim_sq ) {

			/* Noise
			 * noise and noise_meansq are both in photons (^2) */
			noise += val / aduph;
			noise_counts++;
			noise_meansq += pow(val/aduph, 2.0);

		} else if ( fs*fs + ss*ss < lim_sq ) {

			/* Peak */
			pixel_counts++;
			total += val;
			fsct += val*(cfs+fs);
			ssct += val*(css+ss);

		}

	}
	}

	noise_mean = noise / noise_counts;  /* photons */

	/* The centroid is excitingly undefined if there is no intensity */
	if ( centroid && (total != 0) ) {
		*pfs = ((double)fsct / total) + 0.5;
		*pss = ((double)ssct / total) + 0.5;
	} else {
		*pfs = (double)cfs + 0.5;
		*pss = (double)css + 0.5;
	}
	if ( bgsub ) {
		*intensity = total - aduph * pixel_counts*noise_mean; /* ADU */
	} else {
		*intensity = total;  /* ADU */
	}

	if ( in_bad_region(image->det, *pfs, *pss) ) return 1;

	if ( sigma != NULL ) {

		/* First term is standard deviation of background per pixel
		 * sqrt(pixel_counts) - increase of error for integrated value
		 * sqrt(2) - increase of error for background subtraction  */
		*sigma = sqrt(noise_meansq/noise_counts-(noise_mean*noise_mean))
		          * sqrt(2.0*pixel_counts) * aduph;

	}

	if ( pbg != NULL ) {
		*pbg = aduph * (noise / noise_counts);
	}
	if ( pmax != NULL ) {
		*pmax = max;
	}

	return 0;
}


static void search_peaks_in_panel(struct image *image, float threshold,
                                  float min_gradient, float min_snr,
                                  struct panel *p)
{
	int fs, ss, stride;
	float *data;
	double d;
	int idx;
	double f_fs = 0.0;
	double f_ss = 0.0;
	double intensity = 0.0;
	double sigma = 0.0;
	double pbg = 0.0;
	double pmax = 0.0;
	int nrej_dis = 0;
	int nrej_pro = 0;
	int nrej_fra = 0;
	int nrej_bad = 0;
	int nrej_snr = 0;
	int nacc = 0;
	int ncull;
	const int pws = p->peak_sep/2;

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

			for ( s_ss=biggest(mask_ss-pws/2,
			                   p->min_ss);
			      s_ss<=smallest(mask_ss+pws/2,
			                     p->max_ss);
			      s_ss++ ) {
			for ( s_fs=biggest(mask_fs-pws/2,
			                   p->min_fs);
			      s_fs<=smallest(mask_fs+pws/2,
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
		                   &pbg, &pmax, &sigma, 1, 1);

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

		if ( fabs(intensity)/sigma < min_snr ) {
			nrej_snr++;
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
//	       "%i in bad regions, %i with SNR < %g, %i badrow culled.\n",
//	       nacc, nrej_dis, nrej_pro, nrej_fra, nrej_bad,
//	       nrej_snr, min_snr, ncull);
}


void search_peaks(struct image *image, float threshold, float min_gradient,
                  float min_snr)
{
	int i;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	for ( i=0; i<image->det->n_panels; i++ ) {

		struct panel *p = &image->det->panels[i];

		if ( p->no_index ) continue;
		search_peaks_in_panel(image, threshold, min_gradient,
		                      min_snr, p);

	}
}


double peak_lattice_agreement(struct image *image, UnitCell *cell, double *pst)
{
	int i;
	int n_feat = 0;
	int n_sane = 0;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	double min_dist = 0.25;
	double stot = 0.0;

	/* Round towards nearest */
	fesetround(1);

	/* Cell basis vectors for this image */
	cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	/* Loop over peaks, checking proximity to nearest reflection */
	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec q;
		double h,k,l,hd,kd,ld;

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		n_feat++;

		/* Reciprocal space position of found peak */
		q = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);

		/* Decimal and fractional Miller indices of nearest
		 * reciprocal lattice point */
		hd = q.u * ax + q.v * ay + q.w * az;
		kd = q.u * bx + q.v * by + q.w * bz;
		ld = q.u * cx + q.v * cy + q.w * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		/* Check distance */
		if ( (fabs(h - hd) < min_dist) && (fabs(k - kd) < min_dist)
		  && (fabs(l - ld) < min_dist) )
		{
			double sval;
			n_sane++;
			sval = pow(h-hd, 2.0) + pow(k-kd, 2.0) + pow(l-ld, 2.0);
			stot += 1.0 - sval;
			continue;
		}

	}

	*pst = stot;
	return (double)n_sane / (float)n_feat;
}


int peak_sanity_check(struct image *image)
{
	double stot;
	/* 0 means failed test, 1 means passed test */
	return peak_lattice_agreement(image, image->indexed_cell, &stot) >= 0.5;
}


struct integr_ind
{
	signed int h;
	signed int k;
	signed int l;
	double res;
	Reflection *refl;
};


static int compare_resolution(const void *av, const void *bv)
{
	const struct integr_ind *a = av;
	const struct integr_ind *b = bv;

	return a->res > b->res;
}

static struct integr_ind *sort_reflections(RefList *list, UnitCell *cell,
                                           int *np)
{
	struct integr_ind *il;
	Reflection *refl;
	RefListIterator *iter;
	int i, n;

	n = num_reflections(list);
	*np = 0;  /* For now */

	if ( n == 0 ) return NULL;

	il = calloc(n, sizeof(struct integr_ind));
	if ( il == NULL ) return NULL;

	i = 0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		double res;

		get_indices(refl, &h, &k, &l);
		res = resolution(cell, h, k, l);

		il[i].h = h;
		il[i].k = k;
		il[i].l = l;
		il[i].res = res;
		il[i].refl = refl;

		i++;
		assert(i <= n);
	}

	qsort(il, n, sizeof(struct integr_ind), compare_resolution);

	*np = n;
	return il;
}


/* Integrate the list of predicted reflections in "image" */
void integrate_reflections(struct image *image, int use_closer, int bgsub,
                           double min_snr)
{
	struct integr_ind *il;
	int n, i;
	double av = 0.0;
	int first = 1;

	il = sort_reflections(image->reflections, image->indexed_cell, &n);
	if ( il == NULL ) {
		ERROR("Couldn't sort reflections\n");
		return;
	}

	for ( i=0; i<n; i++ ) {

		double fs, ss, intensity;
		double d;
		int idx;
		double bg, max;
		double sigma, snr;
		double pfs, pss;
		int r;
		Reflection *refl;

		refl = il[i].refl;

		get_detector_pos(refl, &pfs, &pss);

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
		                   &intensity, &bg, &max, &sigma, 0, bgsub);

		/* Record intensity and set redundancy to 1 on success */
		if ( r == 0 ) {
			set_intensity(refl, intensity);
			set_esd_intensity(refl, sigma);
			set_redundancy(refl, 1);
		} else {
			set_redundancy(refl, 0);
		}

		snr = intensity / sigma;
		if ( snr > 1.0 ) {
			if ( first ) {
				av = snr;
				first = 0;
			} else {
				av = av + 0.1*(snr - av);
			}
			//STATUS("%5.2f A, %5.2f, av %5.2f\n",
			//       1e10/il[i].res, snr, av);
			//if ( av < 1.0 ) break;
		}
	}

	free(il);
}
