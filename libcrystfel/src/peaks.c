/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
 *   2012      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011      Andrew Martin <andrew.martin@desy.de>
 *   2011      Richard Kirian
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
#include "cell-utils.h"
#include "geometry.h"


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


static void add_crystal_to_mask(struct image *image, struct panel *p,
                                double ir_inn, int w, int h,
                                int *mask, Crystal *cr)
{
	Reflection *refl;
	RefListIterator *iter;

	/* Loop over all reflections */
	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		struct panel *p2;
		double pk2_fs, pk2_ss;
		signed int dfs, dss;
		double pk2_cfs, pk2_css;

		get_detector_pos(refl, &pk2_fs, &pk2_ss);

		/* Determine if reflection is in the same panel */
		p2 = find_panel(image->det, pk2_fs, pk2_ss);
		if ( p2 != p ) continue;

		pk2_cfs = pk2_fs - p->min_fs;
		pk2_css = pk2_ss - p->min_ss;

		for ( dfs=-ir_inn; dfs<=ir_inn; dfs++ ) {
		for ( dss=-ir_inn; dss<=ir_inn; dss++ ) {

			signed int fs, ss;

			/* In peak region for this peak? */
			if ( dfs*dfs + dss*dss > ir_inn*ir_inn ) continue;

			fs = pk2_cfs + dfs;
			ss = pk2_css + dss;

			/* On panel? */
			if ( fs >= w ) continue;
			if ( ss >= h ) continue;
			if ( fs < 0 ) continue;
			if ( ss < 0 ) continue;

			mask[fs + ss*w]++;

		}
		}

	}
}


/* cfs, css relative to panel origin */
int *make_BgMask(struct image *image, struct panel *p, double ir_inn)
{
	int *mask;
	int w, h;
	int i;

	w = p->max_fs - p->min_fs + 1;
	h = p->max_ss - p->min_ss + 1;
	mask = calloc(w*h, sizeof(int));
	if ( mask == NULL ) return NULL;

	if ( image->crystals == NULL ) return mask;

	for ( i=0; i<image->n_crystals; i++ ) {
		add_crystal_to_mask(image, p, ir_inn,
		                    w, h, mask, image->crystals[i]);
	}

	return mask;
}


/* Returns non-zero if peak has been vetoed.
 * i.e. don't use result if return value is not zero. */
static int integrate_peak(struct image *image, int cfs, int css,
                          double *pfs, double *pss,
                          double *intensity, double *sigma,
                          double ir_inn, double ir_mid, double ir_out,
                          int *saturated)
{
	signed int dfs, dss;
	double lim_sq, out_lim_sq, mid_lim_sq;
	double pk_total;
	int pk_counts;
	double fsct, ssct;
	double bg_tot = 0.0;
	int bg_counts = 0;
	struct panel *p;
	double bg_mean, bg_var;
	double bg_tot_sq = 0.0;
	double var;
	double aduph;
	int p_cfs, p_css, p_w, p_h;
	signed int pn;

	pn = find_panel_number(image->det, cfs, css);
	if ( pn == -1 ) return 2;
	p = &image->det->panels[pn];

	if ( saturated != NULL ) *saturated = 0;

	/* Determine regions where there is expected to be a peak */
	p_cfs = cfs - p->min_fs;
	p_css = css - p->min_ss;  /* Panel-relative coordinates */
	p_w = p->max_fs - p->min_fs + 1;
	p_h = p->max_ss - p->min_ss + 1;

	aduph = p->adu_per_eV * ph_lambda_to_eV(image->lambda);

	lim_sq = pow(ir_inn, 2.0);
	mid_lim_sq = pow(ir_mid, 2.0);
	out_lim_sq = pow(ir_out, 2.0);

	/* Estimate the background */
	for ( dfs=-ir_out; dfs<=+ir_out; dfs++ ) {
	for ( dss=-ir_out; dss<=+ir_out; dss++ ) {

		double val;
		int idx;

		/* Restrict to annulus */
		if ( dfs*dfs + dss*dss > out_lim_sq ) continue;
		if ( dfs*dfs + dss*dss < mid_lim_sq ) continue;

		/* Strayed off one panel? */
		if ( (p_cfs+dfs >= p_w) || (p_css+dss >= p_h)
		  || (p_cfs+dfs < 0 ) || (p_css+dss < 0) ) return 4;

		/* Wandered into a bad region? */
		if ( image->bad[pn][p_cfs+dfs + p->w*(p_css+dss)] ) {
			return 14;
		}

		idx = dfs+cfs+image->width*(dss+css);
		val = image->data[idx];

		/* Check if peak contains saturation in bg region */
		if ( (saturated != NULL) && (val > p->max_adu) ) *saturated = 1;

		bg_tot += val;
		bg_tot_sq += pow(val, 2.0);
		bg_counts++;

	}
	}

	if ( bg_counts == 0 ) return 7;
	bg_mean = bg_tot / bg_counts;
	bg_var = (bg_tot_sq/bg_counts) - pow(bg_mean, 2.0);

	/* Measure the peak */
	pk_total = 0.0;
	pk_counts = 0;
	fsct = 0.0;  ssct = 0.0;
	for ( dfs=-ir_inn; dfs<=+ir_inn; dfs++ ) {
	for ( dss=-ir_inn; dss<=+ir_inn; dss++ ) {

		double val;
		int idx;

		/* Inner mask radius */
		if ( dfs*dfs + dss*dss > lim_sq ) continue;

		/* Strayed off one panel? */
		if ( (p_cfs+dfs >= p_w) || (p_css+dss >= p_h)
		  || (p_cfs+dfs < 0 ) || (p_css+dss < 0) ) return 8;

		/* Wandered into a bad region? */
		if ( image->bad[pn][p_cfs+dfs + p->w*(p_css+dss)] ) {
			return 15;
		}

		idx = dfs+cfs+image->width*(dss+css);
		val = image->data[idx];

		/* Check if peak contains saturation */
		if ( (saturated != NULL) && (val > p->max_adu) ) *saturated = 1;

		val -= bg_mean;

		pk_counts++;
		pk_total += val;

		fsct += val*(cfs+dfs);
		ssct += val*(css+dss);

	}
	}

	if ( pk_counts == 0 ) return 11;

	*pfs = ((double)fsct / pk_total) + 0.5;
	*pss = ((double)ssct / pk_total) + 0.5;

	var = pk_counts * bg_var;
	var += aduph * pk_total;
	if ( var < 0.0 ) return 12;

	if ( intensity != NULL ) *intensity = pk_total;
	if ( sigma != NULL ) *sigma = sqrt(var);

	return 0;
}


static void search_peaks_in_panel(struct image *image, float threshold,
                                  float min_gradient, float min_snr,
                                  struct panel *p,
                                  double ir_inn, double ir_mid, double ir_out,
                                  int use_saturated)
{
	int fs, ss, stride;
	float *data;
	double d;
	int idx;
	double f_fs = 0.0;
	double f_ss = 0.0;
	double intensity = 0.0;
	double sigma = 0.0;
	int nrej_dis = 0;
	int nrej_pro = 0;
	int nrej_fra = 0;
	int nrej_fail = 0;
	int nrej_snr = 0;
	int nrej_sat = 0;
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
		int saturated;

		/* Overall threshold */
		if ( data[fs+stride*ss] < threshold ) continue;

		/* Immediate rejection of pixels above max_adu */
		if ( !use_saturated && (data[fs+stride*ss] > p->max_adu) ) {
			continue;
		}

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

			for ( s_ss=biggest(mask_ss-ir_inn, p->min_ss);
			      s_ss<=smallest(mask_ss+ir_inn, p->max_ss);
			      s_ss++ )
			{
			for ( s_fs=biggest(mask_fs-ir_inn, p->min_fs);
			      s_fs<=smallest(mask_fs+ir_inn, p->max_fs);
			      s_fs++ )
			{

				if ( data[s_fs+stride*s_ss] > max ) {
					max = data[s_fs+stride*s_ss];
					mask_fs = s_fs;
					mask_ss = s_ss;
					did_something = 1;
				}

			}
			}

			/* Abort if drifted too far from the foot point */
			if ( distance(mask_fs, mask_ss, fs, ss) > ir_inn )
			{
				break;
			}

		} while ( did_something );

		/* Too far from foot point? */
		if ( distance(mask_fs, mask_ss, fs, ss) > ir_inn ) {
			nrej_dis++;
			continue;
		}

		/* Should be enforced by bounds used above.  Muppet check. */
		assert(mask_fs <= p->max_fs);
		assert(mask_ss <= p->max_ss);
		assert(mask_fs >= p->min_fs);
		assert(mask_ss >= p->min_ss);

		/* Centroid peak and get better coordinates. */
		r = integrate_peak(image, mask_fs, mask_ss,
		                   &f_fs, &f_ss, &intensity, &sigma,
		                   ir_inn, ir_mid, ir_out, &saturated);

		if ( r ) {
			/* Bad region - don't detect peak */
			nrej_fail++;
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
		image_feature_closest(image->features, f_fs, f_ss, &d, &idx,
		                      image->det);
		if ( d < 2.0*ir_inn ) {
			nrej_pro++;
			continue;
		}

		if ( saturated ) {
			image->num_saturated_peaks++;
			if ( !use_saturated ) {
				nrej_sat++;
				continue;
			}
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

	image->num_peaks += nacc;

	//STATUS("%i accepted, %i box, %i proximity, %i outside panel, "
	//       "%i failed integration, %i with SNR < %g, %i badrow culled, "
	//        "%i saturated.\n",
	//       nacc, nrej_dis, nrej_pro, nrej_fra, nrej_fail,
	//       nrej_snr, min_snr, ncull, nrej_sat);

	if ( ncull != 0 ) {
		STATUS("WARNING: %i peaks were badrow culled.  This feature"
		       " should not usually be used.\nConsider setting"
		       " badrow=- in the geometry file.\n", ncull);
	}
}


void search_peaks(struct image *image, float threshold, float min_gradient,
                  float min_snr, double ir_inn, double ir_mid, double ir_out,
                  int use_saturated)
{
	int i;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();
	image->num_peaks = 0;
	image->num_saturated_peaks = 0;

	for ( i=0; i<image->det->n_panels; i++ ) {

		struct panel *p = &image->det->panels[i];

		if ( p->no_index ) continue;
		search_peaks_in_panel(image, threshold, min_gradient,
		                      min_snr, p, ir_inn, ir_mid, ir_out,
		                      use_saturated);

	}
}


int peak_sanity_check(struct image *image, Crystal **crystals, int n_cryst)
{
	int n_feat = 0;
	int n_sane = 0;
	int i;
	const double min_dist = 0.25;

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		struct rvec q;
		double h,k,l,hd,kd,ld;
		int j;

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		n_feat++;

		/* Reciprocal space position of found peak */
		q = get_q(image, f->fs, f->ss, NULL, 1.0/image->lambda);

		for ( j=0; j<n_cryst; j++ ) {

			double ax, ay, az;
			double bx, by, bz;
			double cx, cy, cz;

			cell_get_cartesian(crystal_get_cell(crystals[j]),
			                   &ax, &ay, &az,
			                   &bx, &by, &bz,
			                   &cx, &cy, &cz);

			/* Decimal and fractional Miller indices of nearest
			 * reciprocal lattice point */
			hd = q.u * ax + q.v * ay + q.w * az;
			kd = q.u * bx + q.v * by + q.w * bz;
			ld = q.u * cx + q.v * cy + q.w * cz;
			h = lrint(hd);
			k = lrint(kd);
			l = lrint(ld);

			/* Check distance */
			if ( (fabs(h - hd) < min_dist)
			  && (fabs(k - kd) < min_dist)
			  && (fabs(l - ld) < min_dist) )
			{
				n_sane++;
				continue;
			}

		}


	}

	/* 0 means failed test, 1 means passed test */
	return ((double)n_sane / n_feat) >= 0.5;
}


void validate_peaks(struct image *image, double min_snr,
                    int ir_inn, int ir_mid, int ir_out, int use_saturated,
                    int check_snr)
{
	int i, n;
	ImageFeatureList *flist;
	int n_wtf, n_int, n_dft, n_snr, n_prx, n_sat;

	flist = image_feature_list_new();
	if ( flist == NULL ) return;

	n = image_feature_count(image->features);

	/* Loop over peaks, putting each one through the integrator */
	n_wtf = 0;  n_int = 0;  n_dft = 0;  n_snr = 0;  n_prx = 0;  n_sat = 0;
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		int r;
		double d;
		int idx;
		double f_fs, f_ss;
		double intensity, sigma;
		struct panel *p;
		int saturated;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) {
			n_wtf++;
			continue;
		}

		p = find_panel(image->det, f->fs, f->ss);
		if ( p == NULL ) {
			n_wtf++;
			continue;
		}

		r = integrate_peak(image, f->fs, f->ss,
		                   &f_fs, &f_ss, &intensity, &sigma,
		                   ir_inn, ir_mid, ir_out, &saturated);
		if ( r ) {
			n_int++;
			continue;
		}

		if ( saturated ) {
			if ( !use_saturated ) {
				n_sat++;
				continue;
			}
		}

		/* It is possible for the centroid to fall outside the image */
		if ( (f_fs < p->min_fs) || (f_fs > p->max_fs)
		  || (f_ss < p->min_ss) || (f_ss > p->max_ss) )
		{
			n_dft++;
			continue;
		}

		if ( check_snr && (fabs(intensity)/sigma < min_snr) ) {
			n_snr++;
			continue;
		}

		/* Check for a nearby feature */
		image_feature_closest(flist, f_fs, f_ss, &d, &idx, image->det);
		if ( d < 2.0*ir_inn ) {
			n_prx++;
			continue;
		}

		/* Add using "better" coordinates */
		image_add_feature(flist, f_fs, f_ss, image, intensity, NULL);

	}

	//STATUS("HDF5: %i peaks, validated: %i.  WTF: %i, integration: %i,"
	//       " drifted: %i, SNR: %i, proximity: %i, saturated: %i\n",
	//       n, image_feature_count(flist),
	//       n_wtf, n_int, n_dft, n_snr, n_prx, n_sat);
	image_feature_list_free(image->features);
	image->features = flist;
	image->num_saturated_peaks = n_sat;
	image->num_peaks = image_feature_count(flist);
}
