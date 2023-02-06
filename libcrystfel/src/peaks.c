/*
 * peaks.c
 *
 * Peak search and other image analysis
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
 *   2012      Kenneth Beyerlein <kenneth.beyerlein@desy.de>
 *   2011      Andrew Martin <andrew.martin@desy.de>
 *   2011      Richard Kirian
 *   2017      Valerio Mariani <valerio.mariani@desy.de>
 *   2017-2018 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
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


#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_statistics_int.h>
#include <pthread.h>
#include <fenv.h>

#ifdef HAVE_FDIP
#include "fastDiffractionImageProcessing/adaptions/crystfel/peakFinder9.h"
#include "fastDiffractionImageProcessing/adaptions/crystfel/mask.h"
#include "fastDiffractionImageProcessing/peakList.h"
#endif

#include "image.h"
#include "utils.h"
#include "peaks.h"
#include "detgeom.h"
#include "filters.h"
#include "reflist-utils.h"
#include "cell-utils.h"
#include "geometry.h"
#include "peakfinder8.h"

/** \file peaks.h */

static void add_crystal_to_mask(struct image *image,
                                struct detgeom_panel *p, int pn,
                                double ir_inn, int *mask, Crystal *cr)
{
	Reflection *refl;
	RefListIterator *iter;

	/* Loop over all reflections */
	for ( refl = first_refl(crystal_get_reflections(cr), &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pk2_fs, pk2_ss;
		signed int dfs, dss;

		get_detector_pos(refl, &pk2_fs, &pk2_ss);

		/* Determine if reflection is in the same panel */
		if ( get_panel_number(refl) != pn ) continue;

		for ( dfs=-ir_inn; dfs<=ir_inn; dfs++ ) {
		for ( dss=-ir_inn; dss<=ir_inn; dss++ ) {

			signed int fs, ss;

			/* In peak region for this peak? */
			if ( dfs*dfs + dss*dss > ir_inn*ir_inn ) continue;

			fs = pk2_fs + dfs;
			ss = pk2_ss + dss;

			/* On panel? */
			if ( fs >= p->w ) continue;
			if ( ss >= p->h ) continue;
			if ( fs < 0 ) continue;
			if ( ss < 0 ) continue;

			mask[fs + ss*p->w]++;

		}
		}

	}
}


/* cfs, css relative to panel origin */
int *make_BgMask(struct image *image, struct detgeom_panel *p,
                 int pn, double ir_inn)
{
	int *mask;
	int i;

	mask = calloc(p->w*p->h, sizeof(int));
	if ( mask == NULL ) return NULL;

	if ( image->crystals == NULL ) return mask;

	for ( i=0; i<image->n_crystals; i++ ) {
		add_crystal_to_mask(image, p, pn, ir_inn,
		                    mask, image->crystals[i]);
	}

	return mask;
}


/* Returns non-zero if peak has been vetoed.
 * i.e. don't use result if return value is not zero. */
int integrate_peak(struct image *image,
                   int p_cfs, int p_css, int pn,
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
	double bg_mean, bg_var;
	double bg_tot_sq = 0.0;
	double var;
	double aduph, bias;
	struct detgeom_panel *p;

	if ( saturated != NULL ) *saturated = 0;

	p = &image->detgeom->panels[pn];

	aduph = p->adu_per_photon;
	bias = p->adu_bias;

	lim_sq = pow(ir_inn, 2.0);
	mid_lim_sq = pow(ir_mid, 2.0);
	out_lim_sq = pow(ir_out, 2.0);

	/* Estimate the background */
	for ( dss=-ir_out; dss<=+ir_out; dss++ ) {
	for ( dfs=-ir_out; dfs<=+ir_out; dfs++ ) {

		double val;
		int idx;

		/* Restrict to annulus */
		if ( dfs*dfs + dss*dss > out_lim_sq ) continue;
		if ( dfs*dfs + dss*dss < mid_lim_sq ) continue;

		/* Strayed off one panel? */
		if ( (p_cfs+dfs >= p->w) || (p_css+dss >= p->h)
		  || (p_cfs+dfs < 0 ) || (p_css+dss < 0) ) return 4;

		/* Wandered into a bad region? */
		if ( image->bad[pn][p_cfs+dfs + p->w*(p_css+dss)] ) {
			return 14;
		}

		idx = dfs+p_cfs+p->w*(dss+p_css);
		val = image->dp[pn][idx];

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
	for ( dss=-ir_inn; dss<=+ir_inn; dss++ ) {
	for ( dfs=-ir_inn; dfs<=+ir_inn; dfs++ ) {

		double val;
		int idx;

		/* Inner mask radius */
		if ( dfs*dfs + dss*dss > lim_sq ) continue;

		/* Strayed off one panel? */
		if ( (p_cfs+dfs >= p->w) || (p_css+dss >= p->h)
		  || (p_cfs+dfs < 0 ) || (p_css+dss < 0) ) return 8;

		/* Wandered into a bad region? */
		if ( image->bad[pn][p_cfs+dfs + p->w*(p_css+dss)] ) {
			return 15;
		}

		idx = dfs+p_cfs+p->w*(dss+p_css);
		val = image->dp[pn][idx];

		/* correct raw intensity for the bias */
		val -= bias;

		/* Check if peak contains saturation */
		if ( (saturated != NULL) && (val > p->max_adu) ) *saturated = 1;

		val -= bg_mean;

		pk_counts++;
		pk_total += val;

		fsct += val*(p_cfs+dfs);
		ssct += val*(p_css+dss);

	}
	}

	if ( pk_counts == 0 ) return 11;
	if ( pk_total == 0 ) return 13;

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
                                  float min_sq_gradient, float min_snr, int pn,
                                  double ir_inn, double ir_mid, double ir_out,
                                  int use_saturated)
{
	int fs, ss, stride;
	float *data;
	struct detgeom_panel *p;
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

	p = &image->detgeom->panels[pn];
	data = image->dp[pn];
	stride = p->w;

	for ( ss=1; ss<p->h-1; ss++ ) {
	for ( fs=1; fs<p->w-1; fs++ ) {

		double dx1, dx2, dy1, dy2;
		double dxs, dys;
		double grad;
		int mask_fs, mask_ss;
		int s_fs, s_ss;
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

		/* Calculate overall (squared) gradient */
		grad = dxs + dys;

		if ( grad < min_sq_gradient ) continue;

		mask_fs = fs;
		mask_ss = ss;

		do {

			double max;
			max = data[mask_fs+stride*mask_ss];
			did_something = 0;

			for ( s_ss=biggest(mask_ss-ir_inn, 0);
			      s_ss<=smallest(mask_ss+ir_inn, p->h-1);
			      s_ss++ )
			{
			for ( s_fs=biggest(mask_fs-ir_inn, 0);
			      s_fs<=smallest(mask_fs+ir_inn, p->w-1);
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
		assert(mask_fs <= p->w);
		assert(mask_ss <= p->h);
		assert(mask_fs >= 0);
		assert(mask_ss >= 0);

		/* Centroid peak and get better coordinates. */
		r = integrate_peak(image, mask_fs, mask_ss, pn,
		                   &f_fs, &f_ss, &intensity, &sigma,
		                   ir_inn, ir_mid, ir_out, &saturated);

		if ( r ) {
			/* Bad region - don't detect peak */
			nrej_fail++;
			continue;
		}

		/* It is possible for the centroid to fall outside the image */
		if ( (f_fs < 0) || (f_fs > p->w)
		  || (f_ss < 0) || (f_ss > p->h) ) {
			nrej_fra++;
			continue;
		}

		if ( fabs(intensity)/sigma < min_snr ) {
			nrej_snr++;
			continue;
		}

		/* Check for a nearby feature */
		image_feature_closest(image->features, f_fs, f_ss, pn,
		                      &d, &idx);
		if ( d < 2.0*ir_inn ) {
			nrej_pro++;
			continue;
		}

		if ( saturated && !use_saturated ) {
			nrej_sat++;
			continue;
		}

		/* Add using "better" coordinates */
		image_add_feature(image->features, f_fs, f_ss, pn,
		                  image, intensity, NULL);
		nacc++;

		if ( nacc > 10000 ) {
			ERROR("Too many peaks!  Aborting peak seach "
			      "for panel %s\n", p->name);
			return;
		}

	}
	}

	//STATUS("%i accepted, %i box, %i proximity, %i outside panel, "
	//       "%i failed integration, %i with SNR < %g, %i badrow culled, "
	//        "%i saturated.\n",
	//       nacc, nrej_dis, nrej_pro, nrej_fra, nrej_fail,
	//       nrej_snr, min_snr, nrej_sat);

}


void search_peaks(struct image *image, float threshold, float min_sq_gradient,
                  float min_snr, double ir_inn, double ir_mid,
                  double ir_out, int use_saturated)
{
	int i;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	for ( i=0; i<image->detgeom->n_panels; i++ ) {

		search_peaks_in_panel(image, threshold, min_sq_gradient,
		                      min_snr, i, ir_inn, ir_mid, ir_out,
		                      use_saturated);

	}
}


/**
 * \param image An \ref image structure
 * \param max_n_peaks The maximum number of peaks to be searched for
 * \param threshold The image threshold value, in detector units
 * \param min_snr The minimum signal to noise ratio for a peak
 * \param min_pix_count The minimum number of pixels in a peak
 * \param max_pix_count The maximum number of pixels in a peak
 * \param local_bg_radius The averaging radius for background calculation
 * \param min_res The minimum number of pixels out from the center
 * \param max_res The maximum number of pixels out from the center
 * \param use_saturated Whether saturated peaks should be considered
 *
 * Runs the peakfinder8 peak search algorithm.  This is a thin wrapper which
 * creates an empty \ref ImageFeatureList for \p image, then calls
 * the actual \ref peakfinder8 function, found in \ref peakfinder8.h.
 */
int search_peaks_peakfinder8(struct image *image, int max_n_peaks,
                             float threshold, float min_snr,
                             int min_pix_count, int max_pix_count,
                             int local_bg_radius, int min_res,
                             int max_res, int use_saturated,
                             int fast_mode, void *private_data)
{
	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	return peakfinder8(image, max_n_peaks, threshold, min_snr,
	                   min_pix_count, max_pix_count,
	                   local_bg_radius, min_res,
	                   max_res, use_saturated,
			   fast_mode, private_data);
}


#ifdef HAVE_FDIP

int search_peaks_peakfinder9(struct image *image, float min_snr_biggest_pix,
                             float min_snr_peak_pix, float min_snr_whole_peak,
                             float min_sig, float min_peak_over_neighbour,
                             int window_radius)
{
	peakFinder9_accuracyConstants_t accuracy_consts;
	peakList_t peakList;
	long NpeaksMax = 10000; //more peaks per panel should not appear
	float *data_copy = NULL;
	float *data_copy_new;
	int panel_number;

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();

	accuracy_consts.minSNR_biggestPixel = min_snr_biggest_pix;
	accuracy_consts.minSNR_peakPixel = min_snr_peak_pix;
	accuracy_consts.minSNR_wholePeak = min_snr_whole_peak;
	accuracy_consts.minimumSigma = min_sig;
	accuracy_consts.minimumPeakOversizeOverNeighbours = min_peak_over_neighbour;
	accuracy_consts.windowRadius = window_radius;

	if ( allocatePeakList(&peakList, NpeaksMax) ) return 1;

	for ( panel_number=0; panel_number<image->detgeom->n_panels; panel_number++ ) {

		int w, h;
		int peak_number;
		detectorRawFormat_t det_size_one_panel;

		w = image->detgeom->panels[panel_number].w;
		h = image->detgeom->panels[panel_number].h;

		det_size_one_panel.asic_nx = w;
		det_size_one_panel.asic_ny = h;
		det_size_one_panel.nasics_x = 1;
		det_size_one_panel.nasics_y = 1;
		det_size_one_panel.pix_nx = w;
		det_size_one_panel.pix_ny = h;
		det_size_one_panel.pix_nn = w * h;

		data_copy_new = realloc(data_copy, w*h*sizeof(*data_copy));
		if ( data_copy_new == NULL ) {
			if ( data_copy != NULL ) {
				free(data_copy);
			}
			freePeakList(peakList);
			return 1;
		} else {
			data_copy = data_copy_new;
		}

		mergeMaskAndDataIntoDataCopy(image->dp[panel_number], data_copy,
		                             image->bad[panel_number],
		                             &det_size_one_panel);

		peakList.peakCount = 0;
		peakFinder9_onePanel_noSlab(data_copy, &accuracy_consts,
		                            &det_size_one_panel, &peakList);

		for ( peak_number=0; peak_number<peakList.peakCount; peak_number++) {
			image_add_feature(image->features,
			                  peakList.centerOfMass_rawX[peak_number],
			                  peakList.centerOfMass_rawY[peak_number],
			                  panel_number, image,
			                  peakList.totalIntensity[peak_number],
			                  NULL);
		}

	}

	freePeakList(peakList);
	free(data_copy);
	return 0;
}

#else

int search_peaks_peakfinder9(struct image *image, float min_snr_biggest_pix,
                             float min_snr_peak_pix, float min_snr_whole_peak,
                             float min_sig, float min_peak_over_neighbour,
                             int window_radius)
{
	ERROR("This copy of CrystFEL was compiled without peakfinder9 support.\n");
	return 1;
}

#endif // HAVE_FDIP


/**
 * \param image An \ref image structure
 * \param crystals Pointer to array of pointers to crystals
 * \param n_cryst The number of crystals
 * \param multi_mode Whether the thresholds should be set for multi-lattice indexing
 *
 * Checks whether the peaks in \p image appear to be explained by the crystals
 * provided.
 *
 * Returns 1 if the peaks appear to be well-explained by the crystals.
 * Otherwise, if the indexing solutions appear to be "bad", returns 0.
 */
int indexing_peak_check(struct image *image, Crystal **crystals, int n_cryst,
                        int multi_mode)
{
	int n_feat = 0;
	int n_sane = 0;
	int i;
	const double min_dist = 0.25;

	for ( i=0; i<image_feature_count(image->features); i++ ) {

		struct imagefeature *f;
		double q[3];
		int j;
		int ok = 0;

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(image->features, i);
		if ( f == NULL ) continue;
		n_feat++;

		for ( j=0; j<n_cryst; j++ ) {

			double ax, ay, az;
			double bx, by, bz;
			double cx, cy, cz;
			double dx, dy;
			double h,k,l,hd,kd,ld;

			crystal_get_det_shift(crystals[j], &dx, &dy);

			/* Reciprocal space position of found peak,
			 * based on a calculation including any updates to the
			 * detector position from the refinement of the
			 * current crystal. */
			detgeom_transform_coords(&image->detgeom->panels[f->pn],
			                         f->fs, f->ss, image->lambda,
			                         dx, dy, q);

			cell_get_cartesian(crystal_get_cell(crystals[j]),
			                   &ax, &ay, &az,
			                   &bx, &by, &bz,
			                   &cx, &cy, &cz);

			/* Decimal and fractional Miller indices of nearest
			 * reciprocal lattice point */
			hd = q[0] * ax + q[1] * ay + q[2] * az;
			kd = q[0] * bx + q[1] * by + q[2] * bz;
			ld = q[0] * cx + q[1] * cy + q[2] * cz;
			h = lrint(hd);
			k = lrint(kd);
			l = lrint(ld);

			/* Check distance */
			if ( (fabs(h - hd) < min_dist)
			  && (fabs(k - kd) < min_dist)
			  && (fabs(l - ld) < min_dist) )
			{
				ok = 1;
				break;  /* Don't need to check other crystals */
			}

		}

		n_sane += ok;

	}

	/* 0 means failed test, 1 means passed test */

	if ( multi_mode ) {
		return (n_sane > 70)
		    || ((n_sane > 25) && (n_sane > 0.3*n_feat))
		    || (n_sane > 0.4*n_feat);
	} else {
		return ((double)n_sane / n_feat) >= 0.5;
	}
}


/**
 * Deprecated: use indexing_peak_check instead
 */
int peak_sanity_check(struct image *image, Crystal **crystals, int n_cryst)
{
	return indexing_peak_check(image, crystals, n_cryst, 1);
}


void validate_peaks(struct image *image, double min_snr,
                    int ir_inn, int ir_mid, int ir_out, int use_saturated,
                    int check_snr)
{
	int i, n;
	ImageFeatureList *flist;
	int n_wtf, n_int, n_snr, n_sat;

	flist = image_feature_list_new();
	if ( flist == NULL ) return;

	n = image_feature_count(image->features);

	/* Loop over peaks, putting each one through the integrator */
	n_wtf = 0;  n_int = 0;  n_snr = 0;  n_sat = 0;
	for ( i=0; i<n; i++ ) {

		struct imagefeature *f;
		int r;
		double f_fs, f_ss;
		double intensity, sigma;
		int saturated;

		f = image_get_feature(image->features, i);
		if ( f == NULL ) {
			n_wtf++;
			continue;
		}

		r = integrate_peak(image, f->fs, f->ss, f->pn,
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

		if ( check_snr && (fabs(intensity)/sigma < min_snr) ) {
			n_snr++;
			continue;
		}

		/* Add using "better" coordinates */
		image_add_feature(flist, f->fs, f->ss, f->pn, image,
		                  intensity, NULL);

	}

	//STATUS("HDF5: %i peaks, validated: %i.  WTF: %i, integration: %i, "
	//       "SNR: %i, saturated: %i\n",
	//       n, image_feature_count(flist), n_wtf, n_int, n_snr, n_sat);
	image_feature_list_free(image->features);
	image->features = flist;
}


double estimate_peak_resolution(ImageFeatureList *peaks, double lambda,
                                struct detgeom *det)
{
	int i, npk, ncut;
	double *rns;
	double max_res;

	npk = image_feature_count(peaks);

	/* No peaks -> no resolution! */
	if ( npk == 0 ) return 0.0;

	rns = malloc(npk*sizeof(double));
	if ( rns == NULL ) return -1.0;

	/* Get resolution values for all peaks */
	for ( i=0; i<npk; i++ ) {

		struct imagefeature *f;
		double r[3];

		f = image_get_feature(peaks, i);

		detgeom_transform_coords(&det->panels[f->pn],
		                         f->fs, f->ss,
		                         lambda, 0.0, 0.0, r);
		rns[i] = modulus(r[0], r[1], r[2]);

	}

	/* Slightly horrible outlier removal */
	qsort(rns, npk, sizeof(double), compare_double);
	ncut = npk/50;
	if ( ncut < 2 ) ncut = 0;
	max_res = rns[(npk-1)-ncut];

	free(rns);
	return max_res;
}

const char *str_peaksearch(enum peak_search_method meth)
{
	switch ( meth ) {
	case PEAK_PEAKFINDER9: return "peakfinder9";
	case PEAK_PEAKFINDER8: return "peakfinder8";
	case PEAK_ZAEF: return "zaef";
	case PEAK_HDF5: return "hdf5";
	case PEAK_CXI: return "cxi";
	case PEAK_MSGPACK: return "msgpack";
	case PEAK_NONE: return "none";
	default: return "???";
	}
}

enum peak_search_method parse_peaksearch(const char *arg)
{
	if ( strcmp(arg, "zaef") == 0 ) {
		return PEAK_ZAEF;
	} else if ( strcmp(arg, "peakfinder8") == 0 ) {
		return PEAK_PEAKFINDER8;
	} else if ( strcmp(arg, "hdf5") == 0 ) {
		return PEAK_HDF5;
	} else if ( strcmp(arg, "cxi") == 0 ) {
		return PEAK_CXI;
	} else if ( strcmp(arg, "peakfinder9") == 0 ) {
		return PEAK_PEAKFINDER9;
	} else if ( strcmp(arg, "msgpack") == 0 ) {
		return PEAK_MSGPACK;
	} else if ( strcmp(arg, "none") == 0 ) {
		return PEAK_NONE;
	}

	return PEAK_ERROR;
}
