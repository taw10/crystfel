/*
 * peaks.h
 *
 * Peak search and other image analysis
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2020 Thomas White <taw@physics.org>
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

#ifndef PEAKS_H
#define PEAKS_H

#include <pthread.h>

#include "reflist.h"
#include "crystal.h"
#include "image.h"
#include "detgeom.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \file peaks.h
 * Peak search functions
 */

enum peak_search_method {
	PEAK_PEAKFINDER9,
	PEAK_PEAKFINDER8,
	PEAK_ZAEF,
	PEAK_HDF5,
	PEAK_CXI,
	PEAK_MSGPACK,
	PEAK_NONE,
	PEAK_ERROR
};


struct peak_params {

	enum peak_search_method method;

	float threshold;                /* zaef, pf8 */
	float min_sq_gradient;          /* zaef */
	float min_snr;                  /* zaef, pf8 */
	int min_pix_count;              /* pf8 */
	int max_pix_count;              /* pf8 */
	int local_bg_radius;            /* pf8 */
	int min_res;                    /* pf8 */
	int max_res;                    /* pf8 */
	int peakfinder8_fast;           /* pf8 */

	float min_snr_biggest_pix;      /* pf9 */
	float min_snr_peak_pix;         /* pf9 */
	float min_sig;                  /* pf9 */
	float min_peak_over_neighbour;  /* pf9 */

	float pk_inn;
	float pk_mid;
	float pk_out;
	int half_pixel_shift;           /* cxi, hdf5 */
	int revalidate;
	int noisefilter;
	int median_filter;
	int check_hdf5_snr;
	int use_saturated;
};


extern const char *str_peaksearch(enum peak_search_method meth);

extern enum peak_search_method parse_peaksearch(const char *arg);

extern int *make_BgMask(struct image *image, struct detgeom_panel *p,
                        int pn, double ir_inn);

extern ImageFeatureList *search_peaks(const struct image *image, float threshold,
                                      float min_gradient, float min_snr, double ir_inn,
                                      double ir_mid, double ir_out, int use_saturated);

extern ImageFeatureList *search_peaks_peakfinder9(const struct image *image,
                                                  float min_snr_biggest_pix,
                                                  float min_snr_peak_pix,
                                                  float min_snr_whole_peak, float min_sig,
                                                  float min_peak_over_neighbour,
                                                  int window_radius);

extern int indexing_peak_check(const struct image *image, ImageFeatureList *peaks,
                               Crystal **crystals,
                               int n_cryst, int multi_mode);

extern ImageFeatureList *validate_peaks(const struct image *image,
                                        ImageFeatureList *peaks,
                                        double min_snr,
                                        int ir_inn, int ir_mid, int ir_out,
                                        int use_saturated, int check_snr);

extern double estimate_peak_resolution(ImageFeatureList *peaks,
                                       double lambda,
                                       struct detgeom *det);

#ifdef __cplusplus
}
#endif

#endif	/* PEAKS_H */
