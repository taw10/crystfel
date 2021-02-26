/*
 * robustpeakfinder.c
 *
 * support for Robust Peak Finder for CrystFEL
 * from https://github.com/MarjanHJ/RobustPeakFinder
 * which is based on
 * https://github.com/ARSadri/RobustGaussianFittingLibrary
 * Please cite http://scripts.iucr.org/cgi-bin/paper?S1600576717014340
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 * 	 2017  Marjan Hadian-Jazi <marjan.hadian-jazi@xfel.eu>
 * 	 2020  Alireza Sadri <Alireza.Sadri@desy.de>
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

#ifndef robustpeakfinder_H
#define robustpeakfinder_H

#include "image.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \param img An \ref image structure
 * \param max_n_peaks The maximum number of peaks to be searched for
 * \param detectorsNoiseSTD for photon counting detectors, 1
 * \param min_snr The minimum signal to noise ratio for a peak
 * \param min_pix_count The minimum number of pixels in a peak
 * \param max_pix_count The maximum number of pixels in a peak
 * \param local_bg_radius The averaging radius for background calculation
 * \param min_res The minimum number of pixels out from the center
 * \param max_res The maximum number of pixels out from the center
 * \param use_saturated Whether saturated peaks should be considered
 * \param singlePhotonADU The image singlePhotonADU value, in detector units
 * \param supportGradient use a tilting plane 1 to model the background
 * \param inlier_SNR use this snr far to model the background noise
 * \param search_SNR search for Bragg peaks after this snr
 * \param finiteSampleBias 200 pixels are neded to estimate scale
 * \param n_optIters number of iterations of the model fitting optimization
 * \param topKthPerc the rough guess of portion of data that are inliers, 0.5
 * \param botKthPerc the rough guess of poetion of data that are good inliers 0.3
 * \param maxBackMeanMap the maximum value of the backgournd average
 * \param downSampledSize 200 data points are needed to estimate average.
 * \param highPoissonTh the maximum of acceptable average/variance of background
 * \param lowPoissonTh the minimum of acceptable average/variance of background
 */
extern int robustpeakfinder(struct image *img, 
							int    max_n_peaks,
							float  detectorsNoiseSTD, 
							float  min_snr,
							int    min_pix_count, 
							int    max_pix_count,
							int    local_bg_radius, 
							int    min_res,
							int    max_res, 
							int    use_saturated,
							int    supportGradient,
							float  inlier_SNR,
							float  search_SNR,
							int    finiteSampleBias,
							int    n_optIters,
							float  topKthPerc,
							float  botKthPerc,
							float  maxBackMeanMap,
							int    downSampledSize,
							float  highPoissonTh,
							float  lowPoissonTh);
					   
#ifdef __cplusplus     
}
#endif

#endif // robustpeakfinder_H
