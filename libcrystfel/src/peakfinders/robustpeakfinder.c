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
 *      2017  Marjan Hadian-Jazi <marjan.hadian-jazi@xfel.eu>
 *      2020  Alireza Sadri <Alireza.Sadri@desy.de>
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

#include <stdio.h>
#include "RPFSource.h"
#include "robustpeakfinder.h"
/** \file robustpeakfinder.h */

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
int robustpeakfinder(struct image *img, 
                     int max_n_peaks,
                     float detectorsNoiseSTD, 
                     float min_snr,
                     int min_pix_count, 
                     int max_pix_count,
                     int local_bg_radius, 
                     int min_res,
                     int max_res, 
                     int use_saturated,
                     int supportGradient,
                     float inlier_SNR,
                     float search_SNR,
                     int finiteSampleBias,
                     int n_optIters,
                     float topKthPerc,
                     float botKthPerc,
                     float maxBackMeanMap,
                     int downSampledSize,
                     float highPoissonTh,
                     float lowPoissonTh)
{
                    
    int panelCnt;
    int num_found_peaks;
    int remaining_max_num_peaks;
    char use_Mask = 1;
    int peaks_to_add;
    int pki;
    float minBackMeanMap[1];
    
	float maxBackMeanMap_pt[1];
	maxBackMeanMap_pt[0] = maxBackMeanMap;
    float peakMap[1];
    peakMap[0]=0;
    
    float singlePhotonADU;
    
    if(search_SNR >= min_snr)
        search_SNR = min_snr;
    if(inlier_SNR >= search_SNR)
        inlier_SNR = search_SNR;

    if(local_bg_radius < sqrt(finiteSampleBias)/2)
        local_bg_radius = (int)sqrt(finiteSampleBias)/2;
    
    if ( img-> detgeom == NULL) {
        return 1;
    }
    if ( img->bad == NULL) {
        use_Mask = 0;
    }
    
    remaining_max_num_peaks = max_n_peaks;

	//e.g. for AGIPD-1M, there are 128 panels each 128x64
    for ( panelCnt=0 ; panelCnt < img->detgeom->n_panels ; panelCnt++) {
	
		singlePhotonADU = img -> detgeom -> panels[panelCnt].adu_per_photon;
		if(detectorsNoiseSTD<singlePhotonADU/6.0)
			detectorsNoiseSTD = singlePhotonADU/6.0;
		minBackMeanMap[0] = 3.0*detectorsNoiseSTD;

        num_found_peaks = 0;

		int panelW = img->detgeom->panels[panelCnt].w;
		int panelH = img->detgeom->panels[panelCnt].h;
		float pcnx = img->detgeom->panels[panelCnt].cnx;
		float pcny = img->detgeom->panels[panelCnt].cny;
		float pfsx = img->detgeom->panels[panelCnt].fsx;
		float pfsy = img->detgeom->panels[panelCnt].fsy;
		float pssx = img->detgeom->panels[panelCnt].ssx;
		float pssy = img->detgeom->panels[panelCnt].ssy;
		float rm;
		int rmi;
		int x,y, iss, ifs;

        float *peakList;
        peakList = (float*) malloc(6*max_n_peaks * sizeof(float));

		unsigned char *mask;
		mask = (unsigned char*) malloc(
			panelW*panelH * sizeof(unsigned char));
			
/*			
				Height, panelH
			  slow scan AGIPD 64
			|------------------
	width	|
	fast	|
	scan	|
	AGIPD 128
	panelW
	
*/					
		for ( iss=0 ; iss<panelH ; iss++ ) {
			for ( ifs=0; ifs<panelW; ifs++ ) {

				rmi = ifs + panelW * iss;

				x = (pcnx  + ifs * pfsx + iss * pssx);
				y = (pcny  + ifs * pfsy + iss * pssy);

				rm = sqrt(x * x + y * y);

				mask[rmi] = 0;
				if ( (min_res < rm) && (rm < max_res) ) {
					if (img->bad[panelCnt][rmi] == 0) {
						mask[rmi] = 1;
					}
				}
			}
		}
		num_found_peaks = rpfMain(img->dp[panelCnt],
                                  use_Mask,
                                  mask,
                                  0,
                                  0,
                                  0,
                                  minBackMeanMap, 
                                  0,
                                  maxBackMeanMap_pt, 
                                  0,
                                  peakMap,
                                  peakList,
                                  singlePhotonADU,
                                  max_n_peaks,
                                  min_snr,
                                  2.0,
                                  panelW,
                                  panelH,
                                  local_bg_radius*2,
                                  min_pix_count, 
                                  max_pix_count,
                                  n_optIters,
                                  finiteSampleBias,
                                  downSampledSize,
                                  inlier_SNR,
                                  search_SNR,
                                  highPoissonTh,
                                  lowPoissonTh);
        
		peaks_to_add = num_found_peaks;

        if ( num_found_peaks > remaining_max_num_peaks ) {
            peaks_to_add = remaining_max_num_peaks;
        }

        remaining_max_num_peaks -= peaks_to_add;

        for ( pki=0 ; pki<peaks_to_add ; pki++ ) {

            if ( peakList[6*pki+4] > 
					img->detgeom->panels[panelCnt].max_adu ) {
                if ( !use_saturated ) {
                    continue;
                }
            }
            image_add_feature(img->features,
                              peakList[6*pki+0]+0.5,
                              peakList[6*pki+1]+0.5,
                              panelCnt,
                              img,
                              peakList[6*pki+2],
                              NULL);
        }
        free(peakList);
		free(mask);
    }
    return 0;
} 