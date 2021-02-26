//#################################################################################################
//# This file is part of RobustPeakFinder, a free library WITHOUT ANY WARRANTY                    #
//# Copyright: 2017-2020 LaTrobe University Melbourne, 2019-2020 Deutsches Elektronen-Synchrotron #
//#################################################################################################

#ifndef RPFSOURCE_H
#define RPFSOURCE_H

int rpfMain(float *inData, 
			unsigned char use_Mask,
            unsigned char *inMask, 
			unsigned char use_peakMask,
            unsigned char *peakMask,
            unsigned char minBackMeanHasAMap,
            float *minBackMeanMap, 
            unsigned char maxBackMeanHasAMap,
            float *maxBackMeanMap, 
			unsigned char returnPeakMap,
            float *peakMap,
            float *peakList,
            float singlePhotonADU,
            int MAXIMUM_NUMBER_OF_PEAKS,
            float bckSNR, 
            float pixPAPR,
            int XPIX, 
            int YPIX, 
            int PTCHSZ,    
            int PEAK_MIN_PIX, 
            int PEAK_MAX_PIX,
            int n_optIters, 
            int finiteSampleBias,
            int downSampledSize, 
            float MSSE_LAMBDA, 
            float searchSNR,
            float highPoissonTh, 
            float lowPoissonTh);
#endif