//#################################################################################################
//# This file is part of RobustGaussianFittingLibrary, a free library WITHOUT ANY WARRANTY        #
//# Copyright: 2017-2020 LaTrobe University Melbourne, 2019-2020 Deutsches Elektronen-Synchrotron #
//#################################################################################################

#ifndef RGFLIB_H
#define RGFLIB_H

#define NEGATIVE_MAX -(1e+10)

void islandRemoval(unsigned char* inMask, unsigned char* outMask,
 	 			   unsigned int X, unsigned int Y,
				   unsigned int islandSizeThreshold);

float MSSE(float *error, unsigned int vecLen, float MSSE_LAMBDA, 
			unsigned int k, float minimumResidual);
float MSSEWeighted(float* error, float* weights, unsigned int vecLen,
                   float MSSE_LAMBDA, unsigned int k, float minimumResidual);
	
void fitValue(float* inVec,
			  float* inWeights,
			  float* modelParams,
			  float theta,
			  unsigned int inN,
              float topkPerc,
			  float botkPerc,
              float MSSE_LAMBDA,
			  unsigned char optIters,
              float minimumResidual,
			  unsigned int downSampledSize);

void fitValue2Skewed(float *vec, float *weights, 
					float *modelParams, float theta, unsigned int N,
					float topkPerc, float botkPerc,
					float MSSE_LAMBDA, unsigned char optIters, float minimumResidual);

void medianOfFits(float *vec, float *weights, 
  				  float *modelParams, float theta, unsigned int N,
				  float topkMin, float topkMax, unsigned int numSamples, float samplePerc,
				  float MSSE_LAMBDA, unsigned char optIters, float minimumResidual);

void RobustAlgebraicLineFitting(float* x, float* y, float* mP,
							unsigned int N, float topkPerc, float botkPerc, float MSSE_LAMBDA);

void RobustAlgebraicLineFittingTensor(float *inTensorX, float *inTensorY,
									float *modelParamsMap, unsigned int N,
									unsigned int X, unsigned int Y,
									float topkPerc, float botkPerc,
									float MSSE_LAMBDA);

void RobustAlgebraicPlaneFitting(float* x, float* y, float* z, float* mP, float* mP_Init,
							unsigned int N, float topkPerc, float botkPerc,
							float MSSE_LAMBDA, unsigned char stretch2CornersOpt, 
							float minimumResidual, unsigned char optIters);

void fitValueTensor(float *inTensor, unsigned char* inWeights,
                float *modelParamsMap, unsigned int N, unsigned int X, unsigned int Y,
                float topkPerc, float botkPerc, float MSSE_LAMBDA,
                unsigned char optIters, float minimumResidual,
				unsigned int downSampledSize);

void RSGImage(float* inImage, unsigned char* inMask, float *modelParamsMap,
				unsigned int winX, unsigned int winY,
				unsigned int X, unsigned int Y,
				float topkPerc, float botkPerc,
				float MSSE_LAMBDA, unsigned char stretch2CornersOpt,
				unsigned char numModelParams, unsigned char optIters,
				float minimumResidual);

void RSGImage_by_Image_Tensor(float* inImage_Tensor, unsigned char* inMask_Tensor,
						float *model_mean, float *model_std,
						unsigned int winX, unsigned int winY,
						unsigned int N, unsigned int X, unsigned int Y,
						float topkPerc, float botkPerc,
						float MSSE_LAMBDA, unsigned char stretch2CornersOpt,
						unsigned char numModelParams, unsigned char optIters,
						float minimumResidual);

void fitBackgroundRadially(float* inImage, unsigned char* inMask, 
						   float* modelParamsMap, float* vecMP,
 						   unsigned int minRes, 
						   unsigned int maxRes, 
						   unsigned int shellWidth,
						   unsigned char includeCenter, 
						   unsigned int finiteSampleBias,
						   unsigned int X, unsigned int Y,
						   float topkPerc, float botkPerc, 
						   float MSSE_LAMBDA, 
						   unsigned char optIters,
						   float minimumResidual);
#endif
