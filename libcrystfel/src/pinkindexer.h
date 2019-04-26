/*
 * pinkIndexer.h
 *
 *  Created on: Nov 27, 2017
 *      Author: gevorkov
 */

#ifndef LIBCRYSTFEL_SRC_PINKINDEXER_H_
#define LIBCRYSTFEL_SRC_PINKINDEXER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct pinkIndexer_options {
	unsigned int considered_peaks_count;
	unsigned int angle_resolution;
	unsigned int refinement_type;
	float maxResolutionForIndexing_1_per_A;
	float tolerance;
	int multi;
	int thread_count;
	int min_peaks;

	int no_check_indexed;

	float beamEnergy; //in eV
	float beamBandwidth; //(delta lambda)/lambda
	float detectorDistance; //in m

	float reflectionRadius; //in 1/A
};

#include <stddef.h>
#include "index.h"

extern int run_pinkIndexer(struct image *image, void *ipriv);

extern void *pinkIndexer_prepare(IndexingMethod *indm, UnitCell *cell,
        struct pinkIndexer_options *pinkIndexer_opts);

extern void pinkIndexer_cleanup(void *pp);

extern const char *pinkIndexer_probe(UnitCell *cell);

#endif /* LIBCRYSTFEL_SRC_PINKINDEXER_H_ */
