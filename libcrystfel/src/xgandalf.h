/*
 * xgandalf.h
 *
 *  Created on: 08.08.2017
 *      Author: gevorkov
 */

#ifndef LIBCRYSTFEL_SRC_XGANDALF_H_
#define LIBCRYSTFEL_SRC_XGANDALF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct xgandalf_options {
	unsigned int sampling_pitch;
	unsigned int grad_desc_iteration_selector;
	float tolerance;
	unsigned int no_deviation_from_provided_cell;
	float minLatticeVectorLength_A;
	float maxLatticeVectorLength_A;
};

#include <stddef.h>
#include "index.h"

int run_xgandalf(struct image *image, void *ipriv);

void *xgandalf_prepare(IndexingMethod *indm, UnitCell *cell,
		struct xgandalf_options *xgandalf_opts);

void xgandalf_cleanup(void *pp);
const char *xgandalf_probe(UnitCell *cell);


#endif /* LIBCRYSTFEL_SRC_XGANDALF_H_ */
