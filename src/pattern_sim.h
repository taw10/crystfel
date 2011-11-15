/*
 * pattern_sim.h
 *
 * Simulate diffraction patterns from small crystals
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PATTERN_SIM_H
#define PATTERN_SIM_H


/* Maxmimum index to hold values up to (can be increased if necessary)
 * WARNING: Altering this value constitutes an ABI change, and means you must
 * update data/diffraction.cl then recompile and reinstall everything. */
#define INDMAX 200

/* Array size */
#define IDIM (INDMAX*2 +1)
#define LIST_SIZE (IDIM*IDIM*IDIM)

/* Create functions for storing reflection intensities indexed as h,k,l */
#define LABEL(x) x##_intensity
#define TYPE double
#include "list_tmp.h"

/* CAs above, but for phase values */
#define LABEL(x) x##_phase
#define TYPE double
#include "list_tmp.h"

/* As above, but for (unsigned) integer counts */
#define LABEL(x) x##_count
#define TYPE unsigned int
#include "list_tmp.h"

/* As above, but for simple flags */
#define LABEL(x) x##_flag
#define TYPE unsigned char
#include "list_tmp.h"


#endif	/* PATTERN_SIM_H */
