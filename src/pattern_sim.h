/*
 * pattern_sim.h
 *
 * Simulate diffraction patterns from small crystals
 *
 * Copyright © 2012 Thomas White <taw@physics.org>
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

#ifndef PATTERN_SIM_H
#define PATTERN_SIM_H


/* Maxmimum index to hold values up to (can be increased if necessary)
 * WARNING: Altering this value constitutes an ABI change, and means you must
 * update data/diffraction.cl then recompile and reinstall everything. */
#define INDMAX 120

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
