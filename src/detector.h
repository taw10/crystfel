/*
 * detector.h
 *
 * Detector properties
 *
 * (c) 2007-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DETECTOR_H
#define DETECTOR_H

#include "image.h"

/* Position of central beam for upper and lower CCDs */
#define UPPER_CX (492.8)
#define UPPER_CY (437.6)
#define LOWER_CX (494.0)
#define LOWER_CY (772.1)

extern void record_image(struct image *image, int do_water, int do_poisson,
                         int do_bloom);

#endif	/* DETECTOR_H */
