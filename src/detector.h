/*
 * detector.h
 *
 * Detector properties
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
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
#define UPPER_CX (491.9)
#define UPPER_CY (440.7)
#define LOWER_CX (492.0)
#define LOWER_CY (779.7)

extern void record_image(struct image *image, int do_water, int do_poisson,
                         int do_bloom);

#endif	/* DETECTOR_H */
