/*
 * detector.h
 *
 * Detector properties
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef DETECTOR_H
#define DETECTOR_H

#include "image.h"

extern void record_image(struct image *image, int do_water);

#endif	/* DETECTOR_H */
