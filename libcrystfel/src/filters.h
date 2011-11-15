/*
 * peaks.h
 *
 * Image filtering
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef FILTERS_H
#define FILTERS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


extern void filter_cm(struct image *image);
extern void filter_noise(struct image *image, float *old);


#endif	/* FILTERS_H */
