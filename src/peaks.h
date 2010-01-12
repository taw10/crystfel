/*
 * peaks.h
 *
 * Peak search and other image analysis
 *
 * (c) 2006-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef PEAKS_H
#define PEAKS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


extern int image_fom(struct image *image);
extern void search_peaks(struct image *image, int dump_peaks);


#endif	/* PEAKS_H */
