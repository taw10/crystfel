/*
 * peaks.h
 *
 * Peak search and other image analysis
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifndef PEAKS_H
#define PEAKS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <pthread.h>

#include "reflist.h"

extern void search_peaks(struct image *image, float threshold,
                         float min_gradient, float min_snr);

extern void integrate_reflections(struct image *image,
                                  int polar, int use_closer, int bgsub);

extern int peak_sanity_check(struct image * image);

/* Exported so it can be poked by integration_check */
extern int integrate_peak(struct image *image, int cfs, int css,
                          double *pfs, double *pss, double *intensity,
                          double *pbg, double *pmax, double *sigma,
                          int do_polar, int centroid, int bgsub);

#endif	/* PEAKS_H */
