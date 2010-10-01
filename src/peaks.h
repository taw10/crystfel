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

extern void search_peaks(struct image *image, float threshold);
extern void dump_peaks(struct image *image, pthread_mutex_t *mutex);
extern void output_intensities(struct image *image, UnitCell *cell,
                               pthread_mutex_t *mutex, int polar, int sa,
                               int use_closer, FILE *ofh, int circular_domain,
                               double domain_r);
extern int peak_sanity_check(struct image *image, UnitCell *cell,
                             int circular_domain, double domain_r);
extern int find_projected_peaks(struct image *image, UnitCell *cell,
                                int circular_domain, double domain_r);
extern int integrate_peak(struct image *image, int xp, int yp,
                          float *xc, float *yc, float *intensity,
                          int do_polar, int do_sa);
#endif	/* PEAKS_H */
