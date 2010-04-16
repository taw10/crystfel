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

extern void search_peaks(struct image *image);
extern void dump_peaks(struct image *image, pthread_mutex_t *mutex);
extern void output_intensities(struct image *image, UnitCell *cell,
                               pthread_mutex_t *mutex);

#endif	/* PEAKS_H */
