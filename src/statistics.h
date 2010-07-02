/*
 * statistics.h
 *
 * Structure-factor statistics
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef STATISTICS_H
#define STATISTICS_H

extern double stat_scale_intensity(const double *ref1, const unsigned int *c1,
                                   const double *ref2, const unsigned int *c2);

extern double stat_r2(const double *ref1, const unsigned int *c1,
                      const double *ref2, const unsigned int *c2,
                      double *scalep);

extern double stat_rmerge(const double *ref1, const unsigned int *c1,
                          const double *ref2, const unsigned int *c2,
                          double *scalep);

#endif	/* STATISTICS_H */
