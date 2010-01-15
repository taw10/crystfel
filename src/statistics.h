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

double stat_r2(double *obs, double *calc, unsigned int *c, int size,
               double *scalep);

#endif	/* STATISTICS_H */
