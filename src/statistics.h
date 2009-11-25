/*
 * statistics.h
 *
 * Structure-factor statistics
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * integr_sim - Test relrod integration
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
