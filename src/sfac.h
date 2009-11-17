/*
 * sfac.h
 *
 * Scattering factors
 *
 * (c) 2007-2009 Thomas White <thomas.white@desy.de>
 *
 * pattern_sim - Simulate diffraction patterns from small crystals
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SFAC_H
#define SFAC_H

#include <complex.h>

extern double complex get_sfac(const char *n, double s, double en);

#endif	/* SFAC_H */
