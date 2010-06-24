/*
 * povray.c
 *
 * Invoke POV-ray
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef POVRAY_H
#define POVRAY_H

extern int povray_render_animation(UnitCell *cell, double *ref,
                                   unsigned int *c, unsigned int nproc);

#endif	/* POVRAY_H */
