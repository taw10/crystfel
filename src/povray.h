/*
 * povray.h
 *
 * Invoke POV-ray
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef POVRAY_H
#define POVRAY_H

#include "reflist.h"
#include "cell.h"

extern int povray_render_animation(UnitCell *cell, RefList *list,
                                   unsigned int nproc, const char *sym,
                                   int wght, double boost);

#endif	/* POVRAY_H */
