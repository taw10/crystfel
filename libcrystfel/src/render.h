/*
 * render.h
 *
 * Render a high dynamic range buffer in some sensible way
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2012 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef RENDER_H
#define RENDER_H


enum {
	SCALE_COLOUR,
	SCALE_MONO,
	SCALE_INVMONO,
	SCALE_RATIO,
	SCALE_GEOPTIMISER
};

#ifdef __cplusplus
extern "C" {
#endif

/* Colour scale lookup */
extern void render_scale(double val, double max, int scale,
                         double *rp, double *gp, double *bp);


#ifdef __cplusplus
}
#endif

#endif	/* RENDER_H */
