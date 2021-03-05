/*
 * colscale.h
 *
 * Colour scales
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2009-2020 Thomas White <taw@physics.org>
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

#ifndef COLSCALE_H
#define COLSCALE_H

/**
 * \file colscale.h
 * Colour scales for rendering
 */

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
extern void colscale_lookup(double val, double max, int scale,
                            double *rp, double *gp, double *bp);


#ifdef __cplusplus
}
#endif

#endif	/* COLSCALE_H */
