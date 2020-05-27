/*
 * detgeom.c
 *
 * Utility functions for detgeom structure
 *
 * Copyright © 2019-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

#include <math.h>

#include "detgeom.h"


/**
 * \file detgeom.h
 */


void detgeom_transform_coords(struct detgeom_panel *p,
                              double fs, double ss,
                              double wavelength,
                              double *r)
{
	double ctt, twotheta;
	double xs, ys, zs;
	double az;

	/* Calculate 3D position of given position, in m */
	xs = (p->cnx + fs*p->fsx + ss*p->ssx) * p->pixel_pitch;
	ys = (p->cny + fs*p->fsy + ss*p->ssy) * p->pixel_pitch;
	zs = (p->cnz + fs*p->fsz + ss*p->ssz) * p->pixel_pitch;

	ctt = zs/sqrt(xs*xs + ys*ys + zs*zs);
	twotheta = acos(ctt);
	az = atan2(ys, xs);

	r[0] = sin(twotheta)*cos(az) / wavelength;
	r[1] = sin(twotheta)*sin(az) / wavelength;
	r[2] = (ctt - 1.0) / wavelength;
}
