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
#include <stdlib.h>

#include "detgeom.h"
#include "utils.h"


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


void detgeom_free(struct detgeom *detgeom)
{
	int i;

	for ( i=0; i<detgeom->n_panels; i++ ) {
		free(detgeom->panels[i].name);
	}

	free(detgeom->panels);
	free(detgeom);
}


static double panel_max_res(struct detgeom_panel *p,
                            double wavelength)
{
	double r[3];
	double max_res = 0.0;

	detgeom_transform_coords(p, 0, 0, wavelength, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	detgeom_transform_coords(p, 0, p->h, wavelength, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	detgeom_transform_coords(p, p->w, 0, wavelength, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	detgeom_transform_coords(p, p->w, p->h, wavelength, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	return max_res;
}


double detgeom_max_resolution(struct detgeom *detgeom,
                              double wavelength)
{
	int i;
	double max_res = 0.0;

	for ( i=0; i<detgeom->n_panels; i++ ) {

		double panel_maxres;

		panel_maxres = panel_max_res(&detgeom->panels[i],
		                             wavelength);
		if ( panel_maxres > max_res ) {
			max_res = panel_maxres;
		}
	}

	return max_res;
}


void show_panel(struct detgeom_panel *p)
{
	STATUS("Panel '%s':\n", p->name);
	STATUS("  Size %i x %i px\n", p->w, p->h);
	STATUS("  Transformation [cnx] + [%6.2f %6.2f] [fs] = [x]\n",
	       p->fsx, p->ssx);
	STATUS("                 [cny] + [%6.2f %6.2f] [ss] = [y]\n",
	       p->fsy, p->ssy);
	STATUS("                 [cnz] + [%6.2f %6.2f]      = [z]\n",
	       p->fsz, p->ssz);
	STATUS("  corner x,y,z = %f, %f, %f px\n",
	       p->cnx, p->cny, p->cnz);
	STATUS("               = %f, %f, %f mm\n",
	       p->cnx*p->pixel_pitch*1e3,
	       p->cny*p->pixel_pitch*1e3,
	       p->cnz*p->pixel_pitch*1e3);
	STATUS("  %f adu/photon, max %f adu\n",
	       p->adu_per_photon, p->max_adu);
}
