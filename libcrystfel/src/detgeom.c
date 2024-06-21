/*
 * detgeom.c
 *
 * Utility functions for detgeom structure
 *
 * Copyright © 2019-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "detgeom.h"
#include "utils.h"


/**
 * \file detgeom.h
 */


void detgeom_transform_coords(struct detgeom_panel *p,
                              double fs, double ss,
                              double wavelength,
                              double dx, double dy,
                              double *r)
{
	double xs, ys, zs;
	double fac;

	/* Calculate 3D position of given position, in pixels */
	xs = p->cnx + fs*p->fsx + ss*p->ssx + dx*p->pixel_pitch;
	ys = p->cny + fs*p->fsy + ss*p->ssy + dy*p->pixel_pitch;
	zs = p->cnz + fs*p->fsz + ss*p->ssz;

	fac = wavelength * sqrt(xs*xs + ys*ys + zs*zs);

	r[0] = xs / fac;
	r[1] = ys / fac;
	r[2] = zs / fac - 1.0/wavelength;
}


static void free_group(struct detgeom_panel_group *g)
{
	int i;

	if ( g == NULL ) return;

	for ( i=0; i<g->n_children; i++ ) {
		free_group(g->children[i]);
	}

	cffree(g->name);
	cffree(g->children);
	cffree(g);
}


void detgeom_free(struct detgeom *detgeom)
{
	int i;

	if ( detgeom == NULL ) return;

	for ( i=0; i<detgeom->n_panels; i++ ) {
		cffree(detgeom->panels[i].name);
	}

	free_group(detgeom->top_group);
	cffree(detgeom->panels);
	cffree(detgeom);
}


static double panel_max_res(struct detgeom_panel *p,
                            double wavelength)
{
	double r[3];
	double max_res = 0.0;

	detgeom_transform_coords(p, 0, 0, wavelength, 0.0, 0.0, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	detgeom_transform_coords(p, 0, p->h, wavelength, 0.0, 0.0, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	detgeom_transform_coords(p, p->w, 0, wavelength, 0.0, 0.0, r);
	max_res = biggest(max_res, modulus(r[0], r[1], r[2]));

	detgeom_transform_coords(p, p->w, p->h, wavelength, 0.0, 0.0, r);
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


static double clen(struct detgeom *dg, int i)
{
	return dg->panels[i].cnz * dg->panels[i].pixel_pitch;
}


static double far_corner_clen(struct detgeom *dg, int i)
{
	struct detgeom_panel *p = &dg->panels[i];
	return (p->cnz + p->fsz*p->w + p->ssz*p->h) * p->pixel_pitch;
}


/* Return mean clen in m, or NAN in the following circumstances:
 * 1. If the individual panel distances vary by more than 10% of the average
 * 2. If the tilt of the panel creates a distance variation of more than 10%
 *    of the corner value over the extent of the panel */
double detgeom_mean_camera_length(struct detgeom *dg)
{
	int i;
	double mean;
	double total = 0.0;
	for ( i=0; i<dg->n_panels; i++ ) {
		total += clen(dg, i);
	}
	mean = total / dg->n_panels;

	for ( i=0; i<dg->n_panels; i++ ) {
		if ( !within_tolerance(clen(dg, i), mean, 10.0) ) return NAN;
		if ( !within_tolerance(clen(dg, i), far_corner_clen(dg, i), 10.0) ) return NAN;
	}

	return mean;
}


struct detgeom_panel *detgeom_find_panel(struct detgeom *dg, const char *name)
{
	int i;
	for ( i=0; i<dg->n_panels; i++ ) {
		if ( strcmp(dg->panels[i].name, name) == 0 ) {
			return &dg->panels[i];
		}
	}
	return NULL;
}


static void detgeom_show_group(const struct detgeom_panel_group *group, int level)
{
	int i;

	for ( i=0; i<level; i++ ) STATUS("  ");

	if ( group == NULL ) {
		STATUS("!!!\n");
		return;
	}

	STATUS("%s (serial %i)\n", group->name, group->serial);

	for ( i=0; i<group->n_children; i++ ) {
		detgeom_show_group(group->children[i], level+1);
	}
}


void detgeom_show_hierarchy(const struct detgeom *dg)
{
	detgeom_show_group(dg->top_group, 0);
}


void detgeom_translate_detector_m(struct detgeom *dg, double x, double y, double z)
{
	int i;
	for ( i=0; i<dg->n_panels; i++ ) {
		struct detgeom_panel *p = &dg->panels[i];
		p->cnx += x / p->pixel_pitch;
		p->cny += y / p->pixel_pitch;
		p->cnz += z / p->pixel_pitch;
	}
}


gsl_matrix **make_panel_minvs(struct detgeom *dg)
{
	int i;
	gsl_matrix **Minvs;

	Minvs = cfmalloc(dg->n_panels * sizeof(gsl_matrix *));
	if ( Minvs == NULL ) return NULL;

	for ( i=0; i<dg->n_panels; i++ ) {

		struct detgeom_panel *p = &dg->panels[i];
		gsl_matrix *M = gsl_matrix_calloc(3, 3);

		gsl_matrix_set(M, 0, 0, p->pixel_pitch*p->cnx);
		gsl_matrix_set(M, 0, 1, p->pixel_pitch*p->fsx);
		gsl_matrix_set(M, 0, 2, p->pixel_pitch*p->ssx);
		gsl_matrix_set(M, 1, 0, p->pixel_pitch*p->cny);
		gsl_matrix_set(M, 1, 1, p->pixel_pitch*p->fsy);
		gsl_matrix_set(M, 1, 2, p->pixel_pitch*p->ssy);
		gsl_matrix_set(M, 2, 0, p->pixel_pitch*p->cnz);
		gsl_matrix_set(M, 2, 1, p->pixel_pitch*p->fsz);
		gsl_matrix_set(M, 2, 2, p->pixel_pitch*p->ssz);

		Minvs[i] = matrix_invert(M);
		gsl_matrix_free(M);
		if ( Minvs[i] == NULL ) {
			ERROR("Failed to calculate inverse panel matrix for %s\n",
			      p->name);
			return NULL;
		}

	}

	return Minvs;
}


static void add_point(const struct detgeom_panel *p,
                      int fs, int ss,
                      double *tx, double *ty, double *tz)
{
	*tx += (p->cnx + fs*p->fsx + ss*p->ssx) * p->pixel_pitch;
	*ty += (p->cny + fs*p->fsy + ss*p->ssy) * p->pixel_pitch;
	*tz += (p->cnz + fs*p->fsz + ss*p->ssz) * p->pixel_pitch;
}


int detgeom_group_center(const struct detgeom_panel_group *grp,
                         double *x, double *y, double *z)
{
	if ( grp->n_children == 0 ) {

		const struct detgeom_panel *p = grp->panel;
		if ( p == NULL ) return 1;

		double tx = 0.0;
		double ty = 0.0;
		double tz = 0.0;

		add_point(p, 0, 0, &tx, &ty, &tz);
		add_point(p, p->w, 0, &tx, &ty, &tz);
		add_point(p, 0, p->h, &tx, &ty, &tz);
		add_point(p, p->w, p->h, &tx, &ty, &tz);

		*x = tx / 4.0;
		*y = ty / 4.0;
		*z = tz / 4.0;

		return 0;

	} else {

		int i;
		double tx = 0.0;
		double ty = 0.0;
		double tz = 0.0;

		for ( i=0; i<grp->n_children; i++ ) {
			double gcx, gcy, gcz;
			detgeom_group_center(grp->children[i], &gcx, &gcy, &gcz);
			tx += gcx;
			ty += gcy;
			tz += gcz;
		}

		*x = tx / grp->n_children;
		*y = ty / grp->n_children;
		*z = tz / grp->n_children;

		return 0;

	}
}
