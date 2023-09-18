/*
 * gradient_check_utils.c
 *
 * Check gradients for prediction refinement (common component)
 *
 * Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2023 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <getopt.h>

#include <image.h>
#include <cell.h>
#include <cell-utils.h>
#include <geometry.h>
#include <reflist.h>


double **make_dev_list(struct reflpeak *rps, int n_refls, struct detgeom *det)
{
	int i;
	double **vals;

	vals = malloc(3*sizeof(double *));
	vals[0] = malloc(n_refls*sizeof(double));
	vals[1] = malloc(n_refls*sizeof(double));
	vals[2] = malloc(n_refls*sizeof(double));

	for ( i=0; i<n_refls; i++ ) {
		vals[0][i] = get_exerr(rps[i].refl);
		vals[1][i] = fs_dev(&rps[i], det);
		vals[2][i] = ss_dev(&rps[i], det);
	}

	return vals;
}


static UnitCell *random_rotated_cell(gsl_rng *rng)
{
	UnitCell *cell;
	struct quaternion orientation;

	cell = cell_new_from_parameters(10.0e-9, 10.0e-9, 10.0e-9,
	                                deg2rad(90.0),
	                                deg2rad(90.0),
	                                deg2rad(90.0));
	orientation = random_quaternion(rng);
	return cell_rotate(cell, orientation);
}


static void rot(double *x, double *y, double ang)
{
	double nx, ny;
	nx = (*x)*cos(ang) - (*y)*sin(ang);
	ny = (*x)*sin(ang) + (*y)*cos(ang);
	*x = nx;  *y = ny;
}


struct reflpeak *make_test_image(int *pn_refls, struct image *image)
{
	Crystal *cr;
	gsl_rng *rng;
	RefList *refls;
	Reflection *refl;
	RefListIterator *iter;
	int n_refls;
	int i;
	struct reflpeak *rps;

	image->detgeom = malloc(sizeof(struct detgeom));
	image->detgeom->n_panels = 1;
	image->detgeom->panels = malloc(sizeof(struct detgeom_panel));
	image->detgeom->panels[0].name = "panel";
	image->detgeom->panels[0].adu_per_photon = 1.0;
	image->detgeom->panels[0].max_adu = INFINITY;
	image->detgeom->panels[0].fsx = 1.0;
	image->detgeom->panels[0].fsy = 0.0;
	image->detgeom->panels[0].fsz = 0.0;
	image->detgeom->panels[0].ssx = 0.0;
	image->detgeom->panels[0].ssy = 1.0;
	image->detgeom->panels[0].ssz = 0.0;
	rot(&image->detgeom->panels[0].fsx,
	    &image->detgeom->panels[0].fsy,
	    deg2rad(30));
	rot(&image->detgeom->panels[0].ssx,
	    &image->detgeom->panels[0].ssy,
	    deg2rad(30));
	rot(&image->detgeom->panels[0].fsx,
	    &image->detgeom->panels[0].fsz,
	    deg2rad(15));
	rot(&image->detgeom->panels[0].ssx,
	    &image->detgeom->panels[0].ssz,
	    deg2rad(15));
	image->detgeom->panels[0].cnx = -500.0;
	image->detgeom->panels[0].cny = -500.0;
	image->detgeom->panels[0].cnz = 1000.0; /* pixels */
	image->detgeom->panels[0].w = 1000;
	image->detgeom->panels[0].h = 1000;
	image->detgeom->panels[0].pixel_pitch = 75e-6;

	image->lambda = ph_en_to_lambda(eV_to_J(8000.0));
	image->div = 1e-3;
	image->bw = 0.00001;
	image->filename = malloc(256);
	image->spectrum = spectrum_generate_gaussian(image->lambda, image->bw);

	image->crystals = NULL;
	image->n_crystals = 0;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return NULL;
	}
	crystal_set_mosaicity(cr, 0.0);
	crystal_set_profile_radius(cr, 0.005e9);
	crystal_set_image(cr, image);
	crystal_set_cell(cr, random_rotated_cell(rng));

	refls = predict_to_res(cr, detgeom_max_resolution(image->detgeom, image->lambda));
	crystal_set_reflections(cr, refls);
	n_refls = num_reflections(refls);

	/* Associate a peak with every reflection */
	rps = malloc(n_refls*sizeof(struct reflpeak));
	i = 0;
	for ( refl = first_refl(refls, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double fs, ss;
		int pn;
		struct imagefeature *pk;

		get_detector_pos(refl, &fs, &ss);
		pn = get_panel_number(refl);

		pk = malloc(sizeof(struct imagefeature));

		pk->fs = fs + gsl_ran_gaussian(rng, 5.0);
		pk->ss = ss + gsl_ran_gaussian(rng, 5.0);
		pk->pn = pn;
		pk->intensity = 1.0;

		rps[i].peak = pk;
		rps[i].refl = refl;
		rps[i].Ih = 1.0;
		i++;
	}

	image_add_crystal(image, cr);

	gsl_rng_free(rng);

	*pn_refls = n_refls;
	return rps;
}
