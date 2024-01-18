/*
 * gradient_panel_x.c
 *
 * Check gradients for prediction refinement
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
#include <getopt.h>

#include <image.h>
#include <geometry.h>
#include <predict-refine.h>

#include "gradient_check_utils.h"


int main(int argc, char *argv[])
{
	struct image image;
	struct reflpeak *rps;
	int n_refls;
	double **before;
	double **after;
	int i;
	int n_wrong_r = 0;
	int n_wrong_fs = 0;
	int n_wrong_ss = 0;
	int n_wrong_obsr = 0;
	int fail = 0;
	double step;
	gsl_matrix **panel_matrices;
	int didsomething = 0;
	const double cx = 0.03;  /* Detector is a 7.5 cm side length square */
	const double cy = 0.02;
	const double cz = 0.01;

	rps = make_test_image(&n_refls, &image);
	panel_matrices = make_panel_minvs(image.detgeom);

	before = make_dev_list(rps, n_refls, image.detgeom);

	#ifdef TRANSLATE_PANEL
	struct detgeom_panel *p = &image.detgeom->panels[0];
	step = 0.01e-3;  /* metres */
	image.detgeom->panels[0].THING_TO_MOVE += step/p->pixel_pitch;
	didsomething = 1;
	#endif

	#ifdef ROTATE_PANEL_X
	struct detgeom_panel *p = &image.detgeom->panels[0];
	step = deg2rad(0.01);
	rotate2d(&p->cny, &p->cnz, cy/p->pixel_pitch, cz/p->pixel_pitch, step);
	rotate2d(&p->fsy, &p->fsz, 0, 0, step);
	rotate2d(&p->ssy, &p->ssz, 0, 0, step);
	didsomething = 1;
	#endif

	#ifdef ROTATE_PANEL_Y
	struct detgeom_panel *p = &image.detgeom->panels[0];
	step = deg2rad(0.01);
	rotate2d(&p->cnz, &p->cnx, cz/p->pixel_pitch, cx/p->pixel_pitch, step);
	rotate2d(&p->fsz, &p->fsx, 0, 0, step);
	rotate2d(&p->ssz, &p->ssx, 0, 0, step);
	didsomething = 1;
	#endif

	#ifdef ROTATE_PANEL_Z
	struct detgeom_panel *p = &image.detgeom->panels[0];
	step = deg2rad(0.01);
	rotate2d(&p->cnx, &p->cny, cx/p->pixel_pitch, cy/p->pixel_pitch, step);
	rotate2d(&p->fsx, &p->fsy, 0, 0, step);
	rotate2d(&p->ssx, &p->ssy, 0, 0, step);
	didsomething = 1;
	#endif

	#ifdef CHANGE_CELL
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;
	UnitCell *cell = crystal_get_cell(image.crystals[0].cr);
	step = 0.5e5;
	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);
	THING_TO_MOVE += step;
	cell_set_reciprocal(cell, asx, asy, asz,
	                          bsx, bsy, bsz,
	                          csx, csy, csz);
	didsomething = 1;
	#endif

	if ( !didsomething ) {
		fprintf(stderr, "Nothing changed.  Check the build system.\n");
		return 1;
	}

	update_predictions(image.crystals[0].refls, image.crystals[0].cr, &image);
	after = make_dev_list(rps, n_refls, image.detgeom);

	for ( i=0; i<n_refls; i++ ) {

		float calc[3];
		double obs[3];

		calc[0] = r_gradient(TEST_GPARAM, rps[i].refl,
		                     crystal_get_cell(image.crystals[0].cr),
		                     image.lambda);

		fs_ss_gradient(TEST_GPARAM, rps[i].refl,
		               crystal_get_cell(image.crystals[0].cr),
		               &image.detgeom->panels[rps[i].peak->pn],
		               panel_matrices[rps[i].peak->pn], cx, cy, cz,
		               &calc[1], &calc[2]);

		obs[0] = (after[0][i] - before[0][i]) / step;
		obs[1] = (after[1][i] - before[1][i]) / step;
		obs[2] = (after[2][i] - before[2][i]) / step;

		#ifdef TRANSLATE_PANEL
		if ( fabs(calc[0]) > 1e-12 ) n_wrong_r++;  /* Should be zero */
		if ( fabs(obs[0]) > 1e-12 ) n_wrong_obsr++;  /* Should also be zero */
		if ( fabs(obs[1] - calc[1]) > 10.0 ) n_wrong_fs++;
		if ( fabs(obs[2] - calc[2]) > 10.0 ) n_wrong_ss++;
		#endif

		#if defined(ROTATE_PANEL_X) || defined(ROTATE_PANEL_Y) || defined(ROTATE_PANEL_Z)
		if ( fabs(calc[0]) > 1e-12 ) n_wrong_r++;  /* Should be zero */
		if ( fabs(obs[0]) > 1e-12 ) n_wrong_obsr++;  /* Should also be zero */
		if ( fabs(obs[1] - calc[1]) > 1.0 ) n_wrong_fs++;  /* Units are pixels/rad */
		if ( fabs(obs[2] - calc[2]) > 1.0 ) n_wrong_ss++;  /* (numbers are big) */
		#endif

		#ifdef CHANGE_CELL
		if ( fabs(obs[0] - calc[0]) > 1e-2 ) n_wrong_r++;
		if ( fabs(obs[1] - calc[1]) > 1e-8 ) n_wrong_fs++;
		if ( fabs(obs[2] - calc[2]) > 1e-8 ) n_wrong_ss++;
		#endif

	}

	if ( n_wrong_r > 0 ) {
		fprintf(stderr, "%i out of %i R gradients were wrong.\n",
		        n_wrong_r, n_refls);
		fail = 1;
	}

	if ( n_wrong_fs > 0 ) {
		fprintf(stderr, "%i out of %i fs gradients were wrong.\n",
		        n_wrong_fs, n_refls);
		fail = 1;
	}

	if ( n_wrong_ss > 0 ) {
		fprintf(stderr, "%i out of %i ss gradients were wrong.\n",
		        n_wrong_ss, n_refls);
		fail = 1;
	}

	if ( n_wrong_obsr > 0 ) {
		fprintf(stderr, "%i out of %i observed R gradients were not zero as expected\n",
		        n_wrong_obsr, n_refls);
		fail = 1;
	}

	return fail;
}
