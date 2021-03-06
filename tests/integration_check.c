/*
 * integration_check.c
 *
 * Check reflection integration
 *
 * Copyright © 2013-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2013-2020 Thomas White <taw@physics.org>
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

#include <image.h>
#include <utils.h>
#include <cell.h>
#include <cell-utils.h>
#include <integration.h>

#include "histogram.h"

int main(int argc, char *argv[])
{
	struct image image;
	int fs, ss;
	FILE *fh;
	unsigned long int seed;
	int fail = 0;
	const int w = 128;
	const int h = 128;
	RefList *list;
	Reflection *refl;
	UnitCell *cell;
	struct intcontext *ic;
	const int ir_inn = 2;
	const int ir_mid = 4;
	const int ir_out = 6;
	int i;
	Histogram *hi;
	double esd_sum = 0.0;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	fh = fopen("/dev/urandom", "r");
	if ( fread(&seed, sizeof(seed), 1, fh) == 1 ) {
		gsl_rng_set(rng, seed);
	} else {
		ERROR("Failed to seed RNG\n");
	}
	fclose(fh);

	image.lambda = ph_eV_to_lambda(9000.0);

	image.detgeom = calloc(1, sizeof(struct detgeom));
	image.detgeom->n_panels = 1;
	image.detgeom->panels = calloc(1, sizeof(struct detgeom_panel));

	image.detgeom->panels[0].w = w;
	image.detgeom->panels[0].h = h;
	image.detgeom->panels[0].fsx = 1.0;
	image.detgeom->panels[0].fsy = 0.0;
	image.detgeom->panels[0].ssx = 0.0;
	image.detgeom->panels[0].ssy = 1.0;
	image.detgeom->panels[0].cnx = -w/2;
	image.detgeom->panels[0].cny = -h/2;
	image.detgeom->panels[0].cnz = 60.0e-3 / 100e-6;
	image.detgeom->panels[0].pixel_pitch = 100e-6;  /* 10 px per mm */
	image.detgeom->panels[0].adu_per_photon = 10;
	image.detgeom->panels[0].max_adu = +INFINITY;  /* No cutoff */

	image.dp = malloc(sizeof(float *));
	image.dp[0] = malloc(w*h*sizeof(float));
	memset(image.dp[0], 0, w*h*sizeof(float));
	image.bad = malloc(sizeof(int *));
	image.bad[0] = malloc(w*h*sizeof(int));
	memset(image.bad[0], 0, w*h*sizeof(int));
	image.sat = NULL;

	image.n_crystals = 0;
	image.crystals = NULL;

	hi = histogram_init();

	for ( i=0; i<300; i++ ) {

		for ( fs=0; fs<w; fs++ ) {
		for ( ss=0; ss<h; ss++ ) {
			image.dp[0][fs+w*ss] = 10.0*poisson_noise(rng, 40);
			if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) > 2 ) continue;
			//image.dp[0][fs+w*ss] += 10.0*poisson_noise(10);
		}
		}

		list = reflist_new();
		refl = add_refl(list, 0, 0, 0);
		set_detector_pos(refl, 64, 64);
		set_panel_number(refl, 0);
		cell = cell_new();
		cell_set_lattice_type(cell, L_CUBIC);
		cell_set_centering(cell, 'P');
		cell_set_parameters(cell, 800.0e-10, 800.0e-10, 800.0e-10,
		                    deg2rad(90.0), deg2rad(90.0), deg2rad(90.0));
		cell = cell_rotate(cell, random_quaternion(rng));

		ic = intcontext_new(&image, cell, INTEGRATION_RINGS,
		                    ir_inn, ir_mid, ir_out, NULL);
		if ( ic == NULL ) {
			ERROR("Failed to initialise integration.\n");
			return 1;
		}

		integrate_rings_once(refl, ic, 0);

		cell_free(cell);

		histogram_add_value(hi, get_intensity(refl));
		esd_sum += get_esd_intensity(refl);

	}
	printf("Mean calculated sigma(I) = %.2f\n", esd_sum / 300);

	histogram_show(hi);

	histogram_free(hi);
	detgeom_free(image.detgeom);
	free(image.dp[0]);
	free(image.dp);

	if ( fail ) return 1;

	return 0;
}
