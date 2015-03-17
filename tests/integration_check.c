/*
 * integration_check.c
 *
 * Check reflection integration
 *
 * Copyright © 2013 Thomas White <taw@physics.org>
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


#include <stdlib.h>
#include <stdio.h>

#include <image.h>
#include <utils.h>
#include <histogram.h>

#include "../libcrystfel/src/integration.c"


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
	struct intcontext ic;
	const int ir_inn = 2;
	const int ir_mid = 4;
	const int ir_out = 6;
	int i;
	Histogram *hi;
	double esd_sum = 0.0;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	fh = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fh);
	fclose(fh);
	gsl_rng_set(rng, seed);

	image.flags = NULL;
	image.beam = NULL;
	image.lambda = ph_eV_to_lambda(9000.0);

	image.det = calloc(1, sizeof(struct detector));
	image.det->n_panels = 1;
	image.det->panels = calloc(1, sizeof(struct panel));

	image.det->panels[0].min_fs = 0;
	image.det->panels[0].max_fs = w;
	image.det->panels[0].min_ss = 0;
	image.det->panels[0].max_ss = h;
	image.det->panels[0].w = w;
	image.det->panels[0].h = h;
	image.det->panels[0].fsx = 1.0;
	image.det->panels[0].fsy = 0.0;
	image.det->panels[0].ssx = 0.0;
	image.det->panels[0].ssy = 1.0;
	image.det->panels[0].xfs = 1.0;
	image.det->panels[0].yfs = 0.0;
	image.det->panels[0].xss = 0.0;
	image.det->panels[0].yss = 1.0;
	image.det->panels[0].cnx = -w/2;
	image.det->panels[0].cny = -h/2;
	image.det->panels[0].clen = 60.0e-3;
	image.det->panels[0].res = 100000;  /* 10 px per mm */
	image.det->panels[0].adu_per_eV = 10.0/9000.0; /* 10 adu/ph */
	image.det->panels[0].max_adu = +INFINITY;  /* No cutoff */

	image.width = w;
	image.height = h;
	image.dp = malloc(sizeof(float *));
	image.dp[0] = malloc(w*h*sizeof(float));
	memset(image.dp[0], 0, w*h*sizeof(float));
	image.bad = malloc(sizeof(int *));
	image.bad[0] = malloc(w*h*sizeof(int));
	memset(image.bad[0], 0, w*h*sizeof(int));

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
		cell = cell_new();
		cell_set_lattice_type(cell, L_CUBIC);
		cell_set_centering(cell, 'P');
		cell_set_parameters(cell, 800.0e-10, 800.0e-10, 800.0e-10,
		                    deg2rad(90.0), deg2rad(90.0), deg2rad(90.0));
		cell = cell_rotate(cell, random_quaternion(rng));

		ic.halfw = ir_out;
		ic.image = &image;
		ic.k = 1.0/image.lambda;
		ic.n_saturated = 0;
		ic.n_implausible = 0;
		ic.cell = cell;
		ic.ir_inn = ir_inn;
		ic.ir_mid = ir_mid;
		ic.ir_out = ir_out;
		ic.meth = INTEGRATION_RINGS;
		ic.int_diag = INTDIAG_NONE;
		ic.masks = NULL;
		if ( init_intcontext(&ic) ) {
			ERROR("Failed to initialise integration.\n");
			return 1;
		}
		setup_ring_masks(&ic, ir_inn, ir_mid, ir_out);

		integrate_rings_once(refl, &image, &ic, cell, 0);

		cell_free(cell);

		histogram_add_value(hi, get_intensity(refl));
		esd_sum += get_esd_intensity(refl);

	}
	printf("Mean calculated sigma(I) = %.2f\n", esd_sum / 300);

	histogram_show(hi);

	histogram_free(hi);
	free(image.beam);
	free(image.det->panels);
	free(image.det);
	free(image.dp[0]);
	free(image.dp);

	if ( fail ) return 1;

	return 0;
}
