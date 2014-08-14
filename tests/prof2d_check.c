/*
 * prof2d_check.c
 *
 * Check 2D profile fitting
 *
 * Copyright © 2013-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2013-2014 Thomas White <taw@physics.org>
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


#define ADD_PX(fs, ss, val) \
	if ( ((fs)>0) && ((ss)>0) && ((fs)<w) && ((ss)<h) ) { \
		image.dp[0][(fs)+w*(ss)] += (val); \
	}

int main(int argc, char *argv[])
{
	struct image image;
	int fs, ss;
	FILE *fh;
	unsigned long int seed;
	int fail = 0;
	const int w = 1024;
	const int h = 1024;
	RefList *list;
	RefListIterator *iter;
	Reflection *refl;
	UnitCell *cell;
	Crystal *cr;
	const int ir_inn = 2;
	const int ir_mid = 4;
	const int ir_out = 6;
	Histogram *hi;
	double esd_sum = 0.0;
	int n = 0;
	int n_strong = 0;
	int n_weak = 0;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	fh = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fh);
	fclose(fh);
	gsl_rng_set(rng, seed);

	image.flags = NULL;
	image.beam = NULL;
	image.lambda = ph_eV_to_lambda(9000.0);
	image.bw = 0.000001;
	image.div = 0.0;

	image.det = calloc(1, sizeof(struct detector));
	image.det->n_panels = 1;
	image.det->panels = calloc(1, sizeof(struct panel));

	image.dp = calloc(1, sizeof(float *));
	image.bad = calloc(1, sizeof(int *));

	image.width = w;
	image.height = h;
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

	image.det->furthest_out_panel = &image.det->panels[0];
	image.det->furthest_out_fs = 0;
	image.det->furthest_out_ss = 0;

	image.dp[0] = malloc(w*h*sizeof(float));
	memset(image.dp[0], 0, w*h*sizeof(float));
	image.bad[0] = malloc(w*h*sizeof(int));
	memset(image.bad[0], 0, w*h*sizeof(int));

	cell = cell_new();
	cell_set_lattice_type(cell, L_CUBIC);
	cell_set_centering(cell, 'P');
	cell_set_parameters(cell, 800.0e-10, 800.0e-10, 800.0e-10,
	                    deg2rad(90.0), deg2rad(90.0), deg2rad(90.0));
	cell = cell_rotate(cell, random_quaternion(rng));

	cr = crystal_new();
	crystal_set_profile_radius(cr, 0.001e9);
	crystal_set_mosaicity(cr, 0.0);  /* radians */
	crystal_set_image(cr, &image);
	crystal_set_cell(cr, cell);

	image.n_crystals = 1;
	image.crystals = &cr;

	list = find_intersections(&image, cr, PMODEL_SCSPHERE);
	crystal_set_reflections(cr, list);

	for ( fs=0; fs<w; fs++ ) {
	for ( ss=0; ss<h; ss++ ) {
		image.dp[0][fs+w*ss] = 10.0*poisson_noise(rng, 40);
	}
	}

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		double pfs, pss;
		int fs, ss;
		signed int hj, kj, lj;

		get_detector_pos(refl, &pfs, &pss);
		fs = pfs;  ss = pss;

		get_indices(refl, &hj, &kj, &lj);
		if ( lj % 2 ) {
			const int pk_ph = 1000;
			ADD_PX(fs, ss, 10.0*poisson_noise(rng, pk_ph));
			ADD_PX(fs-1, ss, 10.0*poisson_noise(rng, pk_ph));
			ADD_PX(fs+1, ss, 10.0*poisson_noise(rng, pk_ph));
			ADD_PX(fs, ss-1, 10.0*poisson_noise(rng, pk_ph));
			ADD_PX(fs, ss+1, 10.0*poisson_noise(rng, pk_ph));
			n_strong++;
		} else {
			/* Absent peak */
			n_weak++;
		}
	}

	STATUS("%i strong, %i weak\n", n_strong, n_weak);

	integrate_prof2d(INTEGRATION_PROF2D, cr, &image,
	                 INTDIAG_NONE, 0, 0, 0, ir_inn, ir_mid, ir_out, 0,
	                 NULL);

	printf("Weak reflections:\n");
	hi = histogram_init();
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		if ( get_redundancy(refl) == 0 ) continue;
		if ( l % 2 ) continue;  /* Ignore strong reflections */

		histogram_add_value(hi, get_intensity(refl));
		esd_sum += get_esd_intensity(refl);
		n++;

	}
	printf("Mean calculated sigma(I) = %.2f (%i measurements)\n",
	       esd_sum / n, n);
	histogram_show(hi);
	histogram_free(hi);

	printf("Strong reflections:\n");
	hi = histogram_init();
	esd_sum = 0.0;
	n = 0;
	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		if ( get_redundancy(refl) == 0 ) continue;
		if ( l % 2 == 0 ) continue;  /* Ignore weak reflections */

		histogram_add_value(hi, get_intensity(refl));
		esd_sum += get_esd_intensity(refl);
		n++;

	}
	printf("Mean calculated sigma(I) = %.2f (%i measurements)\n",
	       esd_sum / n, n);
	histogram_show(hi);
	histogram_free(hi);

	free(image.beam);
	free(image.det->panels);
	free(image.det);
	free(image.dp[0]);
	free(image.dp);
	gsl_rng_free(rng);

	if ( fail ) return 1;

	return 0;
}
