/*
 * polarisation_check.c
 *
 * Check polarisation correction
 *
 * Copyright © 2019 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019 Thomas White <taw@physics.org>
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
#include <cell.h>
#include <cell-utils.h>
#include <geometry.h>


int main(int argc, char *argv[])
{
	struct image image;
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
	gsl_rng *rng;
	struct polarisation p;
	int i;
	double *map;
	double *nmap;
	const int ntrial = 1000;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	fh = fopen("/dev/urandom", "r");
	if ( fread(&seed, sizeof(seed), 1, fh) == 1 ) {
		gsl_rng_set(rng, seed);
	} else {
		ERROR("Failed to seed RNG\n");
	}
	fclose(fh);

	image.beam = NULL;
	image.lambda = ph_eV_to_lambda(9000.0);
	image.bw = 0.000001;
	image.div = 0.0;
	image.spectrum = spectrum_generate_gaussian(image.lambda, image.bw);

	image.det = calloc(1, sizeof(struct detector));
	image.det->n_panels = 1;
	image.det->panels = calloc(1, sizeof(struct panel));

	image.dp = calloc(1, sizeof(float *));
	image.bad = calloc(1, sizeof(int *));

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
	image.det->panels[0].clen = 50.0e-3;
	image.det->panels[0].res = 10000;  /* 10 micron pixels */
	image.det->panels[0].adu_per_eV = 10.0/9000.0; /* 10 adu/ph */
	image.det->panels[0].max_adu = +INFINITY;  /* No cutoff */
	image.det->panels[0].orig_min_fs = 0;
	image.det->panels[0].orig_min_ss = 0;
	image.det->panels[0].orig_max_fs = w-1;
	image.det->panels[0].orig_max_ss = h-1;

	image.det->furthest_out_panel = &image.det->panels[0];
	image.det->furthest_out_fs = 0;
	image.det->furthest_out_ss = 0;

	image.dp[0] = malloc(w*h*sizeof(float));
	memset(image.dp[0], 0, w*h*sizeof(float));
	image.bad[0] = malloc(w*h*sizeof(int));
	memset(image.bad[0], 0, w*h*sizeof(int));
	image.sat = NULL;

	cell = cell_new();
	cell_set_lattice_type(cell, L_CUBIC);
	cell_set_centering(cell, 'P');
	cell_set_parameters(cell, 50.0e-10, 50.0e-10, 50.0e-10,
	                    deg2rad(90.0), deg2rad(90.0), deg2rad(90.0));

	cr = crystal_new();
	crystal_set_profile_radius(cr, 0.001e9);
	crystal_set_mosaicity(cr, 0.0);  /* radians */
	crystal_set_image(cr, &image);
	crystal_set_cell(cr, cell);

	image.n_crystals = 1;
	image.crystals = &cr;

	map = malloc(w*h*sizeof(double));
	nmap = malloc(w*h*sizeof(double));
	for ( i=0; i<w*h; i++ ) {
		map[i] = 0.0;
		nmap[i] = 0.0;
	}
	for ( i=0; i<ntrial; i++ ) {

		UnitCell *ncell;

		ncell = cell_rotate(cell, random_quaternion(rng));
		crystal_set_cell(cr, ncell);

		list = predict_to_res(cr, largest_q(&image));
		crystal_set_reflections(cr, list);

		for ( refl = first_refl(list, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			set_intensity(refl, 1.0);
		}

		p.angle = deg2rad(105.0);
		p.fraction = 1.0;
		polarisation_correction(list, ncell, p);

		for ( refl = first_refl(list, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			double fs, ss;
			int nfs, nss;
			get_detector_pos(refl, &fs, &ss);
			nfs = fs;  nss = ss;  /* Explicit truncation */

			/* Intensity in reflist is corrected,
			 * but we want "un-correction" */
			map[nfs + nss*w] += 1.0/get_intensity(refl);
			nmap[nfs + nss*w] += 1.0;
		}

		cell_free(ncell);
		reflist_free(list);

		progress_bar(i+1, ntrial, "Calculating");

	}

	for ( i=0; i<w*h; i++ ) {
		image.dp[0][i] = 1000.0 * map[i] / nmap[i];
		if ( isnan(image.dp[0][i]) ) image.dp[0][i] = 0.0;
	}
	hdf5_write_image("test.h5", &image, "/data/data");

	free(image.beam);
	free(image.det->panels);
	free(image.det);
	free(image.dp[0]);
	free(image.dp);
	gsl_rng_free(rng);

	if ( fail ) return 1;

	return 0;
}
