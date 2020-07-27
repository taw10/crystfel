/*
 * polarisation_check.c
 *
 * Check polarisation correction
 *
 * Copyright © 2019-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2019-2020 Thomas White <taw@physics.org>
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

	image.lambda = ph_eV_to_lambda(9000.0);
	image.bw = 0.000001;
	image.div = 0.0;
	image.spectrum = spectrum_generate_gaussian(image.lambda, image.bw);

	image.detgeom = calloc(1, sizeof(struct detgeom));
	image.detgeom->n_panels = 1;
	image.detgeom->panels = calloc(1, sizeof(struct detgeom_panel));

	image.dp = calloc(1, sizeof(float *));
	image.bad = calloc(1, sizeof(int *));

	image.detgeom->panels[0].w = w;
	image.detgeom->panels[0].h = h;
	image.detgeom->panels[0].fsx = 1.0;
	image.detgeom->panels[0].fsy = 0.0;
	image.detgeom->panels[0].ssx = 0.0;
	image.detgeom->panels[0].ssy = 1.0;
	image.detgeom->panels[0].cnx = -w/2;
	image.detgeom->panels[0].cny = -h/2;
	image.detgeom->panels[0].cnz = 50.0e-3 / 10e-6;
	image.detgeom->panels[0].pixel_pitch = 10e-6;
	image.detgeom->panels[0].adu_per_photon = 1.0;
	image.detgeom->panels[0].max_adu = +INFINITY;  /* No cutoff */

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

		list = predict_to_res(cr, detgeom_max_resolution(image.detgeom,
		                                                 image.lambda));
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
	//hdf5_write_image("test.h5", &image, "/data/data");

	detgeom_free(image.detgeom);
	free(image.dp[0]);
	free(image.dp);
	gsl_rng_free(rng);

	if ( fail ) return 1;

	return 0;
}
