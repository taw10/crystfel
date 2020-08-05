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
	DataTemplate *dtempl;
	struct image *image;
	FILE *fh;
	unsigned long int seed;
	int fail = 0;
	RefList *list;
	RefListIterator *iter;
	Reflection *refl;
	UnitCell *cell;
	Crystal *cr;
	gsl_rng *rng;
	struct polarisation p;
	int i;
	double *map;
	int *nmap;
	const int ntrial = 1000;

	/* NB must match polarisation_check.geom */
	const int w = 512;
	const int h = 512;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	fh = fopen("/dev/urandom", "r");
	if ( fread(&seed, sizeof(seed), 1, fh) == 1 ) {
		gsl_rng_set(rng, seed);
	} else {
		ERROR("Failed to seed RNG\n");
	}
	fclose(fh);

	dtempl = data_template_new_from_file(argv[1]);
	if ( dtempl == NULL ) return 1;

	image = image_create_for_simulation(dtempl);
	if ( image == NULL ) return 1;

	cell = cell_new();
	cell_set_lattice_type(cell, L_CUBIC);
	cell_set_centering(cell, 'P');
	cell_set_parameters(cell, 50.0e-10, 50.0e-10, 50.0e-10,
	                    deg2rad(90.0), deg2rad(90.0), deg2rad(90.0));

	cr = crystal_new();
	crystal_set_profile_radius(cr, 0.001e9);
	crystal_set_mosaicity(cr, 0.0);  /* radians */
	crystal_set_image(cr, image);
	crystal_set_cell(cr, cell);

	image->n_crystals = 1;
	image->crystals = &cr;

	map = malloc(w*h*sizeof(double));
	nmap = malloc(w*h*sizeof(int));
	for ( i=0; i<w*h; i++ ) {
		map[i] = 0.0;
		nmap[i] = 0;
	}
	for ( i=0; i<ntrial; i++ ) {

		UnitCell *ncell;

		ncell = cell_rotate(cell, random_quaternion(rng));
		crystal_set_cell(cr, ncell);

		list = predict_to_res(cr, detgeom_max_resolution(image->detgeom,
		                                                 image->lambda));
		crystal_set_reflections(cr, list);

		for ( refl = first_refl(list, &iter);
		      refl != NULL;
		      refl = next_refl(refl, iter) )
		{
			set_intensity(refl, 1.0);
		}

		p.angle = deg2rad(105.0);
		p.fraction = 1.0;
		p.disable = 0;
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
			nmap[nfs + nss*w] += 1;
		}

		cell_free(ncell);
		reflist_free(list);

		crystal_set_cell(cr, NULL);
		crystal_set_reflections(cr, NULL);

		progress_bar(i+1, ntrial, "Calculating");

	}

	for ( i=0; i<w*h; i++ ) {
		if ( nmap[i] > 0 ) {
			image->dp[0][i] = 1000.0 * map[i] / nmap[i];
		}
	}
	image_write(image, dtempl, "test.h5");

	image->crystals = NULL;
	image->n_crystals = 0;

	data_template_free(dtempl);
	image_free(image);
	gsl_rng_free(rng);

	if ( fail ) return 1;

	return 0;
}
