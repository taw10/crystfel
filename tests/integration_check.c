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
#include <beam-parameters.h>

#include "../libcrystfel/src/integration.c"



int main(int argc, char *argv[])
{
	struct image image;
	double fsp, ssp, intensity, sigma;
	int fs, ss;
	FILE *fh;
	unsigned int seed;
	int fail = 0;
	int r, npx;
	double ex;

	fh = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fh);
	fclose(fh);
	srand(seed);

	image.data = malloc(128*128*sizeof(float));
	image.flags = NULL;
	image.beam = NULL;
	image.lambda = ph_eV_to_lambda(1000.0);

	image.det = calloc(1, sizeof(struct detector));
	image.det->n_panels = 1;
	image.det->panels = calloc(1, sizeof(struct panel));

	image.det->panels[0].min_fs = 0;
	image.det->panels[0].max_fs = 128;
	image.det->panels[0].min_ss = 0;
	image.det->panels[0].max_ss = 128;
	image.det->panels[0].fsx = 1.0;
	image.det->panels[0].fsy = 0.0;
	image.det->panels[0].ssx = 0.0;
	image.det->panels[0].ssy = 1.0;
	image.det->panels[0].xfs = 1.0;
	image.det->panels[0].yfs = 0.0;
	image.det->panels[0].xss = 0.0;
	image.det->panels[0].yss = 1.0;
	image.det->panels[0].cnx = -64.0;
	image.det->panels[0].cny = -64.0;
	image.det->panels[0].clen = 1.0;
	image.det->panels[0].res = 1.0;
	image.det->panels[0].adu_per_eV = 1.0/1000.0;  /* -> 1 adu per photon */
	image.det->panels[0].max_adu = +INFINITY;  /* No cutoff */

	image.width = 128;
	image.height = 128;
	memset(image.data, 0, 128*128*sizeof(float));

	image.n_crystals = 0;
	image.crystals = NULL;

	/* Uniform peak on uniform background */
	npx = 0;
	for ( fs=0; fs<image.width; fs++ ) {
	for ( ss=0; ss<image.height; ss++ ) {
		image.data[fs+image.width*ss] = 1000.0;
		if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) > 9*9 ) continue;
		image.data[fs+image.width*ss] += 1000.0;
		npx++;
	}
	}

	r = integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity, &sigma,
	                   10.0, 15.0, 17.0, NULL, NULL);
	if ( r ) {
		ERROR("   Fifth check: integrate_peak() returned %i (wrong).\n",
		      r);
		fail = 1;
	} else {

		STATUS("  Fifth check: intensity = %.2f, sigma = %.2f\n",
		       intensity, sigma);

		ex = npx*1000.0;
		if ( within_tolerance(ex, intensity, 1.0) == 0 ) {
			ERROR("Intensity should be close to %f\n", ex);
			fail = 1;
		}

		ex = sqrt(npx*1000.0);
		if ( within_tolerance(ex, sigma, 1.0) == 0 ) {
			ERROR("Sigma should be roughly %f.\n", ex);
			fail = 1;
		}

	}


	free(image.beam);
	free(image.det->panels);
	free(image.det);
	free(image.data);

	if ( fail ) return 1;

	return 0;
}
