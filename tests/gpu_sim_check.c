/*
 * gpu_sim_check.c
 *
 * Check that GPU simulation agrees with CPU version
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2015 Thomas White <taw@physics.org>
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

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "../src/diffraction.h"
#include "../src/diffraction-gpu.h"
#include <detector.h>
#include <utils.h>
#include <symmetry.h>
#include <cell-utils.h>


#ifdef HAVE_CLOCK_GETTIME

static double get_hires_seconds()
{
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);
	return (double)tp.tv_sec + ((double)tp.tv_nsec/1e9);
}

#else

/* Fallback version of the above.  The time according to gettimeofday() is not
 * monotonic, so measuring intervals based on it will screw up if there's a
 * timezone change (e.g. daylight savings) while the program is running. */
static double get_hires_seconds()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double)tp.tv_sec + ((double)tp.tv_usec/1e6);
}

#endif


int main(int argc, char *argv[])
{
	struct gpu_context *gctx;
	struct image gpu_image;
	struct image cpu_image;
	UnitCell *cell;
	UnitCell *cell_raw;
	struct detector *det;
	int i;
	double gpu_min, gpu_max, gpu_tot;
	double cpu_min, cpu_max, cpu_tot;
	double dev, perc;
	const double sep = 20.0;
	double start, end;
	double gpu_time, cpu_time;
	SymOpList *sym;
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	gctx = setup_gpu(1, NULL, NULL, NULL, 0);
	if ( gctx == NULL ) {
		ERROR("Couldn't set up GPU.\n");
		return 1;
	}

	cell_raw = cell_new_from_parameters(28.1e-9, 28.1e-9, 16.5e-9,
	                          deg2rad(90.0), deg2rad(90.0), deg2rad(120.0));

	cell = cell_rotate(cell_raw, random_quaternion(rng));

	det = calloc(1, sizeof(struct detector));
	det->n_panels = 2;
	det->panels = calloc(2, sizeof(struct panel));

	det->panels[0].orig_min_fs = 0;
	det->panels[0].orig_max_fs = 1023;
	det->panels[0].orig_min_ss = 0;
	det->panels[0].orig_max_ss = 511;
	det->panels[0].w = 1024;
	det->panels[0].h = 512;
	det->panels[0].fsx = 1.0;
	det->panels[0].fsy = 0.0;
	det->panels[0].fsz = 0.4;
	det->panels[0].ssx = 0.0;
	det->panels[0].ssy = 1.0;
	det->panels[0].ssz = 0.0;
	det->panels[0].xfs = 1.0;
	det->panels[0].yfs = 0.0;
	det->panels[0].xss = 0.0;
	det->panels[0].yss = 1.0;
	det->panels[0].cnx = -512.0;
	det->panels[0].cny = -512.0-sep;
	det->panels[0].clen = 100.0e-3;
	det->panels[0].res = 9090.91;
	det->panels[0].adu_per_eV = 1.0;
	det->panels[0].data = NULL;

	det->panels[1].orig_min_fs = 0;
	det->panels[1].orig_max_fs = 1023;
	det->panels[1].orig_min_ss = 512;
	det->panels[1].orig_max_ss = 1023;
	det->panels[1].w = 1024;
	det->panels[1].h = 512;
	det->panels[1].fsx = 1.0;
	det->panels[1].fsy = 0.0;
	det->panels[1].fsz = 0.0;
	det->panels[1].ssx = 0.0;
	det->panels[1].ssy = 1.0;
	det->panels[1].ssz = 1.4;
	det->panels[1].xfs = 1.0;
	det->panels[1].yfs = 0.0;
	det->panels[1].xss = 0.0;
	det->panels[1].yss = 1.0;
	det->panels[1].cnx = -512.0;
	det->panels[1].cny = sep;
	det->panels[1].clen = 100.0e-3;
	det->panels[1].res = 9090.91;
	det->panels[1].adu_per_eV = 1.0;
	det->panels[1].data = NULL;

	cpu_image.det = det;
	gpu_image.det = det;
	cpu_image.beam = NULL;
	gpu_image.beam = NULL;
	cpu_image.spectrum = NULL;
	gpu_image.spectrum = NULL;

	cpu_image.lambda = ph_en_to_lambda(eV_to_J(6000));
	gpu_image.lambda = ph_en_to_lambda(eV_to_J(6000));
	cpu_image.bw = 1.0 / 100.0;
	gpu_image.bw = 1.0 / 100.0;

	cpu_image.nsamples = 10;
	gpu_image.nsamples = 10;
	cpu_image.spectrum = generate_tophat(&cpu_image);
	gpu_image.spectrum = generate_tophat(&gpu_image);

	start = get_hires_seconds();
	if ( get_diffraction_gpu(gctx, &gpu_image, 8, 8, 8, cell, 1) ) {
		return 1;
	}
	end = get_hires_seconds();
	gpu_time = end - start;

	sym = get_pointgroup("1");

	cpu_image.dp = malloc(det->n_panels * sizeof(float *));
	if ( cpu_image.dp == NULL ) {
		ERROR("Couldn't allocate memory for result.\n");
		return 1;
	}
	for ( i=0; i<det->n_panels; i++ ) {
		struct panel *p = &det->panels[i];
		cpu_image.dp[i] = calloc(p->w * p->h, sizeof(float));
		if ( cpu_image.dp[i] == NULL ) {
			ERROR("Couldn't allocate memory for panel %i\n", i);
			return 1;
		}
	}

	start = get_hires_seconds();
	get_diffraction(&cpu_image, 8, 8, 8, NULL, NULL, NULL, cell,
	                GRADIENT_MOSAIC, sym, 1);
	end = get_hires_seconds();
	cpu_time = end - start;

	free_symoplist(sym);

	STATUS("The GPU version was %5.2f times faster.\n", cpu_time/gpu_time);

	gpu_min = +INFINITY;  gpu_max = -INFINITY;  gpu_tot = 0.0;
	cpu_min = +INFINITY;  cpu_max = -INFINITY;  cpu_tot = 0.0;
	dev = 0.0;
	for ( i=0; i<det->n_panels; i++ ) {

		int j;
		struct panel *p = &det->panels[i];

		for ( j=0; j<p->w*p->h; j++ ) {

			const double cpu = cpu_image.dp[i][j];
			const double gpu = gpu_image.dp[i][j];

			if ( cpu > cpu_max ) cpu_max = cpu;
			if ( cpu < cpu_min ) cpu_min = cpu;
			if ( gpu > gpu_max ) gpu_max = gpu;
			if ( gpu < gpu_min ) gpu_min = gpu;
			gpu_tot += gpu;
			cpu_tot += cpu;
			dev += fabs(gpu - cpu);

		}

	}
	perc = 100.0*dev/cpu_tot;

	STATUS("GPU: min=%8e, max=%8e, total=%8e\n", gpu_min, gpu_max, gpu_tot);
	STATUS("CPU: min=%8e, max=%8e, total=%8e\n", cpu_min, cpu_max, cpu_tot);
	STATUS("dev = %8e (%5.2f%% of CPU total)\n", dev, perc);

	if ( perc > 1.0 ) {

		STATUS("Test failed!  I'm writing cpu-sim.h5 and gpu-sim.h5"
		       " for you to inspect.\n");

		hdf5_write_image("cpu-sim.h5", &cpu_image, NULL);
		hdf5_write_image("gpu-sim.h5", &gpu_image, NULL);

		return 1;

	}

	gsl_rng_free(rng);
	cell_free(cell);
	free_detector_geometry(det);


	return 0;
}
