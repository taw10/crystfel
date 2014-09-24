/*
 * gpu_sim_check.c
 *
 * Check that GPU simulation agrees with CPU version
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2014 Thomas White <taw@physics.org>
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

	gpu_image.width = 1024;
	gpu_image.height = 1024;
	cpu_image.width = 1024;
	cpu_image.height = 1024;
	det = calloc(1, sizeof(struct detector));
	det->n_panels = 2;
	det->panels = calloc(2, sizeof(struct panel));

	det->panels[0].min_fs = 0;
	det->panels[0].max_fs = 1023;
	det->panels[0].min_ss = 0;
	det->panels[0].max_ss = 511;
	det->panels[0].fsx = 1;
	det->panels[0].fsy = 0;
	det->panels[0].ssx = 0;
	det->panels[0].ssy = 1;
	det->panels[0].xfs = 1;
	det->panels[0].yfs = 0;
	det->panels[0].xss = 0;
	det->panels[0].yss = 1;
	det->panels[0].cnx = -512.0;
	det->panels[0].cny = -512.0-sep;
	det->panels[0].clen = 100.0e-3;
	det->panels[0].res = 9090.91;
	det->panels[0].adu_per_eV = 1.0;

	det->panels[1].min_fs = 0;
	det->panels[1].max_fs = 1023;
	det->panels[1].min_ss = 512;
	det->panels[1].max_ss = 1023;
	det->panels[1].fsx = 1;
	det->panels[1].fsy = 0;
	det->panels[1].ssx = 0;
	det->panels[1].ssy = 1;
	det->panels[1].xfs = 1;
	det->panels[1].yfs = 0;
	det->panels[1].xss = 0;
	det->panels[1].yss = 1;
	det->panels[1].cnx = -512.0;
	det->panels[1].cny = sep;
	det->panels[1].clen = 100.0e-3;
	det->panels[1].res = 9090.91;
	det->panels[0].adu_per_eV = 1.0;

	cpu_image.det = det;
	gpu_image.det = det;
	cpu_image.beam = NULL;
	gpu_image.beam = NULL;

	cpu_image.lambda = ph_en_to_lambda(eV_to_J(6000));
	gpu_image.lambda = ph_en_to_lambda(eV_to_J(6000));
	cpu_image.bw = 1.0 / 100.0;
	gpu_image.bw = 1.0 / 100.0;

	cpu_image.nsamples = 10;
	gpu_image.nsamples = 10;
	cpu_image.spectrum = generate_tophat(&cpu_image);
	gpu_image.spectrum = generate_tophat(&gpu_image);

	start = get_hires_seconds();
	get_diffraction_gpu(gctx, &gpu_image, 8, 8, 8, cell, 1);
	end = get_hires_seconds();
	gpu_time = end - start;

	sym = get_pointgroup("1");

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
	for ( i=0; i<1024*1024; i++ ) {

		const double cpu = cpu_image.data[i];
		const double gpu = gpu_image.data[i];

		if ( cpu > cpu_max ) cpu_max = cpu;
		if ( cpu < cpu_min ) cpu_min = cpu;
		if ( gpu > gpu_max ) gpu_max = gpu;
		if ( gpu < gpu_min ) gpu_min = gpu;
		gpu_tot += gpu;
		cpu_tot += cpu;
		dev += fabs(gpu - cpu);

	}
	perc = 100.0*dev/cpu_tot;

	STATUS("GPU: min=%8e, max=%8e, total=%8e\n", gpu_min, gpu_max, gpu_tot);
	STATUS("CPU: min=%8e, max=%8e, total=%8e\n", cpu_min, cpu_max, cpu_tot);
	STATUS("dev = %8e (%5.2f%% of CPU total)\n", dev, perc);

	cell_free(cell);
	free_detector_geometry(det);

	if ( perc > 1.0 ) {

		STATUS("Test failed!  I'm writing cpu-sim.h5 and gpu-sim.h5"
		       " for you to inspect.\n");

		hdf5_write("cpu-sim.h5", cpu_image.data, cpu_image.width,
		            cpu_image.height, H5T_NATIVE_FLOAT);

		hdf5_write("gpu-sim.h5", gpu_image.data, gpu_image.width,
		            gpu_image.height, H5T_NATIVE_FLOAT);

		return 1;

	}

	gsl_rng_free(rng);

	return 0;
}
