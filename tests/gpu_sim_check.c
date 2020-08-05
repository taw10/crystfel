/*
 * gpu_sim_check.c
 *
 * Check that GPU simulation agrees with CPU version
 *
 * Copyright © 2012-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2012-2019 Thomas White <taw@physics.org>
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

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#define CL_TARGET_OPENCL_VERSION 220

#include "../src/diffraction.h"
#include "../src/diffraction-gpu.h"
#include "../src/cl-utils.h"

#include <datatemplate.h>
#include <utils.h>
#include <image.h>
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
	struct image *gpu_image;
	struct image *cpu_image;
	DataTemplate *dtempl;
	UnitCell *cell;
	UnitCell *cell_raw;
	int i;
	double gpu_min, gpu_max, gpu_tot;
	double cpu_min, cpu_max, cpu_tot;
	double dev, perc;
	double start, end;
	double gpu_time, cpu_time;
	SymOpList *sym;
	gsl_rng *rng;

	if ( have_gpu_device() == 0 ) {
		ERROR("No GPU device found - skipping test.\n");
		return 0;
	}

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	gctx = setup_gpu(1, NULL, NULL, NULL, 0);
	if ( gctx == NULL ) {
		ERROR("Couldn't set up GPU.\n");
		return 1;
	}

	cell_raw = cell_new_from_parameters(28.1e-9, 28.1e-9, 16.5e-9,
	                          deg2rad(90.0), deg2rad(90.0), deg2rad(120.0));

	cell = cell_rotate(cell_raw, random_quaternion(rng));

	dtempl = data_template_new_from_file(argv[1]);
	if ( dtempl == NULL ) return 1;

	cpu_image = image_create_for_simulation(dtempl);
	gpu_image = image_create_for_simulation(dtempl);

	start = get_hires_seconds();
	if ( get_diffraction_gpu(gctx, gpu_image, 8, 8, 8, cell, 0, 0, 10) ) {
		return 1;
	}
	end = get_hires_seconds();
	gpu_time = end - start;

	sym = get_pointgroup("1");

	cpu_image->dp = malloc(cpu_image->detgeom->n_panels * sizeof(float *));
	if ( cpu_image->dp == NULL ) {
		ERROR("Couldn't allocate memory for result.\n");
		return 1;
	}
	for ( i=0; i<cpu_image->detgeom->n_panels; i++ ) {
		struct detgeom_panel *p = &cpu_image->detgeom->panels[i];
		cpu_image->dp[i] = calloc(p->w * p->h, sizeof(float));
		if ( cpu_image->dp[i] == NULL ) {
			ERROR("Couldn't allocate memory for panel %i\n", i);
			return 1;
		}
	}

	start = get_hires_seconds();
	get_diffraction(cpu_image, 8, 8, 8, NULL, NULL, NULL, cell,
	                GRADIENT_MOSAIC, sym, 0, 0, 10);
	end = get_hires_seconds();
	cpu_time = end - start;

	free_symoplist(sym);

	STATUS("The GPU version was %5.2f times faster.\n", cpu_time/gpu_time);

	gpu_min = +INFINITY;  gpu_max = -INFINITY;  gpu_tot = 0.0;
	cpu_min = +INFINITY;  cpu_max = -INFINITY;  cpu_tot = 0.0;
	dev = 0.0;
	for ( i=0; i<cpu_image->detgeom->n_panels; i++ ) {

		int j;
		struct detgeom_panel *p = &cpu_image->detgeom->panels[i];

		for ( j=0; j<p->w*p->h; j++ ) {

			const double cpu = cpu_image->dp[i][j];
			const double gpu = gpu_image->dp[i][j];

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

	if ( perc > 1.2 ) {

		STATUS("Test failed!  I'm writing cpu-sim.h5 and gpu-sim.h5"
		       " for you to inspect.\n");

		image_write(cpu_image, dtempl, "cpu-sim.h5");
		image_write(gpu_image, dtempl, "gpu-sim.h5");

		return 1;

	}

	gsl_rng_free(rng);
	cell_free(cell);
	image_free(cpu_image);
	image_free(gpu_image);
	data_template_free(dtempl);

	return 0;
}
