/*
 * gpu_sim_check.c
 *
 * Check that GPU simulation agrees with CPU version
 *
 * (c) 2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>

#include "../src/diffraction.h"
#include "../src/diffraction-gpu.h"
#include "../src/detector.h"
#include "../src/beam-parameters.h"
#include "../src/utils.h"


int main(int argc, char *argv[])
{
	struct gpu_context *gctx;
	struct image gpu_image;
	struct image cpu_image;
	UnitCell *cell;
	UnitCell *cell_raw;
	struct detector *det;
	struct beam_params *beam;
	int i;
	double gpu_min, gpu_max, gpu_tot;
	double cpu_min, cpu_max, cpu_tot;
	double dev, perc;
	const double sep = 20.0;

	gctx = setup_gpu(1, NULL, NULL, NULL, 0);
	if ( gctx == NULL ) {
		ERROR("Couldn't set up GPU.\n");
		return 1;
	}

	cell_raw = cell_new_from_parameters(28.1e-9, 28.1e-9, 16.5e-9,
	                          deg2rad(90.0), deg2rad(90.0), deg2rad(120.0));

	cell = cell_rotate(cell_raw, random_quaternion());

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

	cpu_image.det = det;
	gpu_image.det = det;

	beam = calloc(1, sizeof(struct beam_params));
	beam->fluence = 1.0e15;
	beam->beam_radius = 1.0e-6;
	beam->photon_energy = 9000.0;
	beam->bandwidth = 0.1 / 100.0;
	beam->divergence = 0.0;
	beam->dqe = 1.0;
	beam->adu_per_photon = 1.0;
	cpu_image.beam = beam;
	gpu_image.beam = beam;

	cpu_image.lambda = ph_en_to_lambda(eV_to_J(beam->photon_energy));
	gpu_image.lambda = ph_en_to_lambda(eV_to_J(beam->photon_energy));

	get_diffraction_gpu(gctx, &gpu_image, 8, 8, 8, cell);
	get_diffraction(&cpu_image, 8, 8, 8, NULL, NULL, NULL, cell,
	                GRADIENT_MOSAIC, "1");

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
	free(beam);

	if ( perc > 1.0 ) return 1;

	return 0;
}
