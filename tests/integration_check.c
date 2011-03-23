/*
 * integration_check.c
 *
 * Check peak integration
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

#include "../src/image.h"
#include "../src/peaks.h"
#include "../src/utils.h"
#include "../src/beam-parameters.h"


int main(int argc, char *argv[])
{
	struct image image;
	double fsp, ssp, intensity, bg, max, sigma;
	int fs, ss;

	image.data = malloc(128*128*sizeof(float));
	image.flags = NULL;

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
	image.det->panels[0].integr_radius = 10.0;

	image.width = 128;
	image.height = 128;
	memset(image.data, 0, 128*128*sizeof(float));

	image.beam = calloc(1, sizeof(struct beam_params));
	image.beam->adu_per_photon = 100.0;

	/* First check: no intensity -> zero intensity and bg */
	integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity, &bg, &max, &sigma, 0, 1);
	STATUS(" First check: intensity = %.2f bg = %.2f max = %.2f\n",
	       intensity, bg, max);
	if ( intensity != 0.0 ) {
		ERROR("Intensity should be zero.\n");
		return 1;
	}
	if ( bg != 0.0 ) {
		ERROR("Background should be zero\n");
		return 1;
	}

	/* Second check: uniform peak gives correct value and no bg */
	for ( fs=0; fs<image.width; fs++ ) {
	for ( ss=0; ss<image.height; ss++ ) {
		if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) > 9*9 ) continue;
		image.data[fs+image.width*ss] = 1000.0;
	}
	}
	integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity, &bg, &max, &sigma, 0, 1);
	STATUS("Second check: intensity = %.2f, bg = %.2f, max = %.2f\n",
	       intensity, bg, max);
	if ( fabs(intensity - M_PI*9.0*9.0*1000.0) > 4000.0 ) {
		ERROR("Intensity should be close to 1000*pi*integr_r^2\n");
		return 1;
	}
	if ( bg != 0.0 ) {
		ERROR("Background should be zero\n");
		return 1;
	}
	if ( max != 1000.0 ) {
		ERROR("Max should be 1000\n");
		return 1;
	}

	/* Third check: flat Poisson background should hoover up bg only */
	for ( fs=0; fs<image.width; fs++ ) {
	for ( ss=0; ss<image.height; ss++ ) {
		image.data[fs+image.width*ss] = poisson_noise(10.0) - 10.0;
	}
	}
	integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity, &bg, &max, &sigma, 0, 1);
	STATUS(" Third check: intensity = %.2f, bg = %.2f, max = %.2f\n",
	       intensity, bg, max);
	if ( fabs(intensity) > 100.0 ) {
		ERROR("Intensity should be close to zero\n");
		return 1;
	}
	if ( fabs(bg - sqrt(10.0)) > 1.0 ) {
		ERROR("Background should be close to sqrt(10)\n");
		return 1;
	}

	for ( fs=0; fs<image.width; fs++ ) {
	for ( ss=0; ss<image.height; ss++ ) {
		if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) > 9*9 ) continue;
		image.data[fs+image.width*ss] = 1000.0;
	}
	}
	integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity, &bg, &max, &sigma, 0, 1);
	STATUS("Fourth check: intensity = %.2f, bg = %.2f, max = %.2f\n",
	       intensity, bg, max);
	if ( fabs(intensity - M_PI*9.0*9.0*1000.0) > 4000.0 ) {
		ERROR("Intensity should be close to 1000*pi*peak_r^2\n");
		return 1;
	}
	if ( fabs(bg - sqrt(10.0)) > 1.0 ) {
		ERROR("Background should be close to sqrt(10)\n");
		return 1;
	}

	free(image.beam);
	free(image.det->panels);
	free(image.det);
	free(image.data);

	return 0;
}
