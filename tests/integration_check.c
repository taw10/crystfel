/*
 * integration_check.c
 *
 * Check peak integration
 *
 * (c) 2011 Thomas White <taw@physics.org>
 * (c) 2011 Andrew Martin <andrew.martin@desy.de>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>

#include <image.h>
#include <peaks.h>
#include <utils.h>
#include <beam-parameters.h>


/* The third integration check draws a Poisson background and checks that, on
 * average, it gets subtracted by the background subtraction. */
static void third_integration_check(struct image *image, int n_trials,
                                    int *fail)
{
	double mean_intensity = 0.0;
	double mean_bg = 0.0;
	double mean_max = 0.0;
	double mean_sigma = 0.0;
	int i;
	int fs, ss;

	for ( i=0; i<n_trials; i++ ) {

		double intensity, bg, max, sigma;
		double fsp, ssp;

		for ( fs=0; fs<image->width; fs++ ) {
		for ( ss=0; ss<image->height; ss++ ) {
			image->data[fs+image->width*ss] = poisson_noise(10.0);
		}
		}
		integrate_peak(image, 64, 64, &fsp, &ssp, &intensity,
			       &bg, &max, &sigma, 0, 1, 1);

		mean_intensity += intensity;
		mean_bg += bg;
		mean_max += max;
		mean_sigma += sigma;

	}
	mean_intensity /= n_trials;
	mean_bg /= n_trials;
	mean_max /= n_trials;
	mean_sigma /= n_trials;

	STATUS("  Third check (mean values): intensity = %.2f, bg = %.2f,"
	       " max = %.2f, sigma = %.2f\n",
	       mean_intensity, mean_bg, mean_max, mean_sigma);
	if ( fabs(mean_intensity) > 5.0 ) {
		ERROR("Mean intensity should be close to zero.\n");
		*fail = 1;
	}
	if ( fabs(mean_bg-10.0) > 0.3 ) {
		ERROR("Mean background should be close to ten.\n");
		*fail = 1;
	}
	if ( fabs(mean_intensity) > mean_sigma ) {
		ERROR("Mean intensity should be less than mean sigma.\n");
		*fail = 1;
	}
}


/* The fourth integration check draws a Poisson background and draws a peak on
 * top of it, then checks that the intensity of the peak is correctly recovered
 * accounting for the background. */
static void fourth_integration_check(struct image *image, int n_trials,
                                     int *fail)
{
	double mean_intensity = 0.0;
	double mean_bg = 0.0;
	double mean_max = 0.0;
	double mean_sigma = 0.0;
	int i;
	int fs, ss;
	int pcount = 0;

	for ( i=0; i<n_trials; i++ ) {

		double intensity, bg, max, sigma;
		double fsp, ssp;

		for ( fs=0; fs<image->width; fs++ ) {
		for ( ss=0; ss<image->height; ss++ ) {
			int idx = fs+image->width*ss;
			image->data[idx] = 0.0;
			if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) > 9*9 ) {
				image->data[idx] = poisson_noise(10.0);
			} else {
				image->data[idx] += 1000.0;
				pcount++;
			}
		}
		}
		integrate_peak(image, 64, 64, &fsp, &ssp, &intensity,
			       &bg, &max, &sigma, 0, 1, 1);

		mean_intensity += intensity;
		mean_bg += bg;
		mean_max += max;
		mean_sigma += sigma;

	}
	mean_intensity /= n_trials;
	mean_bg /= n_trials;
	mean_max /= n_trials;
	mean_sigma /= n_trials;
	pcount /= n_trials;

	STATUS(" Fourth check (mean values): intensity = %.2f, bg = %.2f,"
	       " max = %.2f, sigma = %.2f\n",
	       mean_intensity, mean_bg, mean_max, mean_sigma);
	if ( fabs(mean_intensity - pcount*1000.0) > 4000.0 ) {
		ERROR("Mean intensity should be close to %f\n", pcount*1000.0);
		*fail = 1;
	}
	if ( fabs(mean_bg-10.0) > 0.3 ) {
		ERROR("Mean background should be close to ten.\n");
		*fail = 1;
	}
	if ( fabs(mean_intensity) < mean_sigma ) {
		ERROR("Mean intensity should be greater than mean sigma.\n");
		*fail = 1;
	}
}


/* The fifth check integrates a Poisson background with background subtraction
 * switched off, and checks that the result is what would be expected. */
static void fifth_integration_check(struct image *image, int n_trials,
                                    int *fail)
{
	double mean_intensity = 0.0;
	double mean_bg = 0.0;
	double mean_max = 0.0;
	double mean_sigma = 0.0;
	int i;
	int fs, ss;
	int pcount = 0;

	for ( i=0; i<n_trials; i++ ) {

		double intensity, bg, max, sigma;
		double fsp, ssp;

		for ( fs=0; fs<image->width; fs++ ) {
		for ( ss=0; ss<image->height; ss++ ) {
			int idx = fs+image->width*ss;
			image->data[idx] = poisson_noise(10.0);
			if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) < 10*10 ) {
				pcount++;
			}
		}
		}
		integrate_peak(image, 64, 64, &fsp, &ssp, &intensity,
			       &bg, &max, &sigma, 0, 1, 0);

		mean_intensity += intensity;
		mean_bg += bg;
		mean_max += max;
		mean_sigma += sigma;

	}
	mean_intensity /= n_trials;
	mean_bg /= n_trials;
	mean_max /= n_trials;
	mean_sigma /= n_trials;
	pcount /= n_trials;

	STATUS("  Fifth check (mean values): intensity = %.2f, bg = %.2f,"
	       " max = %.2f, sigma = %.2f\n",
	       mean_intensity, mean_bg, mean_max, mean_sigma);
	double s = pcount*10.0;
	if ( fabs(mean_intensity - s) > 5.0 ) {
		ERROR("Mean intensity should be close to %f.\n", pcount*10.0);
		*fail = 1;
	}
	if ( fabs(mean_bg-10.0) > 0.3 ) {
		ERROR("Mean background should be close to ten.\n");
		*fail = 1;
	}
}


/* The sixth check is like the fourth check, except that the background
 * subtraction is switched off */
static void sixth_integration_check(struct image *image, int n_trials,
                                    int *fail)
{
	double mean_intensity = 0.0;
	double mean_bg = 0.0;
	double mean_max = 0.0;
	double mean_sigma = 0.0;
	int i;
	int fs, ss;
	int pcount = 0;
	int npcount = 0;

	for ( i=0; i<n_trials; i++ ) {

		double intensity, bg, max, sigma;
		double fsp, ssp;

		for ( fs=0; fs<image->width; fs++ ) {
		for ( ss=0; ss<image->height; ss++ ) {
			int idx = fs+image->width*ss;
			double r = (fs-64)*(fs-64) + (ss-64)*(ss-64);
			image->data[idx] = poisson_noise(10.0);
			if ( r < 9*9 ) {
				image->data[idx] += 1000.0;
				pcount++;
			} else if ( r < 10*10 ) {
				npcount++;
			}
		}
		}
		integrate_peak(image, 64, 64, &fsp, &ssp, &intensity,
			       &bg, &max, &sigma, 0, 1, 0);

		mean_intensity += intensity;
		mean_bg += bg;
		mean_max += max;
		mean_sigma += sigma;

	}
	mean_intensity /= n_trials;
	mean_bg /= n_trials;
	mean_max /= n_trials;
	mean_sigma /= n_trials;
	pcount /= n_trials;
	npcount /= n_trials;

	STATUS("  Sixth check (mean values): intensity = %.2f, bg = %.2f,"
	       " max = %.2f, sigma = %.2f\n",
	       mean_intensity, mean_bg, mean_max, mean_sigma);

	double s = pcount*1010.0 + npcount*10.0;
	if ( fabs(mean_intensity - s) > 4000.0 ) {
		ERROR("Mean intensity should be close to %f.\n", s);
		*fail = 1;
	}
	if ( fabs(mean_bg-10.0) > 0.3 ) {
		ERROR("Mean background should be close to ten.\n");
		*fail = 1;
	}

	STATUS("            (Absolute value of mean sigma not (yet) tested)\n");
}


int main(int argc, char *argv[])
{
	struct image image;
	double fsp, ssp, intensity, bg, max, sigma;
	int fs, ss;
	FILE *fh;
	unsigned int seed;
	int fail = 0;
	const int n_trials = 1000;

	fh = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fh);
	fclose(fh);
	srand(seed);

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
	integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity,
	               &bg, &max, &sigma, 0, 1, 1);
	STATUS("  First check: intensity = %.2f, bg = %.2f, max = %.2f,"
	       " sigma = %.2f\n", intensity, bg, max, sigma);
	if ( intensity != 0.0 ) {
		ERROR("Intensity should be zero.\n");
		fail = 1;
	}
	if ( bg != 0.0 ) {
		ERROR("Background should be zero.\n");
		fail = 1;
	}

	/* Second check: uniform peak gives correct value and no bg */
	for ( fs=0; fs<image.width; fs++ ) {
	for ( ss=0; ss<image.height; ss++ ) {
		if ( (fs-64)*(fs-64) + (ss-64)*(ss-64) > 9*9 ) continue;
		image.data[fs+image.width*ss] = 1000.0;
	}
	}
	integrate_peak(&image, 64, 64, &fsp, &ssp, &intensity,
	               &bg, &max, &sigma, 0, 1, 1);
	STATUS(" Second check: intensity = %.2f, bg = %.2f, max = %.2f,"
	       " sigma = %.2f\n", intensity, bg, max, sigma);
	if ( fabs(intensity - M_PI*9.0*9.0*1000.0) > 4000.0 ) {
		ERROR("Intensity should be close to 1000*pi*integr_r^2\n");
		fail = 1;
	}
	if ( bg != 0.0 ) {
		ERROR("Background should be zero.\n");
		fail = 1;
	}
	if ( max != 1000.0 ) {
		ERROR("Max should be 1000.\n");
		fail = 1;
	}
	if ( sigma != 0.0 ) {
		ERROR("Sigma should be zero.\n");
		fail = 1;
	}

	/* Third check: Poisson background should get mostly subtracted */
	third_integration_check(&image, n_trials, &fail);

	/* Fourth check: peak on Poisson background */
	fourth_integration_check(&image, n_trials, &fail);

	/* Fifth check: Like third check but without background subtraction */
	fifth_integration_check(&image, n_trials, &fail);

	/* Sixth check: Like fourth check but without background subtraction */
	sixth_integration_check(&image, n_trials, &fail);

	free(image.beam);
	free(image.det->panels);
	free(image.det);
	free(image.data);

	if ( fail ) return 1;

	return 0;
}
