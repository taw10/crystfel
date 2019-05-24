/*
 * spectrum_check.c
 *
 * Check that Spectrum object works
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
#include <stdarg.h>

#include <spectrum.h>
#include <utils.h>

static int check_integral(Spectrum *s, int nsamp)
{
	double min, max, step;
	int i;
	double area = 0.0;

	spectrum_get_range(s, &min, &max);
	fprintf(stderr, "bounds: %e %e\n", min, max);
	step = (max-min)/nsamp;
	for ( i=0; i<=nsamp; i++ ) {
		double x = min+i*step;
		double y = spectrum_get_density_at_k(s, x);
		area += y * step;
	}
	fprintf(stderr, "Total area + %f\n", area);
	if ( area > 1.1 ) return 1;
	if ( area < 0.9 ) return 1;
	return 0;
}


static void plot_spectrum(Spectrum *s)
{
	double min, max, step;
	int i;
	const int nsamp = 100;

	spectrum_get_range(s, &min, &max);
	step = (max-min)/nsamp;
	for ( i=0; i<=nsamp; i++ ) {
		double x = min+i*step;
		double y = spectrum_get_density_at_k(s, x);
		printf("%e %e\n", x, y);
	}
}

int main(int argc, char *argv[])
{
	Spectrum *s;
	struct gaussian gauss;
	gsl_rng *rng;
	int r = 0;

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	s = spectrum_new();
	gauss.kcen = ph_eV_to_k(9000);
	gauss.sigma = ph_eV_to_k(100);
	gauss.height = 1.0;
	spectrum_set_gaussians(s, &gauss, 1);
	r += check_integral(s, 100);
	spectrum_free(s);

	s = spectrum_generate_sase(ph_eV_to_lambda(9000), 0.01, 0.0001, rng);
	r += check_integral(s, 100);
	plot_spectrum(s);
	spectrum_free(s);

	s = spectrum_generate_gaussian(ph_eV_to_lambda(9000), 0.01);
	r += check_integral(s, 100);
	spectrum_free(s);

	s = spectrum_generate_tophat(ph_eV_to_lambda(9000), 0.01);
	r += check_integral(s, 100);
	spectrum_free(s);

	s = spectrum_generate_twocolour(ph_eV_to_lambda(9000), 0.001, ph_eV_to_k(100));
	r += check_integral(s, 100);
	spectrum_free(s);

	gsl_rng_free(rng);

	return r;
}
