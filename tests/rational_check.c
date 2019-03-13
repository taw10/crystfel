/*
 * rational_check.c
 *
 * Check that rational numbers work
 *
 * Copyright © 2012-2019 Deutsches Elektronen-Synchrotron DESY,
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <gsl/gsl_rng.h>

#include <rational.h>
#include <utils.h>


static Rational gen_rtnl(gsl_rng *rng, double *dptr, int sz)
{
	Rational r;
	signed int num, den;
	do {
		num = gsl_rng_uniform_int(rng, 2*sz) - sz;
		den = gsl_rng_uniform_int(rng, 2*sz) - sz;
	} while ( den == 0 );
	r = rtnl(num, den);
	*dptr = (double)num/den;
	return r;
}


static int test_rational(gsl_rng *rng)
{
	Rational r1, r2, r3;
	double rd1, rd2, rd3;
	char op;

	r1 = gen_rtnl(rng, &rd1, 100);
	r2 = gen_rtnl(rng, &rd2, 100);

	switch ( gsl_rng_uniform_int(rng, 3) ) {

		case 0:
		r3 = rtnl_mul(r1, r2);
		rd3 = rd1*rd2;
		op = '*';
		break;

		case 1:
		r3 = rtnl_div(r1, r2);
		rd3 = rd1/rd2;
		op = '/';
		break;

		case 2:
		r3 = rtnl_add(r1, r2);
		rd3 = rd1+rd2;
		op = '+';
		break;

		case 3:
		r3 = rtnl_sub(r1, r2);
		rd3 = rd1-rd2;
		op = '-';
		break;

		default:
		abort();

	}

	if ( isinf(rd3) && isinf(rtnl_as_double(r3)) ) return 0;

	if ( fabs(rtnl_as_double(r3) - rd3) > 0.001 ) {
		ERROR("%10s %c %10s = %10s", rtnl_format(r1), op,
		      rtnl_format(r2), rtnl_format(r3));
		ERROR(" --->  %10f %c %10f = %10f\n", rd1, op, rd2, rd3);
		return 1;
	}

	return 0;
}


static int test_rational_matrix(gsl_rng *rng)
{
	const int size = 3;
	RationalMatrix *rm;
	Rational *rvec;
	Rational *rans;
	gsl_matrix *gm;
	gsl_vector *gvec;
	gsl_vector *gans;
	int i, j;
	int filt;
	int fail = 0;

	rm = rtnl_mtx_new(size, size);
	gm = gsl_matrix_alloc(size, size);
	rvec = malloc(size*sizeof(Rational));
	gvec = gsl_vector_alloc(size);
	rans = malloc(size*sizeof(Rational));
	if ( (rm==NULL) || (gm==NULL) || (rvec==NULL)
	  || (gvec==NULL) || (rans==NULL) ) return 1;

	/* Find a matrix equation which actually works */
	do {
		for ( i=0; i<size; i++ ) {
			double d;
			Rational r;
			for ( j=0; j<size; j++ ) {
				r = gen_rtnl(rng, &d, 5);
				rtnl_mtx_set(rm, i, j, r);
				gsl_matrix_set(gm ,i, j, d);
			}
			r = gen_rtnl(rng, &d, 5);
			rvec[i] = r;
			gsl_vector_set(gvec ,i, d);
		}

		gans = solve_svd(gvec, gm, &filt, 0);
	} while ( (gans==NULL) || (filt!=0) || (isnan(gsl_vector_get(gans,0))) );

	transform_fractional_coords_rtnl(rm, rvec, rans);

	for ( i=0; i<size; i++ ) {
		if ( fabs(rtnl_as_double(rans[i]) - gsl_vector_get(gans, i) > 0.001) )
		{
			STATUS("%25s %10f %10f\n", rtnl_format(rans[i]),
			                           rtnl_as_double(rans[i]),
			                           gsl_vector_get(gans, i));
			fail = 1;
		}
	}

	gsl_matrix_free(gm);
	rtnl_mtx_free(rm);
	free(rans);
	free(rvec);
	gsl_vector_free(gvec);
	gsl_vector_free(gans);
	return fail;
}


int main(int argc, char *argv[])
{
	int fail = 0;
	gsl_rng *rng;
	int i;

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_set_error_handler_off();

	STATUS("Arithmetic test...\n");
	for ( i=0; i<1000; i++ ) {
		fail += test_rational(rng);
	}

	STATUS("Matrix test...\n");
	for ( i=0; i<1000; i++ ) {
		fail += test_rational_matrix(rng);
	}

	gsl_rng_free(rng);

	if ( fail != 0 ) return 1;
	return 0;
}
