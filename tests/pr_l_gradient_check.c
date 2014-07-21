/*
 * pr_l_gradient_check.c
 *
 * Check Lorentz factor gradients for post refinement
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
#include <gsl/gsl_statistics.h>
#include <getopt.h>

#include <image.h>
#include <cell.h>
#include <cell-utils.h>
#include <geometry.h>
#include <reflist.h>
#include "../src/post-refinement.h"


static void scan_partialities(RefList *reflections, RefList *compare,
                              int *valid, long double *vals[3], int idx)
{
	int i;
	Reflection *refl;
	RefListIterator *iter;

	i = 0;
	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{
		signed int h, k, l;
		Reflection *refl2;

		get_indices(refl, &h, &k, &l);
		refl2 = find_refl(compare, h, k, l);
		if ( refl2 == NULL ) {
			valid[i] = 0;
			i++;
			continue;
		}

		vals[idx][i] = get_lorentz(refl2);
		i++;
	}
}


static void shift_parameter(struct image *image, int k, double shift)
{
	switch ( k )
	{
		case REF_DIV :
		image->div += shift;
		break;

		default:
		ERROR("This test can't check parameter %i\n", k);
		exit(1);
	}
}


static void calc_either_side(Crystal *cr, double incr_val,
                             int *valid, long double *vals[3], int refine,
                             PartialityModel pmodel)
{
	RefList *compare;
	struct image *image = crystal_get_image(cr);
	struct image im_moved;

	im_moved = *image;
	shift_parameter(&im_moved, refine, -incr_val);
	compare = find_intersections(&im_moved, cr, pmodel);
	scan_partialities(crystal_get_reflections(cr), compare,
	                  valid, vals, 0);
	reflist_free(compare);

	im_moved = *image;
	shift_parameter(&im_moved, refine, +incr_val);
	compare = find_intersections(&im_moved, cr, pmodel);
	scan_partialities(crystal_get_reflections(cr), compare,
	                  valid, vals, 2);
	reflist_free(compare);
}



static double test_gradients(Crystal *cr, double incr_val, int refine,
                             const char *str, const char *file,
                             PartialityModel pmodel, int quiet, int plot)
{
	Reflection *refl;
	RefListIterator *iter;
	long double *vals[3];
	int i;
	int *valid;
	int nref;
	int n_good, n_invalid, n_small, n_nan, n_bad;
	RefList *reflections;
	FILE *fh;
	int ntot = 0;
	double total = 0.0;
	char tmp[32];
	double *vec1;
	double *vec2;
	int n_line;
	double cc;

	reflections = find_intersections(crystal_get_image(cr), cr, pmodel);
	crystal_set_reflections(cr, reflections);

	nref = num_reflections(reflections);
	if ( nref < 10 ) {
		ERROR("Too few reflections found.  Failing test by default.\n");
		return 0.0;
	}

	vals[0] = malloc(nref*sizeof(long double));
	vals[1] = malloc(nref*sizeof(long double));
	vals[2] = malloc(nref*sizeof(long double));
	if ( (vals[0] == NULL) || (vals[1] == NULL) || (vals[2] == NULL) ) {
		ERROR("Couldn't allocate memory.\n");
		return 0.0;
	}

	valid = malloc(nref*sizeof(int));
	if ( valid == NULL ) {
		ERROR("Couldn't allocate memory.\n");
		return 0.0;
	}
	for ( i=0; i<nref; i++ ) valid[i] = 1;

	scan_partialities(reflections, reflections, valid, vals, 1);

	calc_either_side(cr, incr_val, valid, vals, refine, pmodel);

	if ( plot ) {
		snprintf(tmp, 32, "gradient-test-%s.dat", file);
		fh = fopen(tmp, "w");
	}

	vec1 = malloc(nref*sizeof(double));
	vec2 = malloc(nref*sizeof(double));
	if ( (vec1 == NULL) || (vec2 == NULL) ) {
		ERROR("Couldn't allocate memory.\n");
		return 0.0;
	}

	n_invalid = 0;  n_good = 0;
	n_nan = 0;  n_small = 0;  n_bad = 0;  n_line = 0;
	i = 0;
	for ( refl = first_refl(reflections, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) )
	{

		long double grad1, grad2, grad;
		double cgrad;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		if ( !valid[i] ) {
			n_invalid++;
			i++;
		} else {

			double r1, r2, p;
			int cl, ch;

			grad1 = (vals[1][i] - vals[0][i]) / incr_val;
			grad2 = (vals[2][i] - vals[1][i]) / incr_val;
			grad = (grad1 + grad2) / 2.0;
			i++;

			cgrad = l_gradient(cr, refine, refl, pmodel);

			get_partial(refl, &r1, &r2, &p, &cl, &ch);

			if ( isnan(cgrad) ) {
				n_nan++;
				continue;
			}

			if ( plot ) {
				fprintf(fh, "%e %Le\n", cgrad, grad);
			}

			vec1[n_line] = cgrad;
			vec2[n_line] = grad;
			n_line++;

			if ( (fabs(cgrad) < 5e-8) && (fabs(grad) < 5e-8) ) {
				n_small++;
				continue;
			}

			total += fabs(cgrad - grad);
			ntot++;

			if ( !within_tolerance(grad, cgrad, 5.0)
			  || !within_tolerance(cgrad, grad, 5.0) )
			{

				if ( !quiet ) {
					STATUS("!- %s %3i %3i %3i"
					       " %10.2Le %10.2e ratio = %5.2Lf"
					       " %10.2e %10.2e\n",
					       str, h, k, l, grad, cgrad,
					       cgrad/grad, r1, r2);
				}
				n_bad++;

			} else {

				//STATUS("OK %s %3i %3i %3i"
				//       " %10.2Le %10.2e ratio = %5.2Lf"
				//       " %10.2e %10.2e\n",
				//       str, h, k, l, grad, cgrad, cgrad/grad,
				//       r1, r2);

				n_good++;

			}

		}

	}

	STATUS("%3s: %3i within 5%%, %3i outside, %3i nan, %3i invalid, "
	       "%3i small. ", str, n_good, n_bad, n_nan, n_invalid, n_small);

	if ( plot ) {
		fclose(fh);
	}

	cc = gsl_stats_correlation(vec1, 1, vec2, 1, n_line);
	STATUS("CC = %+f\n", cc);
	return cc;
}


int main(int argc, char *argv[])
{
	struct image image;
	const double incr_frac = 1.0/1000000.0;
	double incr_val;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	UnitCell *cell;
	Crystal *cr;
	struct quaternion orientation;
	int i;
	int fail = 0;
	int quiet = 0;
	int plot = 0;
	int c;
	gsl_rng *rng;

	const struct option longopts[] = {
		{"quiet",       0, &quiet,        1},
		{"plot",        0, &plot,         1},
		{0, 0, NULL, 0}
	};

	while ((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (c) {

			case 0 :
			break;

			case '?' :
			break;

			default :
			ERROR("Unhandled option '%c'\n", c);
			break;

		}

	}

	image.width = 1024;
	image.height = 1024;
	image.det = simple_geometry(&image);
	image.det->panels[0].res = 13333.3;
	image.det->panels[0].clen = 80e-3;
	image.det->panels[0].coffset = 0.0;

	image.lambda = ph_en_to_lambda(eV_to_J(8000.0));
	image.div = 1e-3;
	image.bw = 0.01;
	image.filename = malloc(256);

	cr = crystal_new();
	if ( cr == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return 1;
	}
	crystal_set_mosaicity(cr, 0.0);
	crystal_set_profile_radius(cr, 0.005e9);
	crystal_set_image(cr, &image);

	cell = cell_new_from_parameters(10.0e-9, 10.0e-9, 10.0e-9,
	                                deg2rad(90.0),
	                                deg2rad(90.0),
	                                deg2rad(90.0));

	rng = gsl_rng_alloc(gsl_rng_mt19937);

	for ( i=0; i<3; i++ ) {

		UnitCell *rot;
		double val;
		PartialityModel pmodel;

		if ( i == 0 ) {
			pmodel = PMODEL_SPHERE;
			STATUS("Testing flat sphere model:\n");
		} else if ( i == 1 ) {
			pmodel = PMODEL_GAUSSIAN;
			STATUS("Testing Gaussian model:\n");
		} else if ( i == 2 ) {
			pmodel = PMODEL_THIN;
			STATUS("No need to test thin Ewald sphere model.\n");
			continue;
		} else {
			ERROR("WTF?\n");
			return 1;
		}

		orientation = random_quaternion(rng);
		rot = cell_rotate(cell, orientation);
		crystal_set_cell(cr, rot);

		cell_get_reciprocal(rot,
			            &ax, &ay, &az, &bx, &by,
			            &bz, &cx, &cy, &cz);

		incr_val = incr_frac * image.div;
		val =  test_gradients(cr, incr_val, REF_DIV, "div", "div",
		                      pmodel, quiet, plot);
		if ( val < 0.99 ) fail = 1;

	}

	gsl_rng_free(rng);

	return fail;
}
