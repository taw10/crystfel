/*
 * prediction_gradient_check.c
 *
 * Check partiality gradients for prediction refinement
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
#include <gsl/gsl_statistics.h>
#include <getopt.h>

#include <image.h>
#include <cell.h>
#include <cell-utils.h>
#include <geometry.h>
#include <reflist.h>
#include "../src/predict-refine.c"


int checkrxy;


static void scan(RefList *reflections, RefList *compare,
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
		double rlow, rhigh, p;
		double fs, ss, xh, yh;
		struct panel *panel;

		get_indices(refl, &h, &k, &l);
		refl2 = find_refl(compare, h, k, l);
		if ( refl2 == NULL ) {
			valid[i] = 0;
			i++;
			continue;
		}

		get_partial(refl2, &rlow, &rhigh, &p);
		get_detector_pos(refl2, &fs, &ss);
		panel = get_panel(refl2);
		twod_mapping(fs, ss, &xh, &yh, panel);

		switch ( checkrxy ) {

			case 0 :
			vals[idx][i] = (rlow + rhigh)/2.0;
			break;

			case 1 :
			vals[idx][i] = xh;
			break;

			case 2 :
			vals[idx][i] = yh;
			break;
		}

		i++;
	}
}


static UnitCell *new_shifted_cell(UnitCell *input, int k, double shift)
{
	UnitCell *cell;
	double asx, asy, asz;
	double bsx, bsy, bsz;
	double csx, csy, csz;

	cell = cell_new();
	cell_get_reciprocal(input, &asx, &asy, &asz, &bsx, &bsy, &bsz,
	                    &csx, &csy, &csz);
	switch ( k )
	{
		case GPARAM_ASX :  asx += shift;  break;
		case GPARAM_ASY :  asy += shift;  break;
		case GPARAM_ASZ :  asz += shift;  break;
		case GPARAM_BSX :  bsx += shift;  break;
		case GPARAM_BSY :  bsy += shift;  break;
		case GPARAM_BSZ :  bsz += shift;  break;
		case GPARAM_CSX :  csx += shift;  break;
		case GPARAM_CSY :  csy += shift;  break;
		case GPARAM_CSZ :  csz += shift;  break;
	}
	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	return cell;
}


static Crystal *new_shifted_crystal(Crystal *cr, int refine, double incr_val)
{
	Crystal *cr_new;
	UnitCell *cell;

	cr_new = crystal_copy(cr);
	if ( cr_new == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return NULL;
	}

	crystal_set_image(cr_new, crystal_get_image(cr));

	switch ( refine ) {

		case GPARAM_ASX :
		case GPARAM_ASY :
		case GPARAM_ASZ :
		case GPARAM_BSX :
		case GPARAM_BSY :
		case GPARAM_BSZ :
		case GPARAM_CSX :
		case GPARAM_CSY :
		case GPARAM_CSZ :
		cell = new_shifted_cell(crystal_get_cell(cr), refine,
		                        incr_val);
		crystal_set_cell(cr_new, cell);
		break;

		default:
		ERROR("Can't shift %i\n", refine);
		break;

	}

	return cr_new;
}


static void calc_either_side(Crystal *cr, double incr_val,
                             int *valid, long double *vals[3], int refine)
{
	RefList *compare;
	struct image *image = crystal_get_image(cr);
	Crystal *cr_new;

	cr_new = new_shifted_crystal(cr, refine, -incr_val);
	compare = find_intersections(image, cr_new, PMODEL_SCSPHERE);
	scan(crystal_get_reflections(cr), compare, valid, vals, 0);
	cell_free(crystal_get_cell(cr_new));
	crystal_free(cr_new);
	reflist_free(compare);

	cr_new = new_shifted_crystal(cr, refine, +incr_val);
	compare = find_intersections(image, cr_new, PMODEL_SCSPHERE);
	scan(crystal_get_reflections(cr), compare, valid, vals, 2);
	cell_free(crystal_get_cell(cr_new));
	crystal_free(cr_new);
	reflist_free(compare);
}


static double test_gradients(Crystal *cr, double incr_val, int refine,
                             const char *str, const char *file,
                             int quiet, int plot)
{
	Reflection *refl;
	RefListIterator *iter;
	long double *vals[3];
	int i;
	int *valid;
	int nref;
	int n_good, n_invalid, n_small, n_nan, n_bad;
	RefList *reflections;
	FILE *fh = NULL;
	int ntot = 0;
	double total = 0.0;
	char tmp[32];
	double *vec1;
	double *vec2;
	int n_line;
	double cc;

	reflections = find_intersections(crystal_get_image(cr), cr,
	                                 PMODEL_SCSPHERE);
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

	scan(reflections, reflections, valid, vals, 1);

	calc_either_side(cr, incr_val, valid, vals, refine);

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
		struct reflpeak rp;

		get_indices(refl, &h, &k, &l);

		if ( !valid[i] ) {
			n_invalid++;
			i++;
		} else {

			double r1, r2, p;

			grad1 = (vals[1][i] - vals[0][i]) / incr_val;
			grad2 = (vals[2][i] - vals[1][i]) / incr_val;
			grad = (grad1 + grad2) / 2.0;
			i++;

			if ( checkrxy == 0 ) {

				cgrad = r_gradient(crystal_get_cell(cr), refine,
				                   refl, crystal_get_image(cr));

			} else {

				struct imagefeature pk;
				struct image *image;

				pk.fs = 0.0;
				pk.ss = 0.0;

				image = crystal_get_image(cr);
				rp.refl = refl;
				rp.peak = &pk;
				rp.panel = &image->det->panels[0];

				if ( checkrxy == 1 ) {
					cgrad = x_gradient(refine, &rp,
					          crystal_get_image(cr)->det,
					          crystal_get_image(cr)->lambda,
					          crystal_get_cell(cr));
				} else {
					cgrad = y_gradient(refine, &rp,
					          crystal_get_image(cr)->det,
					          crystal_get_image(cr)->lambda,
					          crystal_get_cell(cr));
				}
			}

			get_partial(refl, &r1, &r2, &p);

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

			if ( (fabs(cgrad) < 5e-12) && (fabs(grad) < 5e-12) ) {
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
	const double incr_frac = 1.0/100000.0;
	double incr_val;
	double ax, ay, az;
	double bx, by, bz;
	double cx, cy, cz;
	UnitCell *cell;
	Crystal *cr;
	struct quaternion orientation;
	int fail = 0;
	int quiet = 0;
	int plot = 0;
	int c;
	gsl_rng *rng;
	UnitCell *rot;
	double val;

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

	for ( checkrxy=0; checkrxy<3; checkrxy++ ) {


		switch ( checkrxy ) {
			case 0 :
			STATUS("Excitation error:\n");
			break;
			case 1:
			STATUS("x coordinate:\n");
			break;
			default:
			case 2:
			STATUS("y coordinate:\n");
			break;
			STATUS("WTF??\n");
			break;
		}

		orientation = random_quaternion(rng);
		rot = cell_rotate(cell, orientation);
		crystal_set_cell(cr, rot);

		cell_get_reciprocal(rot, &ax, &ay, &az,
		                    &bx, &by, &bz, &cx, &cy, &cz);

		if ( checkrxy != 2 ) {

			incr_val = incr_frac * ax;
			val = test_gradients(cr, incr_val, GPARAM_ASX,
			                     "ax*", "ax", quiet, plot);
			if ( val < 0.99 ) fail = 1;
			incr_val = incr_frac * bx;
			val = test_gradients(cr, incr_val, GPARAM_BSX,
			                     "bx*", "bx", quiet, plot);
			if ( val < 0.99 ) fail = 1;
			incr_val = incr_frac * cx;
			val = test_gradients(cr, incr_val, GPARAM_CSX,
			                     "cx*", "cx", quiet, plot);
			if ( val < 0.99 ) fail = 1;

		}

		if ( checkrxy != 1 ) {

			incr_val = incr_frac * ay;
			val = test_gradients(cr, incr_val, GPARAM_ASY,
			                     "ay*", "ay", quiet, plot);
			if ( val < 0.99 ) fail = 1;
			incr_val = incr_frac * by;
			val = test_gradients(cr, incr_val, GPARAM_BSY,
			                     "by*", "by", quiet, plot);
			if ( val < 0.99 ) fail = 1;
			incr_val = incr_frac * cy;
			val = test_gradients(cr, incr_val, GPARAM_CSY,
			                     "cy*", "cy", quiet, plot);
			if ( val < 0.99 ) fail = 1;

		}

		incr_val = incr_frac * az;
		val = test_gradients(cr, incr_val, GPARAM_ASZ, "az*", "az",
		                     quiet, plot);
		if ( val < 0.99 ) fail = 1;
		incr_val = incr_frac * bz;
		val = test_gradients(cr, incr_val, GPARAM_BSZ, "bz*", "bz",
		                     quiet, plot);
		if ( val < 0.99 ) fail = 1;
		incr_val = incr_frac * cz;
		val = test_gradients(cr, incr_val, GPARAM_CSZ, "cz*", "cz",
		                     quiet, plot);
		if ( val < 0.99 ) fail = 1;

	}

	gsl_rng_free(rng);

	return fail;
}
