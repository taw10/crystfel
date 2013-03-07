/*
 * pr_gradient_check.c
 *
 * Check gradients for post refinement
 *
 * Copyright © 2012-2013 Thomas White <taw@physics.org>
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
		double r1, r2, p;
		int clamp_low, clamp_high;

		get_indices(refl, &h, &k, &l);
		refl2 = find_refl(compare, h, k, l);
		if ( refl2 == NULL ) {
			valid[i] = 0;
			i++;
			continue;
		}

		get_partial(refl2, &r1, &r2, &p, &clamp_low, &clamp_high);
		if ( clamp_low && clamp_high ) {
			if ( !within_tolerance(p, 1.0, 0.001) ) {

				signed int h, k, l;

				get_indices(refl, &h, &k, &l);

				ERROR("%3i %3i %3i - double clamped but"
				      " partiality not close to 1.0 (%5.2f)\n",
				      h, k, l, p);

			}
			valid[i] = 0;
		}

		vals[idx][i] = p;
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
		case REF_ASX :  asx += shift;  break;
		case REF_ASY :  asy += shift;  break;
		case REF_ASZ :  asz += shift;  break;
		case REF_BSX :  bsx += shift;  break;
		case REF_BSY :  bsy += shift;  break;
		case REF_BSZ :  bsz += shift;  break;
		case REF_CSX :  csx += shift;  break;
		case REF_CSY :  csy += shift;  break;
		case REF_CSZ :  csz += shift;  break;
	}
	cell_set_reciprocal(cell, asx, asy, asz, bsx, bsy, bsz, csx, csy, csz);

	return cell;
}


static void shift_parameter(struct image *image, int k, double shift)
{
	switch ( k )
	{
		case REF_DIV : image->div += shift;  break;
	}
}


static Crystal *new_shifted_crystal(Crystal *cr, int refine, double incr_val)
{
	Crystal *cr_new;
	double r;
	UnitCell *cell;

	cr_new = crystal_copy(cr);
	if ( cr_new == NULL ) {
		ERROR("Failed to allocate crystal.\n");
		return NULL;
	}

	crystal_set_image(cr_new, crystal_get_image(cr));
	r = crystal_get_profile_radius(cr_new);

	switch ( refine ) {

		case REF_ASX :
		case REF_ASY :
		case REF_ASZ :
		case REF_BSX :
		case REF_BSY :
		case REF_BSZ :
		case REF_CSX :
		case REF_CSY :
		case REF_CSZ :
		cell = new_shifted_cell(crystal_get_cell(cr), refine,
		                        incr_val);
		crystal_set_cell(cr_new, cell);
		break;

		case REF_R :
		cell = cell_new_from_cell(crystal_get_cell(cr));
		crystal_set_cell(cr_new, cell);
		crystal_set_profile_radius(cr_new, r + incr_val);
		break;

		default :
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

	if ( (refine != REF_DIV) ) {

		Crystal *cr_new;

		/* Crystal properties */
		cr_new = new_shifted_crystal(cr, refine, -incr_val);
		compare = find_intersections(image, cr_new);
		scan_partialities(crystal_get_reflections(cr), compare, valid,
		                  vals, 0);
		cell_free(crystal_get_cell(cr_new));
		crystal_free(cr_new);
		reflist_free(compare);

		cr_new = new_shifted_crystal(cr, refine, +incr_val);
		compare = find_intersections(image, cr_new);
		scan_partialities(crystal_get_reflections(cr), compare, valid,
		                  vals, 2);
		cell_free(crystal_get_cell(cr_new));
		crystal_free(cr_new);
		reflist_free(compare);

	} else {

		struct image im_moved;

		/* "Image" properties */
		im_moved = *image;
		shift_parameter(&im_moved, refine, -incr_val);
		compare = find_intersections(&im_moved, cr);
		scan_partialities(crystal_get_reflections(cr), compare,
		                  valid, vals, 0);
		reflist_free(compare);

		im_moved = *image;
		shift_parameter(&im_moved, refine, +incr_val);
		compare = find_intersections(&im_moved, cr);
		scan_partialities(crystal_get_reflections(cr), compare,
		                  valid, vals, 2);
		reflist_free(compare);

	}
}



static double test_gradients(Crystal *cr, double incr_val, int refine,
                             const char *str, PartialityModel pmodel)
{
	Reflection *refl;
	RefListIterator *iter;
	long double *vals[3];
	int i;
	int *valid;
	int nref;
	int n_good, n_invalid, n_small, n_nan, n_bad;
	RefList *reflections;
	//FILE *fh;
	int ntot = 0;
	double total = 0.0;

	reflections = find_intersections(crystal_get_image(cr), cr);
	crystal_set_reflections(cr, reflections);

	nref = num_reflections(reflections);
	if ( nref < 10 ) {
		ERROR("Too few reflections found.  Failing test by default.\n");
		return -1;
	}

	vals[0] = malloc(nref*sizeof(long double));
	vals[1] = malloc(nref*sizeof(long double));
	vals[2] = malloc(nref*sizeof(long double));
	if ( (vals[0] == NULL) || (vals[1] == NULL) || (vals[2] == NULL) ) {
		ERROR("Couldn't allocate memory.\n");
		return -1;
	}

	valid = malloc(nref*sizeof(int));
	if ( valid == NULL ) {
		ERROR("Couldn't allocate memory.\n");
		return -1;
	}
	for ( i=0; i<nref; i++ ) valid[i] = 1;

	scan_partialities(reflections, reflections, valid, vals, 1);

	calc_either_side(cr, incr_val, valid, vals, refine);

	//fh = fopen("wrongness.dat", "a");

	n_invalid = 0;  n_good = 0;
	n_nan = 0;  n_small = 0;  n_bad = 0;
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
		} else {

			double r1, r2, p;
			int cl, ch;

			grad1 = (vals[1][i] - vals[0][i]) / incr_val;
			grad2 = (vals[2][i] - vals[1][i]) / incr_val;
			grad = (grad1 + grad2) / 2.0;

			cgrad = gradient(cr, refine, refl, pmodel);

			get_partial(refl, &r1, &r2, &p, &cl, &ch);

			if ( fabs(cgrad) < 5e-8 ) {
				n_small++;
				continue;
			}

			if ( isnan(cgrad) ) {
				n_nan++;
				continue;
			}

			total += fabs(cgrad - grad);
			ntot++;

			if ( !within_tolerance(grad, cgrad, 5.0) )
			{

				STATUS("!- %s %3i %3i %3i"
				       " %10.2Le %10.2e ratio = %5.2Lf"
				       " %10.2e %10.2e\n",
				       str, h, k, l, grad, cgrad, cgrad/grad,
				       r1, r2);
				n_bad++;

			} else {

				//STATUS("OK %s %3i %3i %3i"
				//       " %10.2Le %10.2e ratio = %5.2Lf"
				//       " %10.2e %10.2e\n",
				//       str, h, k, l, grad, cgrad, cgrad/grad,
				//       r1, r2);

				n_good++;

			}

			//fprintf(fh, "%e %f\n",
			        //resolution(image->indexed_cell, h, k, l),
			        //rad2deg(tt),
			//        cgrad,
			//        fabs((grad-cgrad)/grad));

		}

		i++;

	}

	STATUS("%3s: %3i within 5%%, %3i outside, %3i nan, %3i invalid, "
	       "%3i small. ", str, n_good, n_bad, n_nan, n_invalid, n_small);
	//fclose(fh);

	STATUS("Fractional error = %f %%\n", 100.0*total/ntot);
	return total / ntot;
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
	const PartialityModel pmodel = PMODEL_SPHERE;
	double tot = 0.0;

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

	for ( i=0; i<1; i++ ) {

		UnitCell *rot;

		orientation = random_quaternion();
		rot = cell_rotate(cell, orientation);
		crystal_set_cell(cr, rot);

		cell_get_reciprocal(rot,
			            &ax, &ay, &az, &bx, &by,
			            &bz, &cx, &cy, &cz);

		incr_val = incr_frac * image.div;
		tot += test_gradients(cr, incr_val, REF_DIV, "div", pmodel);

		incr_val = incr_frac * crystal_get_profile_radius(cr);
		tot += test_gradients(cr, incr_val, REF_R, "R", pmodel);

		incr_val = incr_frac * ax;
		tot += test_gradients(cr, incr_val, REF_ASX, "ax*", pmodel);
		incr_val = incr_frac * ay;
		tot += test_gradients(cr, incr_val, REF_ASY, "ay*", pmodel);
		incr_val = incr_frac * az;
		tot += test_gradients(cr, incr_val, REF_ASZ, "az*", pmodel);

		incr_val = incr_frac * bx;
		tot += test_gradients(cr, incr_val, REF_BSX, "bx*", pmodel);
		incr_val = incr_frac * by;
		tot += test_gradients(cr, incr_val, REF_BSY, "by*", pmodel);
		incr_val = incr_frac * bz;
		tot += test_gradients(cr, incr_val, REF_BSZ, "bz*", pmodel);

		incr_val = incr_frac * cx;
		tot += test_gradients(cr, incr_val, REF_CSX, "cx*", pmodel);
		incr_val = incr_frac * cy;
		tot += test_gradients(cr, incr_val, REF_CSY, "cy*", pmodel);
		incr_val = incr_frac * cz;
		tot += test_gradients(cr, incr_val, REF_CSZ, "cz*", pmodel);

	}

	STATUS("Mean fractional error = %f %%\n", 100.0*tot/10.0);
	return 0;
}
