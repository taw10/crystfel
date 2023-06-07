/*
 * crystfel-mille.c
 *
 * Interface to Millepede geometry refinement
 *
 * Copyright © 2023 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <assert.h>

#include "image.h"
#include "geometry.h"
#include "cell-utils.h"
#include "predict-refine.h"
#include "profile.h"

#include <mille_c_wrap.h>


static const enum gparam rv[] =
{
	GPARAM_ASX,
	GPARAM_ASY,
	GPARAM_ASZ,
	GPARAM_BSX,
	GPARAM_BSY,
	GPARAM_BSZ,
	GPARAM_CSX,
	GPARAM_CSY,
	GPARAM_CSZ,
};


int mille_label(int hierarchy_level, int member_index, char param)
{
	int label;

	assert(member_index < 1000);

	label = 100000*hierarchy_level + 100*member_index;
	switch ( param ) {
		case 'x' : return label+1;
		case 'y' : return label+2;
		case 'z' : return label+3;
		default : abort();
	}
}


void write_mille(Mille *mille, int n, UnitCell *cell,
                 struct reflpeak *rps, struct image *image,
                 double dx, double dy,
                 const struct detgeom_panel_group *group,
                 int hierarchy_level, int member_index)
{
	int i;
	float local_gradients[9];
	float global_gradients[6];
	int labels[6];

	/* Spot x-position terms */
	for ( i=0; i<n; i++ ) {

		int j;
		for ( j=0; j<9; j++ ) {
			local_gradients[j] = x_gradient(rv[j], rps[i].refl,
			                                cell, rps[i].panel);
		}

		global_gradients[0] = x_gradient(GPARAM_DETX, rps[i].refl,
		                                 cell, rps[i].panel);
		global_gradients[1] = x_gradient(GPARAM_DETY, rps[i].refl,
		                                 cell, rps[i].panel);
		global_gradients[2] = x_gradient(GPARAM_CLEN, rps[i].refl,
		                                 cell, rps[i].panel);
		labels[0] = mille_label(hierarchy_level, member_index, 'x');
		labels[1] = mille_label(hierarchy_level, member_index, 'y');
		labels[2] = mille_label(hierarchy_level, member_index, 'z');

		mille_add_measurement(mille,
		                      9, local_gradients,
		                      3, global_gradients, labels,
		                      x_dev(&rps[i], image->detgeom, dx, dy),
		                      0.65*rps[i].panel->pixel_pitch);
	}

	/* Spot y-position terms */
	for ( i=0; i<n; i++ ) {

		int j;
		for ( j=0; j<9; j++ ) {
			local_gradients[j] = y_gradient(rv[j], rps[i].refl,
			                                cell, rps[i].panel);
		}

		global_gradients[0] = y_gradient(GPARAM_DETX, rps[i].refl,
		                                 cell, rps[i].panel);
		global_gradients[1] = y_gradient(GPARAM_DETY, rps[i].refl,
		                                 cell, rps[i].panel);
		global_gradients[2] = y_gradient(GPARAM_CLEN, rps[i].refl,
		                                 cell, rps[i].panel);
		labels[0] = mille_label(hierarchy_level, member_index, 'x');
		labels[1] = mille_label(hierarchy_level, member_index, 'y');
		labels[2] = mille_label(hierarchy_level, member_index, 'z');

		mille_add_measurement(mille,
		                      9, local_gradients,
		                      3, global_gradients, labels,
		                      y_dev(&rps[i], image->detgeom, dx, dy),
		                      0.65*rps[i].panel->pixel_pitch);
	}

	/* Next level of hierarchy */
	for ( i=0; i<group->n_children; i++ ) {
		write_mille(mille, n, cell, rps, image, dx, dy,
		            group->children[i], hierarchy_level+1, i);
	}

}


Mille *crystfel_mille_new(const char *outFileName,
                          int asBinary,
                          int writeZero)
{
	return mille_new(outFileName, asBinary, writeZero);
}


void crystfel_mille_free(Mille *m)
{
	mille_free(m);
}


void crystfel_mille_delete_last_record(Mille *m)
{
	mille_delete_last_record(m);
}


void crystfel_mille_write_record(Mille *m)
{
	mille_write_record(m);
}
