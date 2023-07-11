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


int mille_label(int hierarchy_level, int member_index, enum gparam param)
{
	int label;

	assert(member_index < 1000);

	label = 100000*hierarchy_level + 100*member_index;
	switch ( param ) {
		case GPARAM_DET_TX : return label+1;  /* x-shift */
		case GPARAM_DET_TY : return label+2;  /* y-shift */
		case GPARAM_DET_TZ : return label+3;  /* z-shift */
		case GPARAM_DET_RX : return label+4;  /* Rotation around x */
		case GPARAM_DET_RY : return label+5;  /* Rotation around y */
		case GPARAM_DET_RZ : return label+6;  /* Rotation around z */
		default : abort();
	}
}


void write_mille(Mille *mille, int n, UnitCell *cell,
                 struct reflpeak *rps, struct image *image,
                 gsl_matrix **Minvs)
{
	int i;

	for ( i=0; i<n; i++ ) {

		float local_gradients_fs[9];
		float local_gradients_ss[9];
		float global_gradients_fs[64];
		float global_gradients_ss[64];
		int labels[6];
		int j;
		const struct detgeom_panel_group *group;

		/* Local gradients */
		for ( j=0; j<9; j++ ) {
			fs_ss_gradient(rv[j], rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], 0, 0, 0,
			               &local_gradients_fs[j],
			               &local_gradients_ss[j]);
		}

		/* Global gradients (for each group, at least level of hierarchy) */
		j = 0;
		group = image->detgeom->panels[rps[i].peak->pn].group;
		while ( group != NULL ) {

			double cx, cy, cz;

			detgeom_group_center(group, &cx, &cy, &cz);

			fs_ss_gradient(GPARAM_DET_TX, rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], 0, 0, 0,
			               &global_gradients_fs[j],
			               &global_gradients_ss[j]);
			labels[j] = mille_label(group->hierarchy_level,
			                        group->member_index,
			                        GPARAM_DET_TX);
			j++;

			fs_ss_gradient(GPARAM_DET_TY, rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], 0, 0, 0,
			               &global_gradients_fs[j],
			               &global_gradients_ss[j]);
			labels[j] = mille_label(group->hierarchy_level,
			                        group->member_index,
			                        GPARAM_DET_TY);
			j++;

			fs_ss_gradient(GPARAM_DET_TZ, rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], 0, 0, 0,
			               &global_gradients_fs[j],
			               &global_gradients_ss[j]);
			labels[j] = mille_label(group->hierarchy_level,
			                        group->member_index,
			                        GPARAM_DET_TZ);
			j++;

			fs_ss_gradient(GPARAM_DET_RX, rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], cx, cy, cz,
			               &global_gradients_fs[j],
			               &global_gradients_ss[j]);
			labels[j] = mille_label(group->hierarchy_level,
			                        group->member_index,
			                        GPARAM_DET_RX);
			j++;

			fs_ss_gradient(GPARAM_DET_RY, rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], cx, cy, cz,
			               &global_gradients_fs[j],
			               &global_gradients_ss[j]);
			labels[j] = mille_label(group->hierarchy_level,
			                        group->member_index,
			                        GPARAM_DET_RX);
			j++;

			fs_ss_gradient(GPARAM_DET_RZ, rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], cx, cy, cz,
			               &global_gradients_fs[j],
			               &global_gradients_ss[j]);
			labels[j] = mille_label(group->hierarchy_level,
			                        group->member_index,
			                        GPARAM_DET_RZ);
			j++;

			group = group->parent;
		}

		/* Add fs measurement */
		mille_add_measurement(mille,
		                      9, local_gradients_fs,
		                      j, global_gradients_fs, labels,
		                      fs_dev(&rps[i], image->detgeom),
		                      0.65*image->detgeom->panels[rps[i].peak->pn].pixel_pitch);

		/* Add ss measurement */
		mille_add_measurement(mille,
		                      9, local_gradients_ss,
		                      j, global_gradients_ss, labels,
		                      ss_dev(&rps[i], image->detgeom),
		                      0.65*image->detgeom->panels[rps[i].peak->pn].pixel_pitch);
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
