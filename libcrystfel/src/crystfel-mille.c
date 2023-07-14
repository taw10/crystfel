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


int mille_label(int group_serial, enum gparam param)
{
	switch ( param ) {
		case GPARAM_DET_TX : return group_serial+1;  /* x-shift */
		case GPARAM_DET_TY : return group_serial+2;  /* y-shift */
		case GPARAM_DET_TZ : return group_serial+3;  /* z-shift */
		case GPARAM_DET_RX : return group_serial+4;  /* Rotation around x */
		case GPARAM_DET_RY : return group_serial+5;  /* Rotation around y */
		case GPARAM_DET_RZ : return group_serial+6;  /* Rotation around z */
		default : abort();
	}
}


void write_mille(Mille *mille, int n, UnitCell *cell,
                 struct reflpeak *rps, struct image *image,
                 gsl_matrix **Minvs)
{
	int i;

	/* Local parameters */
	const enum gparam rvl[] =
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
	const int nl = 9;

	/* Global parameters */
	const enum gparam rvg[] =
	{
		GPARAM_DET_TX,
		GPARAM_DET_TY,
		GPARAM_DET_TZ,
		GPARAM_DET_RX,
		GPARAM_DET_RY,
		GPARAM_DET_RZ,
	};
	const int ng = 6;
	const int max_hierarchy_levels = 8;

	for ( i=0; i<n; i++ ) {

		float local_gradients_fs[nl];
		float local_gradients_ss[nl];
		float global_gradients_fs[ng*max_hierarchy_levels];
		float global_gradients_ss[ng*max_hierarchy_levels];
		int labels[ng*max_hierarchy_levels];
		int j, levels;
		const struct detgeom_panel_group *group;

		/* Local gradients */
		for ( j=0; j<nl; j++ ) {
			fs_ss_gradient(rvl[j], rps[i].refl, cell,
			               &image->detgeom->panels[rps[i].peak->pn],
			               Minvs[rps[i].peak->pn], 0, 0, 0,
			               &local_gradients_fs[j],
			               &local_gradients_ss[j]);
		}

		/* Global gradients for each hierarchy level, starting at the
		 * individual panel and working up to the top level */
		j = 0;
		levels = 0;
		group = image->detgeom->panels[rps[i].peak->pn].group;
		while ( group != NULL ) {

			double cx, cy, cz;
			int g;

			detgeom_group_center(group, &cx, &cy, &cz);

			for ( g=0; g<ng; g++ ) {
				fs_ss_gradient(rvg[g], rps[i].refl, cell,
				               &image->detgeom->panels[rps[i].peak->pn],
				               Minvs[rps[i].peak->pn], 0, 0, 0,
				               &global_gradients_fs[j],
				               &global_gradients_ss[j]);
				labels[j] = mille_label(group->serial, rvg[g]);
				j++;
			}

			levels++;
			group = group->parent;

			if ( levels >= max_hierarchy_levels ) {
				ERROR("Too many nested hierarchy levels for refinement.\n");
				break;
			}
		}

		/* Add fs measurement */
		mille_add_measurement(mille,
		                      nl, local_gradients_fs,
		                      j, global_gradients_fs, labels,
		                      fs_dev(&rps[i], image->detgeom),
		                      0.65*image->detgeom->panels[rps[i].peak->pn].pixel_pitch);

		/* Add ss measurement */
		mille_add_measurement(mille,
		                      nl, local_gradients_ss,
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
