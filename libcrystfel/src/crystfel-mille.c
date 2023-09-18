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


/* Opposite of mille_label(), for decoding labels later */
enum gparam mille_unlabel(int n)
{
	switch ( n ) {
		case 1 : return GPARAM_DET_TX;
		case 2 : return GPARAM_DET_TY;
		case 3 : return GPARAM_DET_TZ;
		case 4 : return GPARAM_DET_RX;
		case 5 : return GPARAM_DET_RY;
		case 6 : return GPARAM_DET_RZ;
		default : abort();
	}
}


struct mille
{
	float *float_arr;
	int *int_arr;
	int max_entries;
	int n;
	FILE *fh;
};

typedef struct mille Mille;


static void mille_add_measurement(Mille *m,
                                  int NLC, float *derLc,
                                  int NGL, float *derGl, int *labels,
                                  float rMeas, float sigma)
{
	int space_required;
	int i;

	if ( m == NULL ) return;

	/* Allocate extra space if necessary */
	space_required = m->n + NLC + NGL + 2;
	if ( space_required > m->max_entries ) {

		float *new_float_arr;
		int *new_int_arr;
		int new_max_entries;

		if ( m->max_entries == 0 ) {
			new_max_entries = 256;
		} else {
			new_max_entries = m->max_entries;
		}

		while ( new_max_entries < space_required ) {
			new_max_entries *= 2;
		}

		new_float_arr = realloc(m->float_arr, new_max_entries*sizeof(float));
		new_int_arr = realloc(m->int_arr, new_max_entries*sizeof(int));
		if ( (new_float_arr == NULL) || (new_int_arr == NULL) ) return;

		m->float_arr = new_float_arr;
		m->int_arr = new_int_arr;
		m->max_entries = new_max_entries;
	}

	/* The measurement */
	m->float_arr[m->n] = rMeas;
	m->int_arr[m->n] = 0;
	m->n++;

	/* Local gradients */
	for ( i=0; i<NLC; i++ ) {
		if ( derLc[i] != 0.0 ) {
			m->float_arr[m->n] = derLc[i];
			m->int_arr[m->n] = i+1;
			m->n++;
		}
	}

	/* The measurement error */
	m->float_arr[m->n] = sigma;
	m->int_arr[m->n] = 0;
	m->n++;

	/* Global gradients */
	for ( i=0; i<NGL; i++ ) {
		if ( (derGl[i] != 0.0) && (labels[i] > 0) ) {
			m->float_arr[m->n] = derGl[i];
			m->int_arr[m->n] = labels[i];
			m->n++;
		}
	}
}


void write_mille(Mille *mille, int n, UnitCell *cell,
                 struct reflpeak *rps, struct image *image,
                 gsl_matrix **Minvs)
{
	int i;

	/* No groups -> no refinement */
	if ( image->detgeom->top_group == NULL ) return;

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
		float local_gradients_r[nl];
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
			local_gradients_r[j] = r_gradient(rvl[j], rps[i].refl,
			                                  cell, image->lambda);
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
		                      fs_dev(&rps[i], image->detgeom), 0.22);

		/* Add ss measurement */
		mille_add_measurement(mille,
		                      nl, local_gradients_ss,
		                      j, global_gradients_ss, labels,
		                      ss_dev(&rps[i], image->detgeom), 0.22);

		/* Add excitation error "measurement" (local-only) */
		mille_add_measurement(mille, nl, local_gradients_r,
		                      0, NULL, NULL, r_dev(&rps[i]), 1.0);
	}
}


Mille *crystfel_mille_new(const char *outFileName)
{
	Mille *m;

	m = malloc(sizeof(Mille));
	if ( m == NULL ) return NULL;

	m->max_entries = 0;
	m->n = 0;
	m->float_arr = NULL;
	m->int_arr = NULL;

	m->fh = fopen(outFileName, "wb");
	if ( m->fh == NULL ) {
		ERROR("Failed to open Mille file '%s'\n", outFileName);
		free(m);
		return NULL;
	}


	return m;
}


void crystfel_mille_free(Mille *m)
{
	if ( m == NULL ) return;
	fclose(m->fh);
	free(m->float_arr);
	free(m->int_arr);
	free(m);
}


void crystfel_mille_delete_last_record(Mille *m)
{
	m->n = 0;
}


void crystfel_mille_write_record(Mille *m)
{
	float nf = 0.0;
	int ni = 0;
	int nw = (m->n * 2)+2;

	fwrite(&nw, sizeof(int), 1, m->fh);

	fwrite(&nf, sizeof(float), 1, m->fh);
	fwrite(m->float_arr, sizeof(float), m->n, m->fh);

	fwrite(&ni, sizeof(int), 1, m->fh);
	fwrite(m->int_arr, sizeof(int), m->n, m->fh);
	m->n = 0;
}
