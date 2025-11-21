/*
 * image-raw.c
 *
 * Image loading from a raw data array
 *
 * Copyright © 2025 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2025 Thomas White <taw@physics.org>
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
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>

#include <image.h>
#include <utils.h>
#include <profile.h>

#include "datatemplate_priv.h"

/* Read the image data from 'data_block' into 'image', according to 'dtempl' */
int image_raw_read(struct image *image,
                   DataTemplate *dtempl,
                   void *data_block,
                   size_t data_block_size,
                   char *format)
{
	int w, h;
	int i;

	/* Must be simple 3D array, no placeholders, dims in right order,
	 * no fs/ss/pixel offsets */
	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( dtempl->panels[i].dims[0] < 0 ) return 1;
		if ( dtempl->panels[i].dims[1] != DIM_SS ) return 1;
		if ( dtempl->panels[i].dims[2] != DIM_FS ) return 1;
		if ( dtempl->panels[i].dims[3] != DIM_UNDEFINED ) return 1;
		if ( dtempl->panels[i].orig_min_fs != 0 ) return 1;
		if ( dtempl->panels[i].orig_min_ss != 0 ) return 1;
	}

	/* Format needs to be uint16 (for now) */
	if ( strcmp(format, "u16<") != 0 ) return 1;

	/* All panels must have same size */
	w = PANEL_WIDTH(&dtempl->panels[0]);
	h = PANEL_HEIGHT(&dtempl->panels[0]);
	for ( i=1; i<dtempl->n_panels; i++ ) {
		if ( PANEL_WIDTH(&dtempl->panels[i]) != w ) return 1;
		if ( PANEL_HEIGHT(&dtempl->panels[i]) != h ) return 1;
	}

	/* Total amount of data must be exactly as expected */
	if ( data_block_size != w * h * dtempl->n_panels * 2 ) return 1;

	uint16_t *in_data = data_block;
	for ( i=0; i<dtempl->n_panels; i++ ) {

		int fs, ss;

		for ( ss=0; ss<h; ss++ ) {
			for ( fs=0; fs<w; fs++ ) {
				image->dp[i][fs+ss*w] = in_data[fs+ss*w+i*w*h];
			}
		}
	}

	return 0;
}
