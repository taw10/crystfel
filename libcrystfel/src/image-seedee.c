/*
 * image-seedee.c
 *
 * Image loading with Seedee
 *
 * Copyright © 2017-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018-2022 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *   2017      Stijn de Graaf
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


#if defined(HAVE_SEEDEE)

#include <seedee/seedee.h>
#include <cjson/cJSON.h>


static int load_seedee_data(struct panel_template *p,
                            struct SeedeeNDArray *array,
                            float *data, int *bad)
{
	int data_size_fs, data_size_ss;

	data_size_ss = array->shape[0];
	data_size_fs = array->shape[1];

	if ( (p->orig_min_fs + PANEL_WIDTH(p) > data_size_fs)
	  || (p->orig_min_ss + PANEL_HEIGHT(p) > data_size_ss) )
	{
		ERROR("Data for panel %s (%i x %i + %i + %i) is outside data "
		      "array bounds (%i x %i)\n",
		      p->name,
		      PANEL_WIDTH(p), PANEL_HEIGHT(p),
		      p->orig_min_fs, p->orig_min_ss,
		      data_size_fs, data_size_ss);
		return 1;
	}

	if ( (array->datatype == 'u')
	  && (array->itemsize == 2)
	  && (array->byteorder == '<') )
	{
		int fs, ss;
		uint16_t *in_data = (uint16_t *)array->data;

		for ( ss=0; ss<PANEL_HEIGHT(p); ss++ ) {
			for ( fs=0; fs<PANEL_WIDTH(p); fs++ ) {
				size_t idx = fs+p->orig_min_fs + (ss+p->orig_min_ss)*data_size_fs;
				data[fs+ss*PANEL_WIDTH(p)] = in_data[idx];
			}
		}

	} else if ( (array->datatype == 'u')
	  && (array->itemsize == 4)
	  && (array->byteorder == '<') )
	{
		int fs, ss;
		uint32_t *in_data = (uint32_t *)array->data;

		for ( ss=0; ss<PANEL_HEIGHT(p); ss++ ) {
			for ( fs=0; fs<PANEL_WIDTH(p); fs++ ) {
				size_t idx = fs+p->orig_min_fs + (ss+p->orig_min_ss)*data_size_fs;
				data[fs+ss*PANEL_WIDTH(p)] = in_data[idx];
			}
		}

	} else if ( (array->datatype == 'f')
	  && (array->itemsize == 8)
	  && (array->byteorder == '<') )
	{
		int fs, ss;
		double *in_data = (double *)array->data;

		for ( ss=0; ss<PANEL_HEIGHT(p); ss++ ) {
			for ( fs=0; fs<PANEL_WIDTH(p); fs++ ) {
				size_t idx = fs+p->orig_min_fs + (ss+p->orig_min_ss)*data_size_fs;
				data[fs+ss*PANEL_WIDTH(p)] = in_data[idx];
			}
		}

	} else {
		ERROR("Unrecognised data type %c%i%c\n",
		      array->datatype, array->itemsize, array->byteorder);
		return 1;
	}
	/* Note: when adding a new data type, don't forget that you are
	 * responsible for marking NaN/inf pixels as bad. */

	return 0;
}


/* Read the image data from 'data_block' into 'image', according to 'dtempl' */
int image_seedee_read(struct image *image,
                      DataTemplate *dtempl,
                      void *data_block,
                      size_t data_block_size,
                      char *meta_data)
{
	struct SeedeeNDArray array;
	int r;
	bool zero_copy;
	int i;
	cJSON *json;
	cJSON *data_format_str;

	json = cJSON_Parse(meta_data);
	if ( json == NULL ) {
		ERROR("Failed to parse JSON\n");
		return 1;
	}

	data_format_str = cJSON_GetObjectItemCaseSensitive(json, "_data_format");
	if ( !cJSON_IsString(data_format_str) ) {
		ERROR("_data_format isn't a string");
		cJSON_Delete(json);
		return 1;
	}

	profile_start("seedee-get-size");
	array.size = seedee_get_data_size(data_format_str->valuestring,
	                                  data_block, data_block_size,
	                                  &zero_copy, &array);
	profile_end("seedee-get-size");
	array.data = malloc(array.size);
	array.shape = malloc(array.ndims*sizeof(int));
	if ( (array.data == NULL) || (array.shape == NULL) ) {
		cJSON_Delete(json);
		free(array.data);
		free(array.shape);
		return 1;
	}

	if ( array.ndims != 2 ) {
		ERROR("Seedee data has unexpected number of dimensions "
		      "(%i, expected 2)\n", array.ndims);
		free(array.data);
		free(array.shape);
		return 1;
	}

	profile_start("seedee-deserialize");
	r = seedee_deserialize_ndarray(data_format_str->valuestring,
	                               data_block, data_block_size,
	                               0, &array);
	profile_end("seedee-deserialize");
	cJSON_Delete(json);
	if ( r < 0 ) {
		ERROR("Seedee deserialiation failed.\n");
		free(array.data);
		free(array.shape);
		return 1;
	}

	profile_start("seedee-panel");
	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( load_seedee_data(&dtempl->panels[i], &array,
		                      image->dp[i], image->bad[i]) )
		{
			ERROR("Failed to load data for panel '%s'\n",
			      dtempl->panels[i].name);
			profile_end("seedee-panel");
			free(array.data);
			free(array.shape);
			return 1;
		}
	}
	profile_end("seedee-panel");

	free(array.data);
	free(array.shape);

	return 0;
}


#else /* defined(HAVE_SEEDEE) */

int image_seedee_read(struct image *image,
                      const DataTemplate *dtempl,
                      void *data,
                      size_t data_size,
                      char *meta_data)
{
	ERROR("Seedee is not supported in this installation (read).\n");
	return 1;
}

#endif /* defined(HAVE_SEEDEE) */
