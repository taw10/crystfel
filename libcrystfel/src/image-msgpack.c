/*
 * image-msgpack.c
 *
 * Image loading, MessagePack parts
 *
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018-2021 Thomas White <taw@physics.org>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <zmq.h>
#include <msgpack.h>

#include <image.h>
#include <utils.h>
#include <msgpack.h>

#include "datatemplate_priv.h"


static msgpack_object *find_msgpack_kv(msgpack_object *obj, const char *key)
{
	int i;

	if ( obj == NULL ) return NULL;
	if ( obj->type != MSGPACK_OBJECT_MAP ) return NULL;

	for ( i=0; i<obj->via.map.size; i++ ) {
		const char *kstr;
		size_t klen;
		assert(obj->via.map.ptr[i].key.type == MSGPACK_OBJECT_STR);
		kstr = obj->via.map.ptr[i].key.via.str.ptr;
		klen = obj->via.map.ptr[i].key.via.str.size;
		if ( strncmp(kstr, key, klen) == 0 ) {
			return &obj->via.map.ptr[i].val;
		}
	}
	return NULL;
}


/**
 * image_msgpack_read_peaks
 * @dtempl: A %DataTemplate
 * @obj: A %msgpack_object containing data
 * @half_pixel_shift: Non-zero if 0.5 should be added to all peak coordinates
 *
 * Get peaks from msgpack_object. The data should be in a map, with the value
 * given by "peak_list" as an array of arrays. The first of these should contain
 * the list of fs positions of the peaks, the second the ss positions, and the
 * third the intensities of the peaks.
 *
 * http://c.msgpack.org/c/ provides documentation on msgpack objects
 *
 * CrystFEL considers all peak locations to be distances from the corner of the
 * detector panel, in pixel units, consistent with its description of detector
 * geometry (see 'man crystfel_geometry').  The software which generates the
 * CXI files, including Cheetah, may instead consider the peak locations to be
 * pixel indices in the data array.  In this case, the peak coordinates should
 * have 0.5 added to them.  This will be done if @half_pixel_shift is non-zero.
 *
 * Returns: a newly-allocated %ImageFeatureList.
 *
 */
ImageFeatureList *image_msgpack_read_peaks(const DataTemplate *dtempl,
                                           msgpack_object *obj,
                                           int half_pixel_shift)
{
	ImageFeatureList *features;
	int num_peaks;
	int pk;
	msgpack_object *peak_list;
	msgpack_object *peak_x;
	msgpack_object *peak_y;
	msgpack_object *peak_i;
	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( obj == NULL ) {
		ERROR("No MessagePack object to get peaks from.\n");
		return NULL;
	}

	/* Object has structure:
	 *   {
	 *    "peak_list": [[peak_x], [peak_y], [peak_i]]
	 *    "key2":val2,
	 *    ...
	 *   }
	 */
	peak_list = find_msgpack_kv(obj, "peak_list");
	peak_x = &peak_list->via.array.ptr[0];
	peak_y = &peak_list->via.array.ptr[1];
	peak_i = &peak_list->via.array.ptr[2];

	/* Length of peak_x  array gives number of peaks */
	num_peaks = peak_x->via.array.size;

	features = image_feature_list_new();

	for ( pk=0; pk<num_peaks; pk++ ) {

		float fs, ss, val;
		int pn;

		/* Retrieve data from peak_list and apply half_pixel_shift,
		 * if appropriate */
		fs = peak_x->via.array.ptr[pk].via.f64 + peak_offset;
		ss = peak_y->via.array.ptr[pk].via.f64 + peak_offset;
		val = peak_i->via.array.ptr[pk].via.f64;

		/* Convert coordinates to panel-relative */
		if ( data_template_file_to_panel_coords(dtempl, &fs, &ss, &pn) ) {
			ERROR("Peak not in panel!\n");
		} else {
			image_add_feature(features, fs, ss, pn,
			                  NULL, val, NULL);
		}
	}

	return features;
}


static int unpack_slab(struct image *image,
                       const DataTemplate *dtempl,
                       double *data,
                       int data_width, int data_height)
{
	int pi;

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate data arrays.\n");
		return 1;
	}

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		int fs, ss;
		int p_w, p_h;

		p = &dtempl->panels[pi];
		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->dp[pi] = malloc(p_w*p_h*sizeof(float));
		if ( image->dp[pi] == NULL ) {
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		if ( (p->orig_min_fs + p_w > data_width)
		  || (p->orig_min_ss + p_h > data_height) )
		{
			ERROR("Panel %s is outside range of data provided\n",
			      p->name);
			return 1;
		}

		for ( ss=0; ss<p_h; ss++) {
		for ( fs=0; fs<p_w; fs++) {

			int idx;
			int cfs, css;

			cfs = fs+p->orig_min_fs;
			css = ss+p->orig_min_ss;
			idx = cfs + css*data_width;

			image->dp[pi][fs+p_w*ss] = data[idx];

		}
		}

	}

	return 0;
}


static double *find_msgpack_data(msgpack_object *obj, int *width, int *height)
{
	FILE *fh = fopen("msgpack.data", "a");
	fprintf(fh, "object %p:\n", obj);
	msgpack_object_print(fh, *obj);
	fprintf(fh, "\n\n\n");
	fclose(fh);

	#if 0
	printf("Data type: %i\n", obj->type);
	if ( obj->type == MSGPACK_OBJECT_POSITIVE_INTEGER ) {
		printf("got an integer: %li\n", obj->via.i64);
	}

	if ( obj->type == MSGPACK_OBJECT_ARRAY ) {

		int i;
		printf("Array %i items\n", obj->via.array.size);

		for ( i=0; i<obj->via.array.size; i++ ) {
			msgpack_object *obj2 = obj->via.array.ptr[i];
			printf("Item %i: type %i\n", i, obj2->type);
			if ( obj2->type == MSGPACK_OBJECT_MAP ) {
				printf("Map: '%s' -> ");
			}
		}
	}
	#endif

	*width = 2068;
	*height = 2162;
	return NULL;
}


/* Unpacks the raw panel data from a msgpack_object, applies panel geometry,
 * and stores the resulting data in an image struct. Object has structure
 * {
 * "corr_data":
 *   {
 *     "data": binary_data,
 *     "shape": [data_height, data_width],
 *           ...
 *           ...
 *   },
 *   "key2": val2,
 *        ...
 *        ...
 * }
 */
struct image *image_msgpack_read(DataTemplate *dtempl,
                                 void *data,
                                 size_t data_size,
                                 int no_image_data,
                                 int no_mask_data)
{
	struct image *image;
	int data_width, data_height;
	double *image_data;
	msgpack_unpacked unpacked;
	int r;

	if ( data == NULL ) {
		ERROR("No MessagePack object!\n");
		return NULL;
	}

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return NULL;
	}

	msgpack_unpacked_init(&unpacked);
	r = msgpack_unpack_next(&unpacked, data, data_size, NULL);
	if ( r != MSGPACK_UNPACK_SUCCESS ) {
		ERROR("Msgpack unpack failed: %i\n", r);
		return NULL;
	}

	image = image_new();
	if ( image == NULL ) {
		ERROR("Couldn't allocate image structure.\n");
		return NULL;
	}

	if ( !no_image_data ) {
		image_data = find_msgpack_data(&unpacked.data,
		                               &data_width, &data_height);
		if ( image_data == NULL ) {
			ERROR("No image data in MessagePack object.\n");
			return NULL;
		}
		unpack_slab(image, dtempl, image_data,
		            data_width, data_height);
	} else {
		image_set_zero_data(image, dtempl);
	}

	image_set_zero_mask(image, dtempl);

	return image;
}
