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


static msgpack_object *find_main_object(msgpack_unpacked *unpacked)
{
	int n_obj;

	if ( unpacked->data.type != MSGPACK_OBJECT_ARRAY ) {
		ERROR("MessagePack data isn't an array - ignoring.\n");
		return NULL;
	}

	n_obj = unpacked->data.via.array.size;
	if ( n_obj < 1 ) {
		ERROR("No array elements in MessagePack object?\n");
		return NULL;
	}
	if ( n_obj > 1 ) {
		ERROR("WARNING: Multiple (%i) items in MessagePack object - "
		      "ignoring all but the first one.\n", n_obj);
	}

	return &unpacked->data.via.array.ptr[0];
}


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


double image_msgpack_get_value(const char *name,
                               void *data_block,
                               size_t data_block_size,
                               char *ptype)
{
	msgpack_unpacked unpacked;
	msgpack_object *value_obj;
	msgpack_object *obj;
	int r;
	float val = NAN;

	*ptype = 'x';

	if ( data_block == NULL ) {
		ERROR("No MsgPack data!\n");
		goto out;
	}

	msgpack_unpacked_init(&unpacked);
	r = msgpack_unpack_next(&unpacked, data_block, data_block_size, NULL);
	if ( r != MSGPACK_UNPACK_SUCCESS ) {
		ERROR("MessagePack unpack failed: %i\n", r);
		goto out;
	}

	obj = find_main_object(&unpacked);
	if ( obj == NULL ) {
		ERROR("Failed to find main MsgPack object.\n");
		msgpack_unpacked_destroy(&unpacked);
		goto out;
	}

	value_obj = find_msgpack_kv(obj, name);
	if ( obj == NULL ) {
		ERROR("Couldn't find '%s' in MessagePack object\n", name);
		msgpack_unpacked_destroy(&unpacked);
		goto out;
	}

	switch ( value_obj->type ) {

		case MSGPACK_OBJECT_FLOAT64:
		case MSGPACK_OBJECT_FLOAT32:
		//case MSGPACK_OBJECT_FLOAT:
		*ptype = 'f';
		val = value_obj->via.f64;
		break;

		case MSGPACK_OBJECT_POSITIVE_INTEGER:
		case MSGPACK_OBJECT_NEGATIVE_INTEGER:
		*ptype = 'i';
		val = value_obj->via.i64;
		break;

		default:
		ERROR("Unrecognised MsgPack type %i\n", value_obj->type);
		break;

	}

	msgpack_unpacked_destroy(&unpacked);

out:
	return val;
}


static int load_msgpack_data(struct panel_template *p,
                             msgpack_object *map_obj,
                             float **data)
{
	msgpack_object *obj;

	obj = find_msgpack_kv(map_obj, p->data);
	if ( obj == NULL ) {
		ERROR("Couldn't find '%s' in MessagePack object\n",
		      p->data);
		return 1;
	}

	if ( obj->type != MSGPACK_OBJECT_MAP ) {
		ERROR("MessagePack object '%s' is not a map\n", p->data);
		return 1;
	}

	//type = get_msgpack_kv_string(obj, "type");
	//printf("data type is '%s'\n", type);

	//get_msgpack_kv_tuple(obj, "shape", &w, &h);

	//data_obj = find_msgpack_kv(obj, "data");

	return 0;
}


/* Read the image data from 'data_block' into 'image', according to 'dtempl' */
int image_msgpack_read(struct image *image,
                       DataTemplate *dtempl,
                       void *data_block,
                       size_t data_block_size)
{
	int i;
	int r;
	msgpack_unpacked unpacked;
	msgpack_object *obj;

	if ( dtempl == NULL ) {
		ERROR("NULL data template!\n");
		return 1;
	}

	if ( data_block == NULL ) {
		ERROR("No MsgPack data!\n");
		return 1;
	}

	msgpack_unpacked_init(&unpacked);
	r = msgpack_unpack_next(&unpacked, data_block, data_block_size, NULL);
	if ( r != MSGPACK_UNPACK_SUCCESS ) {
		ERROR("MessagePack unpack failed: %i\n", r);
		return 1;
	}

	obj = find_main_object(&unpacked);
	if ( obj == NULL ) {
		ERROR("Failed to find main MsgPack object.\n");
		msgpack_unpacked_destroy(&unpacked);
		return 1;
	}

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	if ( image->dp == NULL ) {
		ERROR("Failed to allocate data array.\n");
		msgpack_unpacked_destroy(&unpacked);
		return 1;
	}

	/* Set all pointers to NULL for easier clean-up */
	for ( i=0; i<dtempl->n_panels; i++ ) image->dp[i] = NULL;

	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( load_msgpack_data(&dtempl->panels[i], obj, &image->dp[i]) )
		{
			ERROR("Failed to load panel data\n");
			msgpack_unpacked_destroy(&unpacked);
			return 1;
		}
	}

	msgpack_unpacked_destroy(&unpacked);
	return 1;
}
