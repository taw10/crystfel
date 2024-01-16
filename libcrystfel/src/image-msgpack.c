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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>

#include <image.h>
#include <utils.h>

#include "datatemplate_priv.h"


#if defined(HAVE_MSGPACK)

#include <msgpack.h>

static msgpack_object *find_main_object(msgpack_unpacked *unpacked)
{
	int n_obj;

	if ( unpacked->data.type == MSGPACK_OBJECT_MAP ) {
		return &unpacked->data;
	}

	if ( unpacked->data.type != MSGPACK_OBJECT_ARRAY ) {
		ERROR("MessagePack data isn't an array - ignoring.\n");
		return NULL;
	}

	n_obj = unpacked->data.via.array.size;
	if ( n_obj < 1 ) {
		ERROR("No array elements in MessagePack object?\n");
		return NULL;
	}

	return &unpacked->data.via.array.ptr[0];
}


static msgpack_object *find_msgpack_kv(msgpack_object *obj, const char *key)
{
	int i;

	if ( obj == NULL ) return NULL;
	if ( obj->type != MSGPACK_OBJECT_MAP ) return NULL;

	for ( i=0; i<obj->via.map.size; i++ ) {

		if ( obj->via.map.ptr[i].key.type == MSGPACK_OBJECT_STR ) {

			const char *kstr;
			size_t klen;

			kstr = obj->via.map.ptr[i].key.via.str.ptr;
			klen = obj->via.map.ptr[i].key.via.str.size;
			if ( strncmp(kstr, key, klen) == 0 ) {
				return &obj->via.map.ptr[i].val;
			}

		} else if ( obj->via.map.ptr[i].key.type == MSGPACK_OBJECT_BIN ) {

			const char *kstr;
			size_t klen;

			kstr = obj->via.map.ptr[i].key.via.bin.ptr;
			klen = obj->via.map.ptr[i].key.via.bin.size;
			if ( strncmp(kstr, key, klen) == 0 ) {
				return &obj->via.map.ptr[i].val;
			}
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
                                           void *data_block,
                                           size_t data_block_size,
                                           int half_pixel_shift)
{
	msgpack_unpacked unpacked;
	msgpack_object *obj;
	int r;
	ImageFeatureList *features;
	int num_peaks;
	int pk;
	msgpack_object *peak_list;
	msgpack_object *peak_fs;
	msgpack_object *peak_ss;
	msgpack_object *peak_i;
	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( data_block == NULL ) {
		ERROR("No MsgPack data!\n");
		return NULL;
	}

	if ( dtempl->peak_list == NULL ) {
		ERROR("You must set 'peak_list' in the geometry file.\n");
		return NULL;
	}

	msgpack_unpacked_init(&unpacked);
	r = msgpack_unpack_next(&unpacked, data_block, data_block_size, NULL);
	if ( r != MSGPACK_UNPACK_SUCCESS ) {
		ERROR("MessagePack unpack failed: %i (read_peaks)\n", r);
		return NULL;
	}

	obj = find_main_object(&unpacked);
	if ( obj == NULL ) {
		ERROR("Failed to find main MsgPack object (read_peaks).\n");
		msgpack_unpacked_destroy(&unpacked);
		return NULL;
	}

	peak_list = find_msgpack_kv(obj, dtempl->peak_list);
	if ( peak_list == NULL ) {
		ERROR("Peak list object not found ('%s')\n", dtempl->peak_list);
		msgpack_unpacked_destroy(&unpacked);
		return NULL;
	}

	peak_fs = find_msgpack_kv(peak_list, "fs");
	peak_ss = find_msgpack_kv(peak_list, "ss");
	peak_i = find_msgpack_kv(peak_list, "intensity");

	if ( (peak_fs == NULL)
	  || (peak_ss == NULL)
	  || (peak_i == NULL)
	  || (peak_fs->type != MSGPACK_OBJECT_ARRAY)
	  || (peak_ss->type != MSGPACK_OBJECT_ARRAY)
	  || (peak_i->type != MSGPACK_OBJECT_ARRAY) )
	{
		ERROR("Peak list has wrong structure.\n");
		msgpack_unpacked_destroy(&unpacked);
		return NULL;
	}

	/* Length of peak_x  array gives number of peaks */
	num_peaks = peak_fs->via.array.size;
	if ( (peak_ss->via.array.size != num_peaks)
	  || (peak_i->via.array.size != num_peaks) )
	{
		ERROR("Peak list arrays do not have equal size.\n");
		msgpack_unpacked_destroy(&unpacked);
		return NULL;
	}

	features = image_feature_list_new();

	for ( pk=0; pk<num_peaks; pk++ ) {

		float fs, ss, val;
		int pn;

		/* Retrieve data from peak_list and apply half_pixel_shift,
		 * if appropriate */
		fs = peak_fs->via.array.ptr[pk].via.f64 + peak_offset;
		ss = peak_ss->via.array.ptr[pk].via.f64 + peak_offset;
		val = peak_i->via.array.ptr[pk].via.f64;

		/* Convert coordinates to panel-relative */
		if ( data_template_slabby_file_to_panel_coords(dtempl, &fs, &ss, &pn) ) {
			ERROR("Failed to convert %i,%i to "
			      "panel-relative coordinates\n", fs, ss);
		} else {
			image_add_feature(features, fs, ss, pn, val, NULL);
		}
	}

	msgpack_unpacked_destroy(&unpacked);
	return features;
}


static char *terminate_str(const char *ptr, size_t len)
{
	char *str;
	if ( len < 1 ) return NULL;
	if ( len > 16*1024 ) return NULL;
	str = cfmalloc(len+1);
	if ( str == NULL ) return NULL;
	strncpy(str, ptr, len);
	str[len] = '\0';
	return str;
}


int image_msgpack_read_header_to_cache(struct image *image,
                                       const char *name)
{
	msgpack_unpacked unpacked;
	msgpack_object *value_obj;
	msgpack_object *obj;
	int r;
	char *str;

	if ( image->data_block == NULL ) {
		ERROR("No MsgPack data!\n");
		return 1;
	}

	msgpack_unpacked_init(&unpacked);
	r = msgpack_unpack_next(&unpacked, image->data_block,
	                        image->data_block_size, NULL);
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

	value_obj = find_msgpack_kv(obj, name);
	if ( value_obj == NULL ) {
		ERROR("Couldn't find '%s' in MessagePack object\n", name);
		msgpack_unpacked_destroy(&unpacked);
		return 1;
	}

	switch ( value_obj->type ) {

		case MSGPACK_OBJECT_FLOAT64:
		case MSGPACK_OBJECT_FLOAT32:
		image_cache_header_float(image, name, value_obj->via.f64);
		msgpack_unpacked_destroy(&unpacked);
		return 0;

		case MSGPACK_OBJECT_POSITIVE_INTEGER:
		case MSGPACK_OBJECT_NEGATIVE_INTEGER:
		image_cache_header_int(image, name, value_obj->via.i64);
		msgpack_unpacked_destroy(&unpacked);
		return 0;

		case MSGPACK_OBJECT_STR:
		str = terminate_str(value_obj->via.str.ptr,
		                    value_obj->via.str.size);
		if ( str == NULL ) {
			msgpack_unpacked_destroy(&unpacked);
			return 1;
		}

		image_cache_header_str(image, name, str);
		cffree(str);
		msgpack_unpacked_destroy(&unpacked);
		return 0;

		default:
		ERROR("Unrecognised MsgPack type %i (%s)\n",
		      value_obj->type, name);
		msgpack_unpacked_destroy(&unpacked);
		return 1;

	}
}


static int load_msgpack_data(struct panel_template *p,
                             msgpack_object *map_obj,
                             float *data, int *bad)
{
	msgpack_object *obj;
	msgpack_object *type_obj;
	msgpack_object *shape_obj;
	msgpack_object *data_obj;
	char *dtype;
	int data_size_fs, data_size_ss;
	int i;

	for ( i=0; i<MAX_DIMS; i++ ) {
		if ( (p->dims[i] >= 0) || (p->dims[i] == DIM_PLACEHOLDER) ) {
			ERROR("Only a single 2D array is supported via MsgPack.\n");
			ERROR("Check the geometry file and remove 'dimX = %' and 'dimX = <n>'\n");
			return 1;
		}
	}

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

	type_obj = find_msgpack_kv(obj, "type");
	if ( type_obj == NULL ) {
		ERROR("Data 'type' not found\n");
		return 1;
	}
	if ( type_obj->type != MSGPACK_OBJECT_STR ) {
		ERROR("Data 'type' isn't a string\n");
		return 1;
	}
	dtype = cfstrndup(type_obj->via.str.ptr, type_obj->via.str.size);

	shape_obj = find_msgpack_kv(obj, "shape");
	if ( shape_obj == NULL ) {
		ERROR("Data 'shape' not found\n");
		cffree(dtype);
		return 1;
	}
	if ( shape_obj->type != MSGPACK_OBJECT_ARRAY ) {
		ERROR("Data 'shape' isn't an array\n");
		cffree(dtype);
		return 1;
	}
	if ( shape_obj->via.array.size != 2 ) {
		ERROR("Data 'shape' has wrong number of dimensions (%i)\n",
		      shape_obj->via.array.size);
		cffree(dtype);
		return 1;
	}
	data_size_ss = shape_obj->via.array.ptr[0].via.u64;
	data_size_fs = shape_obj->via.array.ptr[1].via.u64;

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

	data_obj = find_msgpack_kv(obj, "data");
	if ( data_obj == NULL ) {
		ERROR("Data 'data' not found\n");
		cffree(dtype);
		return 1;
	}
	if ( data_obj->type != MSGPACK_OBJECT_BIN ) {
		ERROR("Data 'data' isn't binary\n");
		cffree(dtype);
		return 1;
	}

	if ( strcmp(dtype, "<i4") == 0 ) {

		int fs, ss;
		int32_t *in_data = (int32_t *)data_obj->via.bin.ptr;

		for ( ss=0; ss<PANEL_HEIGHT(p); ss++ ) {
			for ( fs=0; fs<PANEL_WIDTH(p); fs++ ) {
				size_t idx = fs+p->orig_min_fs + (ss+p->orig_min_ss)*data_size_fs;
				data[fs+ss*PANEL_WIDTH(p)] = in_data[idx];
				/* Integer data -> no need to check NaN/inf */
			}
		}

	} else if ( strcmp(dtype, "<f4") == 0 ) {

		int fs, ss;
		float *in_data = (float *)data_obj->via.bin.ptr;

		for ( ss=0; ss<PANEL_HEIGHT(p); ss++ ) {
			for ( fs=0; fs<PANEL_WIDTH(p); fs++ ) {
				size_t idx = fs+p->orig_min_fs + (ss+p->orig_min_ss)*data_size_fs;
				data[fs+ss*PANEL_WIDTH(p)] = in_data[idx];
				if ( !isfinite(in_data[idx]) ) {
					bad[fs+ss*PANEL_WIDTH(p)] = 1;
				}
			}
		}

	} else {
		ERROR("Unrecognised data type '%s'\n", dtype);
	}

	cffree(dtype);
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

	for ( i=0; i<dtempl->n_panels; i++ ) {
		if ( load_msgpack_data(&dtempl->panels[i], obj,
		                       image->dp[i], image->bad[i]) )
		{
			ERROR("Failed to load data for panel '%s'\n",
			      dtempl->panels[i].name);
			msgpack_unpacked_destroy(&unpacked);
			return 1;
		}
	}

	msgpack_unpacked_destroy(&unpacked);
	return 0;
}


#else /* defined(HAVE_MSGPACK) */

int image_msgpack_read(struct image *image,
                       const DataTemplate *dtempl,
                       void *data,
                       size_t data_size)
{
	ERROR("MessagePack is not supported in this installation (read).\n");
	return 1;
}

ImageFeatureList *image_msgpack_read_peaks(const DataTemplate *dtempl,
                                           void *data_block,
                                           size_t data_block_size,
                                           int half_pixel_shift)
{
	ERROR("MessagePack is not supported in this installation (read_peaks).\n");
	return NULL;
}

int image_msgpack_read_header_to_cache(struct image *image,
                                       const char *name)
{
	ERROR("MessagePack is not supported in this installation"
	      " (read_header_to_cache).\n");
	return 1;
}



#endif /* defined(HAVE_MSGPACK) */
