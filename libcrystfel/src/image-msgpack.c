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
		const char *kstr;
		size_t klen;
		if ( (obj->via.map.ptr[i].key.type != MSGPACK_OBJECT_STR)
		  && (obj->via.map.ptr[i].key.type != MSGPACK_OBJECT_BIN) ) continue;
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
		if ( data_template_file_to_panel_coords(dtempl, &fs, &ss, &pn) ) {
			ERROR("Peak not in panel!\n");
		} else {
			image_add_feature(features, fs, ss, pn,
			                  NULL, val, NULL);
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
	str = malloc(len+1);
	if ( str == NULL ) return NULL;
	strncpy(str, ptr, len);
	str[len] = '\0';
	return str;
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
	char *str;

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

		case MSGPACK_OBJECT_STR:
		str = terminate_str(value_obj->via.str.ptr,
		                    value_obj->via.str.size);
		if ( str != NULL ) {
			int ival;
			if ( convert_int(str, &ival) == 0 ) {
				*ptype = 'i';
				val = ival;
			} else {
				ERROR("MsgPack header %s has a string type (%s)"
				      "(need a number, and can't convert).\n",
				      name, str);
				val = NAN;
			}
			free(str);
		} else {
			ERROR("Failed to read MsgPack string (%s)\n", name);
			val = NAN;
		}
		break;

		default:
		ERROR("Unrecognised MsgPack type %i (%s)\n",
		      value_obj->type, name);
		break;

	}

	msgpack_unpacked_destroy(&unpacked);

out:
	return val;
}


static int load_msgpack_data(struct panel_template *p,
                             msgpack_object *map_obj,
                             float **pdata)
{
	msgpack_object *obj;
	msgpack_object *type_obj;
	msgpack_object *shape_obj;
	msgpack_object *data_obj;
	char *dtype;
	int data_size_fs, data_size_ss;
	void *data = NULL;

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
	dtype = strndup(type_obj->via.str.ptr, type_obj->via.str.size);

	shape_obj = find_msgpack_kv(obj, "shape");
	if ( shape_obj == NULL ) {
		ERROR("Data 'shape' not found\n");
		free(dtype);
		return 1;
	}
	if ( shape_obj->type != MSGPACK_OBJECT_ARRAY ) {
		ERROR("Data 'shape' isn't an array\n");
		free(dtype);
		return 1;
	}
	if ( shape_obj->via.array.size != 2 ) {
		ERROR("Data 'shape' has wrong number of dimensions (%i)\n",
		      shape_obj->via.array.size);
		free(dtype);
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
		free(dtype);
		return 1;
	}
	if ( data_obj->type != MSGPACK_OBJECT_BIN ) {
		ERROR("Data 'data' isn't binary\n");
		free(dtype);
		return 1;
	}

	if ( strcmp(dtype, "<i4") == 0 ) {

		int fs, ss;
		float *fdata;

		fdata = malloc(PANEL_WIDTH(p) * PANEL_HEIGHT(p) * sizeof(float));
		if ( fdata == NULL ) return 1;

		for ( ss=0; ss<PANEL_HEIGHT(p); ss++ ) {
			for ( fs=0; fs<PANEL_WIDTH(p); fs++ ) {
				size_t idx = fs+p->orig_min_fs + (ss+p->orig_min_ss)*data_size_fs;
				fdata[fs+ss*PANEL_WIDTH(p)] = data_obj->via.bin.ptr[idx];
			}
		}

		data = fdata;

	} else {
		ERROR("Unrecognised data type '%s'\n", dtype);
	}

	free(dtype);
	*pdata = data;
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

double image_msgpack_get_value(const char *name,
                               void *data_block,
                               size_t data_block_size,
                               char *ptype)
{
	ERROR("MessagePack is not supported in this installation (get_value).\n");
	*ptype = 'f';
	return NAN;
}



#endif /* defined(HAVE_MSGPACK) */
