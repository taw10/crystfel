/*
 * zmq.c
 *
 * ZMQ data interface
 *
 * Copyright © 2017-2020 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018-2020 Thomas White <taw@physics.org>
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

#include <events.h>
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


static struct image *unpack_slab(const DataTemplate *dtempl,
                                 double *data,
                                 int data_width, int data_height)
{
	uint16_t *flags = NULL;
	float *sat = NULL;
	int pi;
	struct image *image;

	image = image_new();
	if ( image == NULL ) return NULL;

	image->dp = malloc(dtempl->n_panels*sizeof(float *));
	image->bad = malloc(dtempl->n_panels*sizeof(int *));
	image->sat = malloc(dtempl->n_panels*sizeof(float *));
	if ( (image->dp == NULL) || (image->bad == NULL) || (image->sat == NULL) ) {
		ERROR("Failed to allocate data arrays.\n");
		image_free(image);
		return NULL;
	}

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {

		struct panel_template *p;
		int fs, ss;
		int p_w, p_h;

		p = &dtempl->panels[pi];
		p_w = p->orig_max_fs - p->orig_min_fs + 1;
		p_h = p->orig_max_ss - p->orig_min_ss + 1;

		image->dp[pi] = malloc(p_w*p_h*sizeof(float));
		image->bad[pi] = malloc(p_w*p_h*sizeof(int));
		image->sat[pi] = malloc(p_w*p_h*sizeof(float));
		if ( (image->dp[pi] == NULL) || (image->bad[pi] == NULL)
		  || (image->sat[pi] == NULL) )
		{
			ERROR("Failed to allocate panel\n");
			return NULL;
		}

		if ( (p->orig_min_fs + p_w > data_width)
		  || (p->orig_min_ss + p_h > data_height) )
		{
			ERROR("Panel %s is outside range of data provided\n",
			      p->name);
			return NULL;
		}

		for ( ss=0; ss<p_h; ss++) {
		for ( fs=0; fs<p_w; fs++) {

			int idx;
			int cfs, css;
			int bad = 0;

			cfs = fs+p->orig_min_fs;
			css = ss+p->orig_min_ss;
			idx = cfs + css*data_width;

			image->dp[pi][fs+p_w*ss] = data[idx];

			if ( sat != NULL ) {
				image->sat[pi][fs+p_w*ss] = sat[idx];
			} else {
				image->sat[pi][fs+p_w*ss] = INFINITY;
			}

			if ( p->bad ) bad = 1;

			if ( data_template_in_bad_region(dtempl, pi,
			                                 fs, ss)
			     || isnan(image->dp[pi][fs+ss*p_w])
			     || isinf(image->dp[pi][fs+ss*p_w]) )
			{
				bad = 1;
			}

			if ( isnan(data[idx]) || isinf(data[idx]) ) bad = 1;

			if ( flags != NULL ) {

				int f;

				f = flags[idx];

				if ( (f & dtempl->mask_good)
				    != dtempl->mask_good ) bad = 1;

				if ( f & dtempl->mask_bad ) bad = 1;

			}
			image->bad[pi][fs+p_w*ss] = bad;
		}
		}

	}

	return 0;
}


static double *find_msgpack_data(msgpack_object *obj, int *width, int *height)
{
	msgpack_object *corr_data_obj;
	msgpack_object *data_obj;
	msgpack_object *shape_obj;
	double *data;

	corr_data_obj = find_msgpack_kv(obj, "corr_data");
	if ( corr_data_obj == NULL ) {
		ERROR("No corr_data MessagePack object found.\n");
		return NULL;
	}

	data_obj = find_msgpack_kv(corr_data_obj, "data");
	if ( data_obj == NULL ) {
		ERROR("No data MessagePack object found inside corr_data.\n");
		return NULL;
	}
	if ( data_obj->type != MSGPACK_OBJECT_STR ) {
		ERROR("corr_data.data isn't a binary object.\n");
		return NULL;
	}
	data = (double *)data_obj->via.str.ptr;

	shape_obj = find_msgpack_kv(corr_data_obj, "shape");
	if ( shape_obj == NULL ) {
		ERROR("No shape MessagePack object found inside corr_data.\n");
		return NULL;
	}
	if ( shape_obj->type != MSGPACK_OBJECT_ARRAY ) {
		ERROR("corr_data.shape isn't an array object.\n");
		return NULL;
	}
	if ( shape_obj->via.array.size != 2 ) {
		ERROR("corr_data.shape is wrong size (%i, should be 2)\n",
		      shape_obj->via.array.size);
		return NULL;
	}
	if ( shape_obj->via.array.ptr[0].type != MSGPACK_OBJECT_POSITIVE_INTEGER ) {
		ERROR("corr_data.shape contains wrong type of element.\n");
		return NULL;
	}
	*height = shape_obj->via.array.ptr[0].via.i64;
	*width = shape_obj->via.array.ptr[1].via.i64;
	return data;
}


static double *zero_array(const DataTemplate *dtempl, int *dw, int *dh)
{
	int max_fs = 0;
	int max_ss = 0;
	int pi;
	double *data;

	for ( pi=0; pi<dtempl->n_panels; pi++ ) {
		if ( dtempl->panels[pi].orig_max_fs > max_fs ) {
			max_fs = dtempl->panels[pi].orig_max_fs;
		}
		if ( dtempl->panels[pi].orig_max_ss > max_ss ) {
			max_ss = dtempl->panels[pi].orig_max_ss;
		}
	}

	data = calloc((max_fs+1)*(max_ss+1), sizeof(double));
	*dw = max_fs+1;
	*dh = max_ss+1;
	return data;
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
                                 msgpack_object *obj,
                                 int no_image_data)
{
	struct image *image;
	int data_width, data_height;
	double *data;

	if ( obj == NULL ) {
		ERROR("No MessagePack object!\n");
		return NULL;
	}

	if ( !no_image_data ) {
		data = find_msgpack_data(obj, &data_width, &data_height);
		if ( data == NULL ) {
			ERROR("No image data in MessagePack object.\n");
			return NULL;
		}
	} else {
		data = zero_array(dtempl, &data_width, &data_height);
	}

	image = unpack_slab(dtempl, data, data_width, data_height);

	if ( image == NULL ) {
		ERROR("Failed to unpack data slab.\n");
		return NULL;
	}

	return image;
}
