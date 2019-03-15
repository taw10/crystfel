/*
 * zmq.c
 *
 * ZMQ data interface
 *
 * Copyright © 2017-2018 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2018 Thomas White <taw@physics.org>
 *   2014 Valerio Mariani
 *   2017 Stijn de Graaf
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
#include <hdf5.h>
#include <assert.h>
#include <unistd.h>
#include <zmq.h>
#include <msgpack.h>

#include "events.h"
#include "image.h"
#include "hdf5-file.h"
#include "utils.h"


struct im_zmq
{
	void *ctx;
	void *socket;
	zmq_msg_t msg;
	msgpack_unpacked unpacked;
	int unpacked_set;
};


struct im_zmq *im_zmq_connect(const char *zmq_address)
{
	struct im_zmq *z;

	z = malloc(sizeof(struct im_zmq));
	if ( z == NULL ) return NULL;

	z->unpacked_set = 0;

	z->ctx = zmq_ctx_new();
	if ( z->ctx == NULL ) return NULL;

	z->socket = zmq_socket(z->ctx, ZMQ_REQ);
	if ( z->socket == NULL ) return NULL;

	STATUS("Connecting to ZMQ at '%s'\n", zmq_address);
	if ( zmq_connect(z->socket, zmq_address) == -1 ) {
		ERROR("ZMQ connection failed: %s\n", zmq_strerror(errno));
		return NULL;
	}
	STATUS("ZMQ connected.\n");

	return z;
}


msgpack_object *im_zmq_fetch(struct im_zmq *z)
{
	int msg_size;
	int r;

	if ( zmq_send(z->socket, "m", 1, 0) == -1 ) {
		ERROR("ZMQ message send failed: %s\n", zmq_strerror(errno));
		return NULL;
	}

	zmq_msg_init(&z->msg);
	msg_size = zmq_msg_recv(&z->msg, z->socket, 0);
	if ( msg_size == -1 ) {
		ERROR("ZMQ recieve failed: %s\n", zmq_strerror(errno));
		zmq_msg_close(&z->msg);
		return NULL;
	}

	msgpack_unpacked_init(&z->unpacked);
	r = msgpack_unpack_next(&z->unpacked, zmq_msg_data(&z->msg),
	                        msg_size, NULL);
	if ( r != MSGPACK_UNPACK_SUCCESS ) {
		ERROR("Msgpack unpack failed: %i\n", r);
		zmq_msg_close(&z->msg);
		return NULL;
	}
	z->unpacked_set = 1;

	return &z->unpacked.data;
}


/* Clean structures ready for next frame */
void im_zmq_clean(struct im_zmq *z)
{
	if ( z->unpacked_set ) {
		msgpack_unpacked_destroy(&z->unpacked);
		zmq_msg_close(&z->msg);
		z->unpacked_set = 0;
	}
}


void im_zmq_shutdown(struct im_zmq *z)
{
	if ( z == NULL ) return;

	zmq_msg_close(&z->msg);
	zmq_close(z->socket);
	zmq_ctx_destroy(z->ctx);
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
 * get_peaks_msgpack:
 * @obj: A %msgpack_object containing data in OnDA format
 * @image: An %image structure
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
 * Returns: non-zero on error, zero otherwise.
 *
 */
int get_peaks_msgpack(msgpack_object *obj, struct image *image,
                      int half_pixel_shift)
{

	int num_peaks;
	int pk;
	msgpack_object *peak_list;
	msgpack_object *peak_x;
	msgpack_object *peak_y;
	msgpack_object *peak_i;
	double peak_offset = half_pixel_shift ? 0.5 : 0.0;

	if ( obj == NULL ) {
		ERROR("No MessagePack object to get peaks from.\n");
		return 1;
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

	if ( image->features != NULL ) {
		image_feature_list_free(image->features);
	}
	image->features = image_feature_list_new();
	image->num_peaks = num_peaks;

	for ( pk=0; pk<num_peaks; pk++ ) {

		float fs, ss, val;
		struct panel *p;

		/* Retrieve data from peak_list and apply half_pixel_shift,
		 * if appropriate */
		fs = peak_x->via.array.ptr[pk].via.f64 + peak_offset;
		ss = peak_y->via.array.ptr[pk].via.f64 + peak_offset;
		val = peak_i->via.array.ptr[pk].via.f64;

		p = find_orig_panel(image->det, fs, ss);
		if ( p == NULL ) continue;
		if ( p->no_index ) continue;

		/* Convert coordinates to panel-relative */
		fs = fs - p->orig_min_fs;
		ss = ss - p->orig_min_ss;

		image_add_feature(image->features, fs, ss, p, image, val, NULL);
	}

	return 0;
}


static void im_zmq_fill_in_clen(struct detector *det)
{
	int i = 0;
	for ( i=0; i<det->n_panels; i++) {
		struct panel *p = &det->panels[i];
		if ( p->clen_from != NULL ) {
			ERROR("Can't get clen over ZMQ yet.\n");
		}
		adjust_centering_for_rail(p);
	}
}


static void im_zmq_fill_in_beam_parameters(struct beam_params *beam,
                                           struct image *image)
{
	double eV;
	if ( beam->photon_energy_from == NULL ) {
		/* Explicit value given */
		eV = beam->photon_energy;
	} else {
		ERROR("Can't get photon energy over ZMQ yet.\n");
		eV = 0.0;
	}
	image->lambda = ph_en_to_lambda(eV_to_J(eV))*beam->photon_energy_scale;
}


static int unpack_slab(struct image *image, double *data,
                       int data_width, int data_height)
{
	uint16_t *flags = NULL;
	float *sat = NULL;
	int pi;

	image->dp = malloc(image->det->n_panels*sizeof(float *));
	image->bad = malloc(image->det->n_panels*sizeof(int *));
	image->sat = malloc(image->det->n_panels*sizeof(float *));
	if ( (image->dp == NULL) || (image->bad == NULL) || (image->sat == NULL) ) {
		ERROR("Failed to allocate data arrays.\n");
		return 1;
	}

	for ( pi=0; pi<image->det->n_panels; pi++ ) {

		struct panel *p;
		int fs, ss;

		p = &image->det->panels[pi];
		image->dp[pi] = malloc(p->w*p->h*sizeof(float));
		image->bad[pi] = malloc(p->w*p->h*sizeof(int));
		image->sat[pi] = malloc(p->w*p->h*sizeof(float));
		if ( (image->dp[pi] == NULL) || (image->bad[pi] == NULL)
		  || (image->sat[pi] == NULL) )
		{
			ERROR("Failed to allocate panel\n");
			return 1;
		}

		if ( (p->orig_min_fs + p->w > data_width)
		  || (p->orig_min_ss + p->h > data_height) )
		{
			ERROR("Panel %s is outside range of data provided\n",
			      p->name);
			return 1;
		}

		for ( ss=0; ss<p->h; ss++) {
		for ( fs=0; fs<p->w; fs++) {

			int idx;
			int cfs, css;
			int bad = 0;

			cfs = fs+p->orig_min_fs;
			css = ss+p->orig_min_ss;
			idx = cfs + css*data_width;

			image->dp[pi][fs+p->w*ss] = data[idx];

			if ( sat != NULL ) {
				image->sat[pi][fs+p->w*ss] = sat[idx];
			} else {
				image->sat[pi][fs+p->w*ss] = INFINITY;
			}

			if ( p->no_index ) bad = 1;

			if ( in_bad_region(image->det, p, cfs, css) ) {
				bad = 1;
			}

			if ( flags != NULL ) {

				int f;

				f = flags[idx];

				if ( (f & image->det->mask_good)
				    != image->det->mask_good ) bad = 1;

				if ( f & image->det->mask_bad ) bad = 1;

			}
			image->bad[pi][fs+p->w*ss] = bad;
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


static double *zero_array(struct detector *det, int *dw, int *dh)
{
	int max_fs = 0;
	int max_ss = 0;
	int pi;
	double *data;

	for ( pi=0; pi<det->n_panels; pi++ ) {
		if ( det->panels[pi].orig_max_fs > max_fs ) {
			max_fs = det->panels[pi].orig_max_fs;
		}
		if ( det->panels[pi].orig_max_ss > max_ss ) {
			max_ss = det->panels[pi].orig_max_ss;
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
int unpack_msgpack_data(msgpack_object *obj, struct image *image,
                        int no_image_data)
{
	int data_width, data_height;
	double *data;

	if ( image->det == NULL ) {
		ERROR("Geometry not available.\n");
		return 1;
	}

	if ( obj == NULL ) {
		ERROR("No MessagePack object!\n");
		return 1;
	}

	if ( !no_image_data ) {
		data = find_msgpack_data(obj, &data_width, &data_height);
		if ( data == NULL ) {
			ERROR("No image data in MessagePack object.\n");
			return 1;
		}
	} else {
		data = zero_array(image->det, &data_width, &data_height);
	}

	if ( unpack_slab(image, data, data_width, data_height) ) {
		ERROR("Failed to unpack data slab.\n");
		return 1;
	}

	if ( image->beam != NULL ) {
		im_zmq_fill_in_beam_parameters(image->beam, image);
		if ( image->lambda > 1000 ) {
			ERROR("Warning: Missing or nonsensical wavelength "
			      "(%e m).\n", image->lambda);
		}
	}
	im_zmq_fill_in_clen(image->det);
	fill_in_adu(image);

	return 0;
}
