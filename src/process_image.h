/*
 * process_image.h
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2023 Thomas White <taw@physics.org>
 *   2014-2017 Valerio Mariani <valerio.mariani@desy.de>
 *   2017-2018 Yaroslav Gevorkov <yaroslav.gevorkov@desy.de>
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

#ifndef PROCESS_IMAGE_H
#define PROCESS_IMAGE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct index_args;

#ifdef HAVE_MSGPACK
#include <msgpack.h>
#endif

#include "integration.h"
#include "im-sandbox.h"
#include "peaks.h"
#include "image.h"
#include "im-asapo.h"
#include "predict-refine.h"


/* Information about the indexing process which is common to all patterns */
struct index_args
{
	/* Input */
	DataTemplate *dtempl;
	signed int wait_for_file; /* -1 means wait forever */
	int no_image_data;
	int no_mask_data;
	float highres;
	DataSourceType data_format;

	/* Peak search */
	struct peak_params peak_search;
	void *pf_private;

	/* Hit finding */
	int min_peaks;

	/* Indexing */
	IndexingPrivate *ipriv;
	UnitCell *cell;
	float tols[6];
	float wavelength_estimate;
	float clen_estimate;
	int n_threads;
	int mille;
	int max_mille_level;

	/* Integration */
	IntegrationMethod int_meth;
	IntDiag int_diag;
	float ir_inn;
	float ir_mid;
	float ir_out;
	signed int int_diag_h;
	signed int int_diag_k;
	signed int int_diag_l;
	float push_res;
	float fix_profile_r;
	float fix_divergence;
	int overpredict;
	int cell_params_only;

	/* Output */
	int stream_flags;
	int stream_nonhits;
};


/* Information about the indexing process for one pattern */
struct pattern_args
{
	/* "Input" */
	char *filename;
	char *event;

	void *zmq_data;
	size_t zmq_data_size;

	char *asapo_data;
	size_t asapo_data_size;
	char *asapo_meta;
};


extern void process_image(const struct index_args *iargs,
                          struct pattern_args *pargs, Stream *st,
                          int cookie, const char *tmpdir, int serial,
                          struct sb_shm *sb_shared, char *last_task,
                          struct im_asapo *asapostuff,
                          Mille *mille, ImageDataArrays *ida);


#endif	/* PROCESS_IMAGE_H */
