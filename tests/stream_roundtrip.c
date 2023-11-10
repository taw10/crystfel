/*
 * stream_roundtrip.c
 *
 * Check that peaks and reflections can be sent on a stream round-trip
 *
 * Copyright © 2020-2022 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2022 Thomas White <taw@physics.org>
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


#include <stdio.h>
#include <unistd.h>

#include <stream.h>
#include <image.h>
#include <datatemplate.h>

#define N_PEAKS (10)

int main(int argc, char *argv[])
{
	gsl_rng *rng;
	int i;
	float peak_fs[N_PEAKS];
	float peak_ss[N_PEAKS];
	float peak_i[N_PEAKS];
	int peak_pn[N_PEAKS] = { 0, 0, 1, 1, 0, 1, 0, 1, 0, 1 };
	struct image *image;
	DataTemplate *dtempl;
	Stream *st;
	int fail = 0;

	/* Create test data ................................................. */

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	for ( i=0; i<N_PEAKS; i++ ) {
		peak_fs[i] = gsl_rng_uniform(rng) * 100.0;
		peak_ss[i] = gsl_rng_uniform(rng) * 100.0;
		peak_i[i] = 10.0*(i+1);
		STATUS("%f %f %i %f\n",
		       peak_fs[i], peak_ss[i], peak_pn[i], peak_i[i]);
	}
	gsl_rng_free(rng);

	/* Write stream ..................................................... */

	dtempl = data_template_new_from_file(argv[1]);
	if ( dtempl == NULL ) {
		ERROR("Failed to load data template\n");
		return 1;
	}

	image = image_create_for_simulation(dtempl);
	if ( image == NULL ) {
		ERROR("Failed to create image\n");
		return 1;
	}

	image->features = image_feature_list_new();
	for ( i=0; i<N_PEAKS; i++ ) {
		image_add_feature(image->features, peak_fs[i], peak_ss[i],
		                  peak_pn[i], peak_i[i], NULL);
	}

	st = stream_open_for_write("stream_roundtrip.stream", dtempl);
	if ( st == NULL ) {
		ERROR("Failed to open stream for writing\n");
		return 1;
	}

	stream_write_geometry_file(st, argv[1]);

	if ( stream_write_chunk(st, image, STREAM_PEAKS) ) {
		ERROR("Failed to write stream chunk\n");
		return 1;
	}

	stream_close(st);
	image_free(image);

	/* Read stream ...................................................... */

	st = stream_open_for_read("stream_roundtrip.stream");
	if ( st == NULL ) {
		ERROR("Failed to open stream for reading\n");
		return 1;
	}

	image = stream_read_chunk(st, STREAM_PEAKS);
	if ( image == NULL ) {
		ERROR("Failed to read stream chunk\n");
		return 1;
	}

	stream_close(st);

	for ( i=0; i<image_feature_count(image->features); i++ ) {
		struct imagefeature *f;
		f = image_get_feature(image->features, i);
		STATUS("%f %f %i %f\n", f->fs, f->ss, f->pn, f->intensity);
		if ( f->pn != peak_pn[i] ) {
			STATUS("Panel number doesn't match (should be %i)\n",
			       peak_pn[i]);
			fail = 1;
		}
		if ( !within_tolerance(f->fs, peak_fs[i], 0.1) ) {
			ERROR("fs doesn't match (should be %f)\n", peak_fs[i]);
			fail = 1;
		}
		if ( !within_tolerance(f->ss, peak_ss[i], 0.1) ) {
			ERROR("ss doesn't match (should be %f)\n", peak_ss[i]);
			fail = 1;
		}
		if ( !within_tolerance(f->intensity, peak_i[i], 0.1) ) {
			ERROR("Intensity doesn't match (should be %f)\n",
			      peak_i[i]);
			fail = 1;
		}
	}

	unlink("stream_roundtrip.stream");

	return fail;
}
