/*
 * process_image.c
 *
 * The processing pipeline for one image
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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
#include <hdf5.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <unistd.h>

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "thread-pool.h"
#include "geometry.h"
#include "stream.h"
#include "reflist-utils.h"
#include "process_image.h"
#include "integration.h"


static int cmpd2(const void *av, const void *bv)
{
	double a, b;

	a = *(double *)av;
	b = *(double *)bv;

	if ( fabs(a) < fabs(b) ) return -1;
	return 1;
}


static double *excitation_errors(UnitCell *cell, ImageFeatureList *flist,
                                 RefList *reflist, int *pnacc)
{
	int i;
	const double min_dist = 0.25;
	double *acc;
	int n_acc = 0;
	int n_notintegrated = 0;
	int max_acc = 1024;

	acc = malloc(max_acc*sizeof(double));
	if ( acc == NULL ) {
		ERROR("Allocation failed when refining radius!\n");
		return NULL;
	}

	for ( i=0; i<image_feature_count(flist); i++ ) {

		struct imagefeature *f;
		double h, k, l, hd, kd, ld;

		/* Assume all image "features" are genuine peaks */
		f = image_get_feature(flist, i);
		if ( f == NULL ) continue;

		double ax, ay, az;
		double bx, by, bz;
		double cx, cy, cz;

		cell_get_cartesian(cell,
		                   &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

		/* Decimal and fractional Miller indices of nearest
		 * reciprocal lattice point */
		hd = f->rx * ax + f->ry * ay + f->rz * az;
		kd = f->rx * bx + f->ry * by + f->rz * bz;
		ld = f->rx * cx + f->ry * cy + f->rz * cz;
		h = lrint(hd);
		k = lrint(kd);
		l = lrint(ld);

		/* Check distance */
		if ( (fabs(h - hd) < min_dist)
		  && (fabs(k - kd) < min_dist)
		  && (fabs(l - ld) < min_dist) )
		{
			double rlow, rhigh, p;
			Reflection *refl;

			/* Dig out the reflection */
			refl = find_refl(reflist, h, k, l);
			if ( refl == NULL ) {
				n_notintegrated++;
				continue;
			}

			get_partial(refl, &rlow, &rhigh, &p);
			acc[n_acc++] = fabs(rlow+rhigh)/2.0;
			if ( n_acc == max_acc ) {
				max_acc += 1024;
				acc = realloc(acc, max_acc*sizeof(double));
				if ( acc == NULL ) {
					ERROR("Allocation failed during"
					      " estimate_resolution!\n");
					return NULL;
				}
			}
		}

	}

	if ( n_acc < 3 ) {
		STATUS("WARNING: Too few peaks to estimate profile radius.\n");
		return NULL;
	}

	*pnacc = n_acc;
	return acc;
}


static void refine_radius(Crystal *cr, ImageFeatureList *flist)
{
	int n = 0;
	int n_acc;
	double *acc;

	acc = excitation_errors(crystal_get_cell(cr), flist,
	                        crystal_get_reflections(cr), &n_acc);
	if ( acc == NULL ) return;

	qsort(acc, n_acc, sizeof(double), cmpd2);
	n = n_acc/50;
	if ( n < 2 ) n = 2; /* n_acc is always >= 2 */
	crystal_set_profile_radius(cr, acc[(n_acc-1)-n]);

	free(acc);
}


void process_image(const struct index_args *iargs, struct pattern_args *pargs,
                   Stream *st, int cookie, const char *tmpdir, int results_pipe,
                   int serial)
{
	float *data_for_measurement;
	size_t data_size;
	int check;
	struct hdfile *hdfile;
	struct image image;
	int i;
	int r;
	int ret;
	char *rn;

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.copyme = iargs->copyme;
	image.id = cookie;
	image.filename = pargs->filename_p_e->filename;
	image.event = pargs->filename_p_e->ev;
	image.beam = iargs->beam;
	image.det = copy_geom(iargs->det);
	image.crystals = NULL;
	image.n_crystals = 0;
	image.serial = serial;
	image.indexed_by = INDEXING_NONE;

	hdfile = hdfile_open(image.filename);
	if ( hdfile == NULL ) {
		ERROR("Couldn't open file: %s\n", image.filename);
		return;
	}

	check = hdf5_read2(hdfile, &image, image.event, 0);
	if ( check ) {
		return;
	}

	/* Take snapshot of image after CM subtraction but before applying
	 * horrible noise filters to it */
	data_size = image.width * image.height * sizeof(float);
	data_for_measurement = malloc(data_size);
	memcpy(data_for_measurement, image.data, data_size);

	if ( iargs->median_filter > 0 ) {
		filter_median(&image, iargs->median_filter);
	}

	if ( iargs->noisefilter ) {
		filter_noise(&image);
	}

	mark_resolution_range_as_bad(&image, iargs->highres, +INFINITY);

	switch ( iargs->peaks ) {

		case PEAK_HDF5:
		if ( get_peaks(&image, hdfile, iargs->hdf5_peak_path) ) {
			ERROR("Failed to get peaks from HDF5 file.\n");
		}
		if ( !iargs->no_revalidate ) {
			validate_peaks(&image, iargs->min_snr,
				       iargs->pk_inn, iargs->pk_mid,
		                       iargs->pk_out, iargs->use_saturated,
				       iargs->check_hdf5_snr);
		}
		break;

		case PEAK_CXI:
		if ( get_peaks_cxi(&image, hdfile, iargs->hdf5_peak_path,
		                   pargs->filename_p_e) ) {
			ERROR("Failed to get peaks from CXI file.\n");
		}
		if ( !iargs->no_revalidate ) {
			validate_peaks(&image, iargs->min_snr,
				       iargs->pk_inn, iargs->pk_mid,
		                       iargs->pk_out, iargs->use_saturated,
				       iargs->check_hdf5_snr);
		}
		break;

		case PEAK_ZAEF:
		search_peaks(&image, iargs->threshold,
		             iargs->min_gradient, iargs->min_snr,
		             iargs->pk_inn, iargs->pk_mid,iargs->pk_out,
		             iargs->use_saturated);
		break;

	}

	/* Get rid of noise-filtered version at this point
	 * - it was strictly for the purposes of peak detection. */
	free(image.data);
	image.data = data_for_measurement;

	rn = getcwd(NULL, 0);

	r = chdir(tmpdir);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		hdfile_close(hdfile);
		return;
	}

	/* Index the pattern */
	index_pattern(&image, iargs->indm, iargs->ipriv);

	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		hdfile_close(hdfile);
		return;
	}
	free(rn);

	pargs->n_crystals = image.n_crystals;
	for ( i=0; i<image.n_crystals; i++ ) {
		crystal_set_image(image.crystals[i], &image);
	}

	/* Set beam/crystal parameters */
	if ( iargs->fix_divergence >= 0.0 ) {
		image.div = iargs->fix_divergence;
	} else {
		image.div = 0.0;
	}
	if ( iargs->fix_bandwidth >= 0.0 ) {
		image.bw = iargs->fix_bandwidth;
	} else {
		image.bw = 0.00000001;
	}
	if ( iargs->fix_profile_r >= 0.0 ) {
		for ( i=0; i<image.n_crystals; i++ ) {
			crystal_set_profile_radius(image.crystals[i],
			                           iargs->fix_profile_r);
			crystal_set_mosaicity(image.crystals[i], 0.0);
		}
	} else {
		for ( i=0; i<image.n_crystals; i++ ) {
			crystal_set_profile_radius(image.crystals[i], 0.02e9);
			crystal_set_mosaicity(image.crystals[i], 0.0);
		}
	}

	/* Integrate all the crystals at once - need all the crystals so that
	 * overlaps can be detected. */
	if ( iargs->fix_profile_r < 0.0 ) {

		integrate_all_4(&image, iargs->int_meth, PMODEL_SCSPHERE,
		                iargs->push_res,
		                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
		                INTDIAG_NONE, 0, 0, 0, results_pipe);

		for ( i=0; i<image.n_crystals; i++ ) {
			refine_radius(image.crystals[i], image.features);
			reflist_free(crystal_get_reflections(image.crystals[i]));
		}

		integrate_all_4(&image, iargs->int_meth, PMODEL_SCSPHERE,
		                iargs->push_res,
		                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
		                iargs->int_diag, iargs->int_diag_h,
		                iargs->int_diag_k, iargs->int_diag_l,
		                results_pipe);
	} else {

		integrate_all_4(&image, iargs->int_meth, PMODEL_SCSPHERE,
		                iargs->push_res,
		                iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
		                iargs->int_diag, iargs->int_diag_h,
		                iargs->int_diag_k, iargs->int_diag_l,
		                results_pipe);

	}

	ret = write_chunk(st, &image, hdfile,
	                  iargs->stream_peaks, iargs->stream_refls,
	                  pargs->filename_p_e->ev);
	if ( ret != 0 ) {
		ERROR("Error writing stream file.\n");
	}

	int n = 0;
	for ( i=0; i<image.n_crystals; i++ ) {
		n += crystal_get_num_implausible_reflections(image.crystals[i]);
	}
	if ( n > 0 ) {
		STATUS("WARNING: %i implausibly negative reflection%s in %s "
		       "%s\n", n, n>1?"s":"", image.filename,
		       get_event_string(image.event));
	}

	for ( i=0; i<image.n_crystals; i++ ) {
		cell_free(crystal_get_cell(image.crystals[i]));
		reflist_free(crystal_get_reflections(image.crystals[i]));
		crystal_free(image.crystals[i]);
	}
	free(image.crystals);

	for ( i=0; i<image.det->n_panels; i++ ) {
		free(image.dp[i]);
		free(image.bad[i]);
	}
	free(image.dp);
	free(image.bad);

	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	image_feature_list_free(image.features);
	free_detector_geometry(image.det);
	hdfile_close(hdfile);
}
