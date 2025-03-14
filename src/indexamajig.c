/*
 * indexamajig.c
 *
 * Index patterns, output hkl+intensity etc.
 *
 * Copyright © 2012-2023 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2023 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
 *   2017      Valerio Mariani <valerio.mariani@desy.de>
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SCHED_SETAFFINITY
#define _GNU_SOURCE
#include <sys/sysinfo.h>
#include <sched.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/mman.h>

#ifdef HAVE_CJSON
#include <cjson/cJSON.h>
#endif

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include <utils.h>
#include <index.h>
#include <peaks.h>
#include <filters.h>
#include <thread-pool.h>
#include <geometry.h>
#include <stream.h>
#include <reflist-utils.h>
#include <cell-utils.h>
#include <integration.h>
#include <image.h>
#include <datatemplate.h>
#include <detgeom.h>
#include <peakfinder8.h>

#include "im-sandbox.h"
#include "im-argparse.h"
#include "im-zmq.h"
#include "im-asapo.h"
#include "version.h"
#include "profile.h"


#ifdef HAVE_CJSON

static double nan_if_neg(double n)
{
	if ( n < 0.0 ) {
		return NAN;
	} else {
		return n;
	}
}


static void write_json_cell(cJSON *gp, const char *name, UnitCell *cell)
{
	double a, b, c, al, be, ga;
	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	if ( cell == NULL ) {
		cJSON_AddNullToObject(gp, name);
	} else {

		cJSON *obj = cJSON_AddObjectToObject(gp, name);
		char tmp[2];

		cJSON_AddStringToObject(obj, "lattice_type", str_lattice(cell_get_lattice_type(cell)));

		tmp[0] = cell_get_centering(cell);   tmp[1] = '\0';
		cJSON_AddStringToObject(obj, "centering", tmp);

		tmp[0] = cell_get_unique_axis(cell);
		cJSON_AddStringToObject(obj, "unique_axis", tmp);

		cJSON_AddNumberToObject(obj, "a_Angstrom", a*1e10);
		cJSON_AddNumberToObject(obj, "b_Angstrom", b*1e10);
		cJSON_AddNumberToObject(obj, "c_Angstrom", c*1e10);
		cJSON_AddNumberToObject(obj, "alpha_deg", rad2deg(al));
		cJSON_AddNumberToObject(obj, "beta_deg", rad2deg(be));
		cJSON_AddNumberToObject(obj, "gamma_deg", rad2deg(ga));
	}
}


static void write_methods(cJSON *gp, const char *name, IndexingPrivate *ipriv)
{
	cJSON *arr = cJSON_AddArrayToObject(gp, name);

	if ( ipriv != NULL ) {
		const IndexingMethod *methods;
		int i, n;
		methods = indexing_methods(ipriv, &n);
		for ( i=0; i<n; i++ ) {
			cJSON_AddItemToArray(arr, cJSON_CreateString(indexer_str(methods[i])));
		}
	}
}

#endif /* HAVE_CJSON */


static void write_harvest_file(struct index_args *args,
                               const char *filename,
                               int if_multi, int if_refine, int if_retry,
                               int if_peaks, int if_checkcell)
{
	#ifdef HAVE_CJSON
	FILE *fh;
	cJSON *harvest;
	cJSON *gp;
	cJSON *radii;

	harvest = cJSON_CreateObject();
	if ( harvest == NULL ) return;

	gp = cJSON_AddObjectToObject(harvest, "input");
	cJSON_AddNumberToObject(gp, "highres", args->highres);

	gp = cJSON_AddObjectToObject(harvest, "peaksearch");
	cJSON_AddStringToObject(gp, "method", str_peaksearch(args->peak_search.method));

	radii = cJSON_AddObjectToObject(gp, "radii");
	cJSON_AddNumberToObject(radii, "inner_px", args->peak_search.pk_inn);
	cJSON_AddNumberToObject(radii, "middle_px",  args->peak_search.pk_mid);
	cJSON_AddNumberToObject(radii, "outer_px",  args->peak_search.pk_out);

	cJSON_AddBoolToObject(gp, "noise_filter", args->peak_search.noisefilter);
	cJSON_AddNumberToObject(gp, "median_filter", args->peak_search.median_filter);
	cJSON_AddNumberToObject(gp, "threshold_adu", args->peak_search.threshold);
	cJSON_AddNumberToObject(gp, "min_squared_gradient_adu2", args->peak_search.min_sq_gradient);
	cJSON_AddNumberToObject(gp, "min_snr", args->peak_search.min_snr);
	cJSON_AddBoolToObject(gp, "check_hdf5_snr", args->peak_search.check_hdf5_snr);
	cJSON_AddBoolToObject(gp, "peakfinder8_fast", args->peak_search.peakfinder8_fast);
	cJSON_AddBoolToObject(gp, "half_pixel_shift", args->peak_search.half_pixel_shift);
	cJSON_AddNumberToObject(gp, "min_res_px", args->peak_search.min_res);
	cJSON_AddNumberToObject(gp, "max_res_px", args->peak_search.max_res);
	cJSON_AddNumberToObject(gp, "min_pixel_count", args->peak_search.min_pix_count);
	cJSON_AddNumberToObject(gp, "max_pixel_count", args->peak_search.max_pix_count);
	cJSON_AddNumberToObject(gp, "local_bg_radius_px", args->peak_search.local_bg_radius);
	cJSON_AddBoolToObject(gp, "use_saturated", args->peak_search.use_saturated);
	cJSON_AddBoolToObject(gp, "revalidate_hdf5", args->peak_search.revalidate);
	cJSON_AddNumberToObject(gp, "min_snr_of_biggest_pixel", args->peak_search.min_snr_biggest_pix);
	cJSON_AddNumberToObject(gp, "min_snr_of_peak_pixel", args->peak_search.min_snr_peak_pix);
	cJSON_AddNumberToObject(gp, "min_sig_adu", args->peak_search.min_sig);
	cJSON_AddNumberToObject(gp, "min_peak_over_neighbour_adu", args->peak_search.min_peak_over_neighbour);

	gp = cJSON_AddObjectToObject(harvest, "hitfinding");
	cJSON_AddNumberToObject(gp, "min_num_peaks", args->min_peaks);

	if ( args->ipriv == NULL ) {
		cJSON_AddNullToObject(harvest, "indexing");
		cJSON_AddNullToObject(harvest, "integration");
	} else {

		gp = cJSON_AddObjectToObject(harvest, "indexing");
		write_methods(gp, "methods", args->ipriv);
		write_json_cell(gp, "target_cell", args->cell);

		cJSON *tol = cJSON_AddObjectToObject(gp, "tolerances");
		cJSON_AddNumberToObject(tol, "a_percent", 100.0*args->tols[0]);
		cJSON_AddNumberToObject(tol, "b_percent", 100.0*args->tols[1]);
		cJSON_AddNumberToObject(tol, "c_percent", 100.0*args->tols[2]);
		cJSON_AddNumberToObject(tol, "alpha_deg", rad2deg(args->tols[3]));
		cJSON_AddNumberToObject(tol, "beta_deg",  rad2deg(args->tols[4]));
		cJSON_AddNumberToObject(tol, "gamma_deg", rad2deg(args->tols[5]));

		cJSON_AddBoolToObject(gp, "multi_lattice", if_multi);
		cJSON_AddBoolToObject(gp, "refine", if_refine);
		cJSON_AddBoolToObject(gp, "retry", if_retry);
		cJSON_AddBoolToObject(gp, "check_peaks", if_peaks);
		cJSON_AddBoolToObject(gp, "check_cell", if_checkcell);
		cJSON_AddNumberToObject(gp, "wavelength_estimate_m", args->wavelength_estimate);
		cJSON_AddNumberToObject(gp, "clen_estimate_m", args->clen_estimate);

		gp = cJSON_AddObjectToObject(harvest, "integration");

		char *tmp = str_integration_method(args->int_meth);
		cJSON_AddStringToObject(gp, "method", tmp);
		free(tmp);

		radii = cJSON_AddObjectToObject(gp, "radii");
		cJSON_AddNumberToObject(radii, "inner_px", args->ir_inn);
		cJSON_AddNumberToObject(radii, "middle_px",  args->ir_mid);
		cJSON_AddNumberToObject(radii, "outer_px",  args->ir_out);

		cJSON_AddNumberToObject(gp, "push_res_invm", args->push_res);
		cJSON_AddNumberToObject(gp, "fix_profile_radius_invm", nan_if_neg(args->fix_profile_r));
		cJSON_AddNumberToObject(gp, "fix_divergence_rad", nan_if_neg(args->fix_divergence));
		cJSON_AddBoolToObject(gp, "overpredict", args->overpredict);
		cJSON_AddBoolToObject(gp, "cell_parameters_only", args->cell_params_only);
	}

	char *json = cJSON_Print(harvest);
	if ( json == NULL ) {
		cJSON_Delete(harvest);
		return;
	}

	fh = fopen(filename, "w");
	if ( fh == NULL ) {
		ERROR("Unable to write parameter harvesting file.\n");
		return;
	}
	fputs(json, fh);
	fclose(fh);

	cJSON_Delete(harvest);
	free(json);

	#else /* HAVE_CJSON */
	ERROR("Cannot write harvest file - cJSON library not present\n");
	#endif /* HAVE_CJSON */
}


static void shuffle_events(struct sb_shm *sb_shared)
{
	int i;

	for ( i=1; i<sb_shared->n_events; i++ ) {
		memcpy(sb_shared->queue[i-1], sb_shared->queue[i], MAX_EV_LEN);
	}
	sb_shared->n_events--;
}


static void pin_to_cpu(int slot)
{
	#ifdef HAVE_SCHED_SETAFFINITY
	cpu_set_t c;

	CPU_ZERO(&c);
	CPU_SET(slot, &c);
	if ( sched_setaffinity(0, sizeof(cpu_set_t), &c) ) {
		fprintf(stderr, "Failed to set CPU affinity for %i\n", slot);
	}
	#endif
}


struct indexamajig_debug_data
{
	struct sb_shm *shared;
	int worker;
};


static void set_last_task_sandbox(const char *task, void *vp)
{
	struct indexamajig_debug_data *db = vp;
	assert(strlen(task) < MAX_TASK_LEN-1);
	pthread_mutex_lock(&db->shared->debug_lock);
	strcpy(db->shared->last_task[db->worker], task);
	pthread_mutex_unlock(&db->shared->debug_lock);
}


static int notify_alive_sandbox(void *vp)
{
	struct indexamajig_debug_data *db = vp;
	pthread_mutex_lock(&db->shared->debug_lock);
	db->shared->pings[db->worker]++;
	pthread_mutex_unlock(&db->shared->debug_lock);
	return 0;
}


static int run_work(struct indexamajig_arguments *args)
{
	int allDone = 0;
	struct im_zmq *zmqstuff = NULL;
	struct im_asapo *asapostuff = NULL;
	Mille *mille;
	ImageDataArrays *ida;
	size_t ll;
	char *tmp;
	struct stat s;
	Stream *st;
	int shm_fd;
	sem_t *queue_sem;
	IndexingFlags flags = 0;
	struct pf8_private_data *pf8_data = NULL;
	struct indexamajig_debug_data debugdata;
	struct sb_shm *shared;

	if ( args->cpu_pin ) pin_to_cpu(args->worker_id);

	ll = 64 + strlen(args->worker_tmpdir);
	tmp = malloc(ll);
	if ( tmp == NULL ) {
		ERROR("Failed to allocate temporary dir\n");
		return 1;
	}

	snprintf(tmp, ll, "%s/worker.%i", args->worker_tmpdir, args->worker_id);

	if ( stat(tmp, &s) == -1 ) {

		int r;

		if ( errno != ENOENT ) {
			ERROR("Failed to stat temporary folder.\n");
			exit(1);
		}

		r = mkdir(tmp, S_IRWXU);
		if ( r ) {
			ERROR("Failed to create temporary folder: %s\n",
			strerror(errno));
			exit(1);
		}
	}

	st = stream_open_fd_for_write(args->fd_stream, args->iargs.dtempl);

	if ( args->profile ) {
		profile_init();
	}

	/* Load unit cell (if given) */
	if ( args->cellfile != NULL ) {
		args->iargs.cell = load_cell_from_file(args->cellfile);
		if ( args->iargs.cell == NULL ) {
			ERROR("Couldn't read unit cell (from %s)\n", args->cellfile);
			return 1;
		}
	} else {
		args->iargs.cell = NULL;
	}

	if ( args->if_checkcell ) {
		flags |= INDEXING_CHECK_CELL;
	}
	if ( args->if_refine ) {
		flags |= INDEXING_REFINE;
	}
	if ( args->if_peaks ) {
		flags |= INDEXING_CHECK_PEAKS;
	}
	if ( args->if_multi ) {
		flags |= INDEXING_MULTI;
	}
	if ( args->if_retry ) {
		flags |= INDEXING_RETRY;
	}

	args->iargs.ipriv = setup_indexing(args->indm_str,
	                                   args->iargs.cell,
	                                   args->iargs.tols,
	                                   flags,
	                                   args->iargs.wavelength_estimate,
	                                   args->iargs.clen_estimate,
	                                   args->iargs.n_threads,
	                                   *args->taketwo_opts_ptr,
	                                   *args->xgandalf_opts_ptr,
	                                   *args->ffbidx_opts_ptr,
	                                   *args->pinkindexer_opts_ptr,
	                                   *args->felix_opts_ptr,
	                                   *args->fromfile_opts_ptr,
	                                   *args->smallcell_opts_ptr,
	                                   *args->asdf_opts_ptr);

	if ( args->iargs.ipriv == NULL ) {
		ERROR("Failed to set up indexing system\n");
		return 1;
	}

	if ( args->worker_id == 0 ) {

		print_indexing_info(args->iargs.ipriv);

		if ( (args->harvest_file != NULL) && (args->serial_start <= 1) ) {
			write_harvest_file(&args->iargs, args->harvest_file,
			                   args->if_multi, args->if_refine, args->if_retry,
			                   args->if_peaks, args->if_checkcell);
		}

	}

	if ( args->iargs.peak_search.method == PEAK_PEAKFINDER8 ) {
		struct detgeom *dg;
		dg = data_template_get_2d_detgeom_if_possible(args->iargs.dtempl);
		if ( dg == NULL ) {
			ERROR("WARNING: Detector geometry is not static.  "
			      "Peak search will be slower than optimal.\n");
		}
		pf8_data = prepare_peakfinder8(dg, args->iargs.peak_search.peakfinder8_fast);
		args->iargs.pf_private = pf8_data;
		detgeom_free(dg);
	}

	/* Set up SHM */
	shm_fd = shm_open(args->shm_name, O_RDWR, 0);
	if ( shm_fd == -1 ) {
		ERROR("SHM setup failed: %s\n", strerror(errno));
		return 1;
	}
	shared = mmap(NULL, sizeof(struct sb_shm), PROT_READ | PROT_WRITE,
			MAP_SHARED, shm_fd, 0);
	if ( shared == MAP_FAILED ) {
		ERROR("SHM setup failed: %s\n", strerror(errno));
		return 1;
	}

	queue_sem = sem_open(args->queue_sem, 0);
	if ( queue_sem == SEM_FAILED ) {
		ERROR("Failed to open semaphore: %s\n", strerror(errno));
		return 1;
	}

	/* Connect via ZMQ */
	if ( args->zmq_params.addr != NULL ) {
		zmqstuff = im_zmq_connect(&args->zmq_params);
		if ( zmqstuff == NULL ) {
			ERROR("ZMQ setup failed.\n");
			return 1;
		}
	}

	if ( args->asapo_params.endpoint != NULL ) {
		asapostuff = im_asapo_connect(&args->asapo_params);
		if ( asapostuff == NULL ) {
			ERROR("ASAP::O setup failed.\n");
			pthread_mutex_lock(&shared->totals_lock);
			shared->should_shutdown = 1;
			pthread_mutex_unlock(&shared->totals_lock);
			return 1;
		}
	}

	debugdata.shared = shared;
	debugdata.worker = args->worker_id;
	set_debug_funcs(set_last_task_sandbox, notify_alive_sandbox, &debugdata);

	mille = crystfel_mille_new_fd(args->fd_mille);

	ida = image_data_arrays_new();

	while ( !allDone ) {

		struct pattern_args pargs;
		int ser;
		char *line;
		size_t len;
		int i;
		char *event_str = NULL;
		char *ser_str = NULL;
		int ok = 1;
		int should_shutdown;

		/* Wait until an event is ready */
		notify_alive();
		set_last_task("wait_event");
		profile_start("wait-queue-semaphore");
		if ( sem_wait(queue_sem) != 0 ) {
			ERROR("Failed to wait on queue semaphore: %s\n",
			      strerror(errno));
		}
		profile_end("wait-queue-semaphore");

		/* Get the event from the queue */
		set_last_task("read_queue");
		pthread_mutex_lock(&shared->totals_lock);
		should_shutdown = shared->should_shutdown;
		pthread_mutex_unlock(&shared->totals_lock);
		pthread_mutex_lock(&shared->queue_lock);
		if ( ((shared->n_events==0) && (shared->no_more))
		   || should_shutdown )
		{
			/* Queue is empty and no more are coming,
			 * or another process has initiated a shutdown.
			 * Either way, it's time to get out of here. */
			pthread_mutex_unlock(&shared->queue_lock);
			allDone = 1;
			continue;
		}
		if ( shared->n_events == 0 ) {
			ERROR("Got the semaphore, but no events in queue!\n");
			ERROR("no_more = %i\n", shared->no_more);
			pthread_mutex_unlock(&shared->queue_lock);
			allDone = 1;
			continue;
		}

		line = strdup(shared->queue[0]);

		len = strlen(line);
		assert(len > 1);
		for ( i=len-1; i>0; i-- ) {
			if ( line[i] == ' ' ) {
				line[i] = '\0';
				ser_str = &line[i+1];
				break;
			}
		}
		len = strlen(line);
		assert(len > 1);
		for ( i=len-1; i>0; i-- ) {
			if ( line[i] == ' ' ) {
				line[i] = '\0';
				event_str = &line[i+1];
				break;
			}
		}
		if ( (ser_str != NULL) && (event_str != NULL) ) {
			if ( sscanf(ser_str, "%i", &ser) != 1 ) {
				STATUS("Invalid serial number '%s'\n",
				       ser_str);
				ok = 0;
			}
		}
		if ( !ok ) {
			STATUS("Invalid event string '%s'\n",
			       shared->queue[0]);
			ok = 0;
		}

		pthread_mutex_lock(&shared->debug_lock);
		memcpy(shared->last_ev[args->worker_id], shared->queue[0],
		       MAX_EV_LEN);
		pthread_mutex_unlock(&shared->debug_lock);

		shuffle_events(shared);
		pthread_mutex_unlock(&shared->queue_lock);

		if ( !ok ) continue;

		pargs.filename = strdup(line);
		pargs.event = safe_strdup(event_str);

		free(line);
		ok = 0;

		/* Default values */
		pargs.zmq_data = NULL;
		pargs.zmq_data_size = 0;
		pargs.asapo_data = NULL;
		pargs.asapo_data_size = 0;
		pargs.asapo_meta = NULL;

		if ( args->zmq_params.addr != NULL ) {

			profile_start("zmq-fetch");
			set_last_task("ZMQ fetch");
			pargs.zmq_data = im_zmq_fetch(zmqstuff,
			                              &pargs.zmq_data_size);
			profile_end("zmq-fetch");

			if ( (pargs.zmq_data != NULL)
			  && (pargs.zmq_data_size > 15) ) ok = 1;

			/* The filename/event, which will be 'fake' values in
			 * this case, still came via the event queue.  More
			 * importantly, the event queue gave us a unique
			 * serial number for this image. */

		} else if ( args->asapo_params.endpoint != NULL ) {

			char *filename;
			char *event;
			int finished = 0;
			int asapo_message_id;

			profile_start("asapo-fetch");
			set_last_task("ASAPO fetch");
			pargs.asapo_data = im_asapo_fetch(asapostuff,
			                                  &pargs.asapo_data_size,
			                                  &pargs.asapo_meta,
			                                  &filename,
			                                  &event,
			                                  &finished,
			                                  &asapo_message_id);
			profile_end("asapo-fetch");
			if ( pargs.asapo_data != NULL ) {
				ok = 1;

				/* ASAP::O provides a meaningful filename, which
				 * replaces the placeholder. */
				free(pargs.filename);
				free(pargs.event);
				pargs.filename = filename;
				pargs.event = event;
				pthread_mutex_lock(&shared->queue_lock);
				shared->end_of_stream[args->worker_id] = 0;
				pthread_mutex_unlock(&shared->queue_lock);

				/* We will also use ASAP::O's serial number
				 * instead of our own. */
				ser = asapo_message_id;
			} else {
				if ( finished ) {
					pthread_mutex_lock(&shared->queue_lock);
					shared->end_of_stream[args->worker_id] = 1;
					pthread_mutex_unlock(&shared->queue_lock);
				}
			}

		} else {
			ok = 1;
		}

		if ( ok ) {
			pthread_mutex_lock(&shared->debug_lock);
			shared->time_last_start[args->worker_id] = get_monotonic_seconds();
			pthread_mutex_unlock(&shared->debug_lock);
			profile_start("process-image");
			process_image(&args->iargs, &pargs, st, args->worker_id,
			              args->worker_tmpdir, ser,
			              shared, asapostuff, mille, ida);
			profile_end("process-image");

			if ( asapostuff != NULL ) {
				im_asapo_finalise(asapostuff, ser);
			}

		}

		/* NB pargs.zmq_data, pargs.asapo_data and  pargs.asapo_meta
		 * will be copied into the image structure, so
		 * that it can be queried for "header" values etc.  They will
		 * eventually be freed by image_free() under process_image(). */

		if ( args->profile ) {
			profile_print_and_reset(args->worker_id);
		}

		free(pargs.filename);
		free(pargs.event);
	}

	stream_close(st);
	free(tmp);
	munmap(shared, sizeof(struct sb_shm));
	sem_close(queue_sem);

	image_data_arrays_free(ida);
	crystfel_mille_free(mille);

	/* These are both no-ops if argument is NULL */
	im_zmq_shutdown(zmqstuff);
	im_asapo_shutdown(asapostuff);

	data_template_free(args->iargs.dtempl);
	if ( pf8_data != NULL ) free_pf8_private_data(pf8_data);
	cleanup_indexing(args->iargs.ipriv);
	cell_free(args->iargs.cell);
	cleanup_indexamajig_args(args);

#ifdef HAVE_FFTW
	fftw_cleanup();
#endif

	return 0;
}


int main(int argc, char *argv[])
{
	FILE *fh = NULL;
	Stream *st;
	struct indexamajig_arguments *args;
	char *tmpdir;  /* e.g. /tmp/indexamajig.12345 */
	char *rn;  /* e.g. /home/taw/indexing */
	int r;
	int timeout = 240;
	double wl_from_dt;
	double clen_from_dt;
	int err = 0;
	char *probed_methods = NULL;
	size_t mille_fn_len;
	char *mille_filename;
	FILE *mille_fh;

	args = parse_indexamajig_args(argc, argv);
	if ( args == NULL ) return 1;

	/* Load data template (new API) */
	args->iargs.dtempl = data_template_new_from_file(args->geom_filename);
	if ( args->iargs.dtempl == NULL ) {
		ERROR("Failed to read detector geometry from '%s'\n",
		      args->geom_filename);
		return 1;
	}

	/* Add any headers we need to copy */
	for ( r=0; r<args->n_copy_headers; r++ ) {
		data_template_add_copy_header(args->iargs.dtempl,
		                              args->copy_headers[r]);
	}


	wl_from_dt = data_template_get_wavelength_if_possible(args->iargs.dtempl);
	if ( !isnan(wl_from_dt) ) {
		if ( !isnan(args->iargs.wavelength_estimate) && !args->worker) {
			ERROR("WARNING: Ignoring your value for "
			      "--wavelength-estimate because the geometry file "
			      "already contains a static value.\n");
		}
		args->iargs.wavelength_estimate = wl_from_dt;
	}

	clen_from_dt = data_template_get_clen_if_possible(args->iargs.dtempl);
	if ( !isnan(clen_from_dt) ) {
		if ( !isnan(args->iargs.clen_estimate) && !args->worker) {
			ERROR("WARNING: Ignoring your value for "
			      "--camera-length-estimate because the geometry file "
			      "already contains a static value.\n");
		}
		args->iargs.clen_estimate = clen_from_dt;
	}

	/* If no integration radii were given, apply the defaults */
	if ( args->iargs.ir_inn < 0 ) {
		args->iargs.ir_inn = 4.0;
		args->iargs.ir_mid = 5.0;
		args->iargs.ir_out = 7.0;
	}

	/* If no peak radii were given, copy the integration radii */
	if ( args->iargs.peak_search.pk_inn < 0.0 ) {
		args->iargs.peak_search.pk_inn = args->iargs.ir_inn;
		args->iargs.peak_search.pk_mid = args->iargs.ir_mid;
		args->iargs.peak_search.pk_out = args->iargs.ir_out;
	}

	if ( args->worker ) {
		/* I am a worker process */
		return run_work(args);
	} /* else I am the main process */

	/* Check for minimal information */
	if ( (args->filename == NULL)
	  && (args->zmq_params.addr == NULL)
	  && (args->asapo_params.endpoint == NULL) ) {
		ERROR("You need to provide the input filename (use -i)\n");
		return 1;
	}
	if ( args->geom_filename == NULL ) {
		ERROR("You need to specify the geometry filename (use -g)\n");
		return 1;
	}
	if ( args->outfile == NULL ) {
		ERROR("You need to specify the output filename (use -o)\n");
		return 1;
	}

	if ( (args->filename != NULL) && (args->zmq_params.addr != NULL) ) {
		ERROR("The options --input and --zmq-input are mutually "
		      "exclusive.\n");
		return 1;
	}

	if ( (args->filename != NULL) && (args->asapo_params.endpoint != NULL) ) {
		ERROR("The options --input and --asapo-endpoint are mutually "
		      "exclusive.\n");
		return 1;
	}

	if ( (args->asapo_params.endpoint != NULL) && (args->zmq_params.addr != NULL) ) {
		ERROR("The options --asapo-endpoint and --zmq-input are mutually "
		      "exclusive.\n");
		return 1;
	}

	if ( (args->zmq_params.request != NULL) && (args->zmq_params.n_subscriptions > 0) ) {
		ERROR("The options --zmq-request and --zmq-subscribe are "
		      "mutually exclusive.\n");
		return 1;
	}

	if ( (args->filename != NULL)
	  && (strcmp(args->filename, "-") != 0)
	  && is_hdf5_file(args->filename, &err) )
	{
		ERROR("Your input file appears to be an HDF5 file.\n");
		ERROR("The input file should be a list of data files, not the "
		      "data file itself.\n");
		ERROR("If you have only one input file, try the following:\n");
		ERROR("  echo %s > files.lst\n", args->filename);
		ERROR("  indexamajig -i files.lst -o %s -g %s ...\n",
		      args->outfile, args->geom_filename);
		return 1;
	} else if ( err ) {
		ERROR("Couldn't open '%s'\n", args->filename);
		return 1;
	}

	/* Open input */
	if ( args->filename != NULL ) {
		if ( strcmp(args->filename, "-") == 0 ) {
			fh = stdin;
		} else {
			fh = fopen(args->filename, "r");
		}
		if ( fh == NULL ) {
			ERROR("Failed to open input file '%s'\n", args->filename);
			return 1;
		}
	}

	/* Check prefix (if given) */
	if ( args->check_prefix ) {
		args->prefix = check_prefix(args->prefix);
	}

	/* Check number of processes */
	if ( args->n_proc == 0 ) {
		ERROR("Invalid number of processes.\n");
		return 1;
	}

	/* Load unit cell (if given) */
	if ( args->cellfile != NULL ) {
		args->iargs.cell = load_cell_from_file(args->cellfile);
		if ( args->iargs.cell == NULL ) {
			ERROR("Couldn't read unit cell (from %s)\n", args->cellfile);
			return 1;
		}
	} else {
		args->iargs.cell = NULL;
	}

	if ( args->iargs.cell != NULL ) {
		STATUS("This is what I understood your unit cell to be:\n");
		cell_print(args->iargs.cell);
	} else {
		STATUS("No reference unit cell provided.\n");
	}

	tmpdir = create_tempdir(args->temp_location);
	if ( tmpdir == NULL ) return 1;

	/* Change into temporary folder, temporarily, to control the crap
	 * dropped by indexing programs during setup */
	rn = getcwd(NULL, 0);
	r = chdir(tmpdir);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		return 1;
	}

	/* Auto-detect indexing methods if 'requested' */
	if ( args->indm_str == NULL ) {

		STATUS("No indexing methods specified.  I will try to ");
		STATUS("automatically detect the available methods.\n");
		STATUS("To disable auto-detection of indexing methods, specify ");
		STATUS("which methods to use with --indexing=<methods>.\n");
		STATUS("Use --indexing=none to disable indexing and integration.\n");

		probed_methods = detect_indexing_methods(args->iargs.cell);
		if ( probed_methods == NULL ) {
			ERROR("No indexing method specified, and no usable indexing ");
			ERROR("methods auto-detected.\n");
			ERROR("Install some indexing programs (mosflm,dirax etc), or ");
			ERROR("try again with --indexing=none.\n");
			return 1;
		} else {
			STATUS("Auto-determined indexing methods: %s\n", probed_methods);
		}

	} else if ( strcmp(args->indm_str, "none") == 0 ) {

		STATUS("Indexing/integration disabled.\n");
		if ( args->iargs.cell != NULL ) {
			STATUS("Ignoring your unit cell.\n");
		}
		args->iargs.ipriv = NULL;

	} else {

		if ( strstr(args->indm_str, "pinkindexer") != NULL ) {
			/* Extend timeout if using pinkIndexer */
			timeout = 3000;
		}

	}

	/* Change back to where we were before.  Sandbox code will create
	 * worker subdirs inside the temporary folder, and process_image will
	 * change into them. */
	r = chdir(rn);
	if ( r ) {
		ERROR("Failed to chdir: %s\n", strerror(errno));
		return 1;
	}
	free(rn);

	/* Open output stream */
	st = stream_open_for_write(args->outfile, args->iargs.dtempl);
	if ( st == NULL ) {
		ERROR("Failed to open stream '%s'\n", args->outfile);
		return 1;
	}

	/* Write audit info */
	stream_write_commandline_args(st, argc, argv);
	stream_write_geometry_file(st, args->geom_filename);
	stream_write_target_cell(st, args->iargs.cell);
	stream_write_indexing_methods(st, args->indm_str);

	r = mkdir(args->milledir, S_IRWXU);
	if ( r ) {
		if ( errno != EEXIST ) {
			ERROR("Failed to create folder for Millepede data: %s\n",
			      strerror(errno));
			return 1;
		}
	}
	mille_fn_len = strlen(args->milledir) +strlen(args->millefile)+2;
	mille_filename = malloc(mille_fn_len);
	if ( mille_filename == NULL ) {
		ERROR("Failed to generate Millepede filename\n");
		return 1;
	}
	snprintf(mille_filename, mille_fn_len, "%s/%s", args->milledir, args->millefile);
	mille_fh = fopen(mille_filename, "wb");
	if ( mille_fh == NULL ) {
		ERROR("Failed to open Millepede file: %s\n", strerror(errno));
		return 1;
	}
	free(mille_filename);

	r = create_sandbox(&args->iargs, args->n_proc, args->prefix, args->basename,
	                   fh, st, tmpdir, args->serial_start,
	                   &args->zmq_params, &args->asapo_params,
	                   timeout, args->profile, args->cpu_pin,
	                   args->no_data_timeout, argc, argv,
			   probed_methods, mille_fh);

	fclose(mille_fh);
	free(tmpdir);
	free(probed_methods);
	data_template_free(args->iargs.dtempl);
	cell_free(args->iargs.cell);
	stream_close(st);
	cleanup_indexamajig_args(args);

	return r;
}
