/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2012 Thomas White <taw@physics.org>
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
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

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <hdf5.h>
#include <gsl/gsl_errno.h>
#include <pthread.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <signal.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "utils.h"
#include "hdf5-file.h"
#include "index.h"
#include "peaks.h"
#include "detector.h"
#include "filters.h"
#include "thread-pool.h"
#include "beam-parameters.h"
#include "geometry.h"
#include "stream.h"
#include "reflist-utils.h"

#include "im-sandbox.h"


/* Write statistics at APPROXIMATELY this interval */
#define STATS_EVERY_N_SECONDS (5)


struct sandbox
{
	pthread_mutex_t lock;

	int n_indexable;
	int n_processed;
	int n_indexable_last_stats;
	int n_processed_last_stats;
	int t_last_stats;

	struct index_args *iargs;

	int n_proc;
	pid_t *pids;
	FILE *ofh;
	FILE **fhs;

	int *running;
	FILE **result_fhs;
	int *filename_pipes;
	int *stream_pipe_read;
	int *stream_pipe_write;
	char **last_filename;
};


/* Horrible global variable for signal handler */
struct sandbox *sb;


static char *get_pattern(FILE *fh, char **use_this_one_instead,
                         int config_basename, const char *prefix)
{
	char *line;
	char *filename;

	do {

		/* Get the next filename */
		if ( *use_this_one_instead != NULL ) {

			line = *use_this_one_instead;
			*use_this_one_instead = NULL;

		} else {

			char *rval;

			line = malloc(1024*sizeof(char));
			rval = fgets(line, 1023, fh);
			if ( rval == NULL ) {
				free(line);
				return NULL;
			}

		}

		chomp(line);

	} while ( strlen(line) == 0 );

	if ( config_basename ) {
		char *tmp;
		tmp = safe_basename(line);
		free(line);
		line = tmp;
	}

	filename = malloc(strlen(prefix)+strlen(line)+1);

	snprintf(filename, 1023, "%s%s", prefix, line);

	free(line);

	return filename;
}


static void process_image(const struct index_args *iargs,
                          struct pattern_args *pargs, FILE *ofh,
                          int cookie)
{
	float *data_for_measurement;
	size_t data_size;
	UnitCell *cell = iargs->cell;
	int config_cmfilter = iargs->config_cmfilter;
	int config_noisefilter = iargs->config_noisefilter;
	int config_verbose = iargs->config_verbose;
	IndexingMethod *indm = iargs->indm;
	struct beam_params *beam = iargs->beam;
	int r, check;
	struct hdfile *hdfile;
	struct image image;

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.indexed_cell = NULL;
	image.det = copy_geom(iargs->det);
	image.copyme = iargs->copyme;
	image.reflections = NULL;
	image.id = cookie;
	image.filename = pargs->filename;
	image.beam = beam;

	hdfile = hdfile_open(image.filename);
	if ( hdfile == NULL ) return;

	r = hdfile_set_first_image(hdfile, "/");
	if ( r ) {
		ERROR("Couldn't select first path\n");
		hdfile_close(hdfile);
		return;
	}

	check = hdf5_read(hdfile, &image, 1);
	if ( check ) {
		hdfile_close(hdfile);
		return;
	}

	if ( (image.width != image.det->max_fs + 1 )
	  || (image.height != image.det->max_ss + 1))
	{
		ERROR("Image size doesn't match geometry size"
			" - rejecting image.\n");
		ERROR("Image size: %i,%i.  Geometry size: %i,%i\n",
		image.width, image.height,
		image.det->max_fs + 1, image.det->max_ss + 1);
		hdfile_close(hdfile);
		free_detector_geometry(image.det);
		return;
	}

	if ( image.lambda < 0.0 ) {
		if ( beam != NULL ) {
			ERROR("Using nominal photon energy of %.2f eV\n",
			beam->photon_energy);
			image.lambda = ph_en_to_lambda(
			eV_to_J(beam->photon_energy));
		} else {
			ERROR("No wavelength in file, so you need to give "
				"a beam parameters file with -b.\n");
			hdfile_close(hdfile);
			free_detector_geometry(image.det);
			return;
		}
	}
	fill_in_values(image.det, hdfile);

	if ( config_cmfilter ) {
		filter_cm(&image);
	}

	/* Take snapshot of image after CM subtraction but before
	 * the aggressive noise filter. */
	data_size = image.width * image.height * sizeof(float);
	data_for_measurement = malloc(data_size);

	if ( config_noisefilter ) {
		filter_noise(&image, data_for_measurement);
	} else {
		memcpy(data_for_measurement, image.data, data_size);
	}

	switch ( iargs->peaks ) {

		case PEAK_HDF5:
		// Get peaks from HDF5
		if (get_peaks(&image, hdfile,
			iargs->hdf5_peak_path)) {
			ERROR("Failed to get peaks from HDF5 file.\n");
		}
		break;

		case PEAK_ZAEF:
		search_peaks(&image, iargs->threshold,
		             iargs->min_gradient, iargs->min_snr,
		             iargs->ir_inn, iargs->ir_mid, iargs->ir_out);
		break;

	}

	/* Get rid of noise-filtered version at this point
	 * - it was strictly for the purposes of peak detection. */
	free(image.data);
	image.data = data_for_measurement;

	/* Calculate orientation matrix (by magic) */
	image.div = beam->divergence;
	image.bw = beam->bandwidth;
	image.profile_radius = 0.0001e9;

	index_pattern(&image, cell, indm, iargs->cellr,
	              config_verbose, iargs->ipriv,
	              iargs->config_insane, iargs->tols);

	if ( image.indexed_cell != NULL ) {
		pargs->indexable = 1;
		image.reflections = find_intersections(&image,
				image.indexed_cell);
		if (image.reflections != NULL) {
			integrate_reflections(&image,
			                      iargs->config_closer,
			                      iargs->config_bgsub,
			                      iargs->min_int_snr,
			                      iargs->ir_inn,
			                      iargs->ir_mid,
			                      iargs->ir_out);
		}
	} else {
		image.reflections = NULL;
	}

	write_chunk(ofh, &image, hdfile, iargs->stream_flags);
	fprintf(ofh, "END\n");
	fflush(ofh);

	/* Only free cell if found */
	cell_free(image.indexed_cell);

	reflist_free(image.reflections);
	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	image_feature_list_free(image.features);
	hdfile_close(hdfile);
	free_detector_geometry(image.det);
}


static void run_work(const struct index_args *iargs,
                     int filename_pipe, int results_pipe, FILE *ofh,
                     int cookie)
{
	int allDone = 0;
	FILE *fh;
	int w;

	fh = fdopen(filename_pipe, "r");
	if ( fh == NULL ) {
		ERROR("Failed to fdopen() the filename pipe!\n");
		return;
	}

	w = write(results_pipe, "\n", 1);
	if ( w < 0 ) {
		ERROR("Failed to send request for first filename.\n");
	}

	while ( !allDone ) {

		struct pattern_args pargs;
		int  c;
		char *line;
		char *rval;
		char buf[1024];

		line = malloc(1024*sizeof(char));
		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) {
			free(line);
			if ( feof(fh) ) {
				allDone = 1;
				STATUS("Exiting!\n");
				continue;
			} else {
				ERROR("Read error!\n");
				break;
			}
		}

		chomp(line);

		if ( strlen(line) == 0 ) {

			allDone = 1;

		} else {

			pargs.filename = line;
			pargs.indexable = 0;

			process_image(iargs, &pargs, ofh, cookie);

			/* Request another image */
			c = sprintf(buf, "%i\n", pargs.indexable);
			w = write(results_pipe, buf, c);
			if ( w < 0 ) {
				ERROR("write P0\n");
			}

		}

		free(line);

	}

	cleanup_indexing(iargs->ipriv);
	free(iargs->indm);
	free(iargs->ipriv);
	free_detector_geometry(iargs->det);
	free(iargs->beam);
	free(iargs->element);
	free(iargs->hdf5_peak_path);
	free_copy_hdf5_field_list(iargs->copyme);
	cell_free(iargs->cell);
	fclose(fh);
}


#ifdef HAVE_CLOCK_GETTIME

static time_t get_monotonic_seconds()
{
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);
	return tp.tv_sec;
}

#else

/* Fallback version of the above.  The time according to gettimeofday() is not
 * monotonic, so measuring intervals based on it will screw up if there's a
 * timezone change (e.g. daylight savings) while the program is running. */
static time_t get_monotonic_seconds()
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec;
}

#endif

static int pump_chunk(FILE *fh, FILE *ofh)
{
	int chunk_started = 0;
	int chunk_finished = 0;

	do {

		char line[1024];
		char *rval;

		rval = fgets(line, 1024, fh);
		if ( rval == NULL ) {

			if ( feof(fh) ) {
				/* Process died */
				if ( chunk_started ) {
					ERROR("EOF during chunk!\n");
					fprintf(ofh, "Chunk is unfinished!\n");
				}
				return 1;
			} else {
				ERROR("fgets() failed: %s\n", strerror(errno));
			}
			chunk_finished = 1;
			continue;

		}

		if ( strcmp(line, "END\n") == 0 ) {
			chunk_finished = 1;
		} else {
			chunk_started = 1;
			fprintf(ofh, "%s", line);
		}

	} while ( !chunk_finished );

	return 0;
}


static void *run_reader(void *sbv)
{
	struct sandbox *sb = sbv;
	int done = 0;

	while ( !done ) {

		int r, i;
		struct timeval tv;
		fd_set fds;
		int fdmax;

		tv.tv_sec = 5;
		tv.tv_usec = 0;

		FD_ZERO(&fds);
		fdmax = 0;
		pthread_mutex_lock(&sb->lock);
		for ( i=0; i<sb->n_proc; i++ ) {

			int fd;

			if ( !sb->running[i] ) continue;

			fd = sb->stream_pipe_read[i];

			FD_SET(fd, &fds);
			if ( fd > fdmax ) fdmax = fd;

		}

		pthread_mutex_unlock(&sb->lock);

		r = select(fdmax+1, &fds, NULL, NULL, &tv);

		if ( r == -1 ) {
			if ( errno != EINTR ) {
				ERROR("select() failed: %s\n", strerror(errno));
			} /* Otherwise no big deal */
			continue;
		}

		if ( r == 0 ) continue; /* Nothing this time.  Try again */

		pthread_mutex_lock(&sb->lock);
		for ( i=0; i<sb->n_proc; i++ ) {

			if ( !sb->running[i] ) continue;

			if ( !FD_ISSET(sb->stream_pipe_read[i], &fds) ) continue;

			if ( pump_chunk(sb->fhs[i], sb->ofh) ) {
				sb->running[i] = 0;
			}

		}

		done = 1;
		for ( i=0; i<sb->n_proc; i++ ) {
			if ( sb->running[i] ) done = 0;
		}
		pthread_mutex_unlock(&sb->lock);

	}

	return NULL;
}


static void start_worker_process(struct sandbox *sb, int slot)
{
	pid_t p;
	int filename_pipe[2];
	int result_pipe[2];

	if ( pipe(filename_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	if ( pipe(result_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	pthread_mutex_lock(&sb->lock);
	p = fork();
	if ( p == -1 ) {
		ERROR("fork() failed!\n");
		pthread_mutex_unlock(&sb->lock);
		return;
	}

	if ( p == 0 ) {

		FILE *sfh;
		int j;
		struct sigaction sa;
		int r;

		/* FIXME: Is lock inherited? */
		pthread_mutex_unlock(&sb->lock);

		/* First, disconnect the signal handler */
		sa.sa_flags = 0;
		sigemptyset(&sa.sa_mask);
		sa.sa_handler = SIG_DFL;
		r = sigaction(SIGCHLD, &sa, NULL);
		if ( r == -1 ) {
			ERROR("Failed to set signal handler!\n");
			return;
		}

		/* Free resources which will not be needed by worker */
		for ( j=0; j<sb->n_proc; j++ ) {
			if ( (j != slot) && (sb->running[j]) ) {
				close(sb->stream_pipe_write[j]);
			}
		}
		for ( j=0; j<sb->n_proc; j++ ) {
			if ( (j != slot) && (sb->running[j]) ) {
				fclose(sb->result_fhs[j]);
				close(sb->filename_pipes[j]);
			}
		}
		free(sb->filename_pipes);
		free(sb->result_fhs);
		free(sb->pids);
		/* Also prefix, use_this_one_instead and fh */

		/* Child process gets the 'read' end of the filename
		 * pipe, and the 'write' end of the result pipe. */
		close(filename_pipe[1]);
		close(result_pipe[0]);

		sfh = fdopen(sb->stream_pipe_write[slot], "w");
		run_work(sb->iargs, filename_pipe[0], result_pipe[1],
		         sfh, slot);
		fclose(sfh);

		//close(filename_pipe[0]);
		close(result_pipe[1]);

		exit(0);

	}

	/* Parent process gets the 'write' end of the filename pipe
	 * and the 'read' end of the result pipe. */
	sb->pids[slot] = p;
	sb->running[slot] = 1;
	close(filename_pipe[0]);
	close(result_pipe[1]);
	sb->filename_pipes[slot] = filename_pipe[1];
	sb->fhs[slot] = fdopen(sb->stream_pipe_read[slot], "r");
	if ( sb->fhs[slot] == NULL ) {
		ERROR("Couldn't fdopen() stream!\n");
		pthread_mutex_unlock(&sb->lock);
		return;
	}

	sb->result_fhs[slot] = fdopen(result_pipe[0], "r");
	if ( sb->result_fhs[slot] == NULL ) {
		ERROR("fdopen() failed.\n");
		pthread_mutex_unlock(&sb->lock);
		return;
	}

	pthread_mutex_unlock(&sb->lock);
}


static void signal_handler(int sig, siginfo_t *si, void *uc_v)
{
	int i, found;

	if ( si->si_signo != SIGCHLD ) {
		ERROR("Unhandled signal %i?\n", si->si_signo);
		return;
	}

	found = 0;
	pthread_mutex_lock(&sb->lock);
	for ( i=0; i<sb->n_proc; i++ ) {
		if ( (sb->running[i]) && (sb->pids[i] == si->si_pid) ) {
			found = 1;
			break;
		}
	}
	pthread_mutex_unlock(&sb->lock);

	if ( !found ) {
		ERROR("SIGCHLD from unknown child %i?\n", si->si_pid);
		return;
	}

	if ( (si->si_code == CLD_TRAPPED) || (si->si_code == CLD_STOPPED)
	  || (si->si_code == CLD_CONTINUED) ) return;

	if ( si->si_code == CLD_EXITED )
	{
		pthread_mutex_lock(&sb->lock);
		sb->running[i] = 0;
		pthread_mutex_unlock(&sb->lock);
		STATUS("Worker process %i exited normally.\n", i);
		return;
	}

	if ( (si->si_code != CLD_DUMPED) && (si->si_code != CLD_KILLED) ) {
		ERROR("Unhandled si_code %i (worker process %i).\n",
		      si->si_code, i);
		return;
	}

	ERROR("Worker process %i exited abnormally!\n", i);
	ERROR("  -> Signal %i, last filename %s.\n",
	      si->si_signo, sb->last_filename[i]);

	start_worker_process(sb, i);
}


void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                    int config_basename, FILE *fh, char *use_this_one_instead,
                    FILE *ofh)
{
	int i;
	int allDone;
	struct sigaction sa;
	int r;
	pthread_t reader_thread;

	sb = calloc(1, sizeof(struct sandbox));
	if ( sb == NULL ) {
		ERROR("Couldn't allocate memory for sandbox.\n");
		return;
	}

	sb->n_indexable = 0;
	sb->n_processed = 0;
	sb->n_indexable_last_stats = 0;
	sb->n_processed_last_stats = 0;
	sb->t_last_stats = get_monotonic_seconds();
	sb->n_proc = n_proc;
	sb->ofh = ofh;
	sb->iargs = iargs;

	pthread_mutex_init(&sb->lock, NULL);

	sb->stream_pipe_read = calloc(n_proc, sizeof(int));
	sb->stream_pipe_write = calloc(n_proc, sizeof(int));
	if ( sb->stream_pipe_read == NULL ) {
		ERROR("Couldn't allocate memory for pipes.\n");
		return;
	}
	if ( sb->stream_pipe_write == NULL ) {
		ERROR("Couldn't allocate memory for pipes.\n");
		return;
	}

	for ( i=0; i<n_proc; i++ ) {

		int stream_pipe[2];

		if ( pipe(stream_pipe) == - 1 ) {
			ERROR("pipe() failed!\n");
			return;
		}

		sb->stream_pipe_read[i] = stream_pipe[0];
		sb->stream_pipe_write[i] = stream_pipe[1];

	}

	pthread_mutex_lock(&sb->lock);
	sb->filename_pipes = calloc(n_proc, sizeof(int));
	sb->result_fhs = calloc(n_proc, sizeof(FILE *));
	sb->pids = calloc(n_proc, sizeof(pid_t));
	sb->running = calloc(n_proc, sizeof(int));
	sb->fhs = calloc(sb->n_proc, sizeof(FILE *));
	if ( sb->filename_pipes == NULL ) {
		ERROR("Couldn't allocate memory for pipes.\n");
		return;
	}
	if ( sb->result_fhs == NULL ) {
		ERROR("Couldn't allocate memory for pipe file handles.\n");
		return;
	}
	if ( sb->pids == NULL ) {
		ERROR("Couldn't allocate memory for PIDs.\n");
		return;
	}
	if ( sb->running == NULL ) {
		ERROR("Couldn't allocate memory for process flags.\n");
		return;
	}

	sb->last_filename = calloc(n_proc, sizeof(char *));
	if ( sb->last_filename == NULL ) {
		ERROR("Couldn't allocate memory for last filename list.\n");
		return;
	}
	if ( sb->fhs == NULL ) {
		ERROR("Couldn't allocate memory for file handles!\n");
		return;
	}
	pthread_mutex_unlock(&sb->lock);

	if ( pthread_create(&reader_thread, NULL, run_reader, (void *)sb) ) {
		ERROR("Failed to create reader thread.\n");
		return;
	}

	/* Set up signal handler to take action if any children die */
	sa.sa_flags = SA_SIGINFO | SA_NOCLDSTOP;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = signal_handler;
	r = sigaction(SIGCHLD, &sa, NULL);
	if ( r == -1 ) {
		ERROR("Failed to set signal handler!\n");
		return;
	}

	/* Fork the right number of times */
	for ( i=0; i<n_proc; i++ ) {
		start_worker_process(sb, i);
	}

	allDone = 0;
	while ( !allDone ) {

		int r, i;
		struct timeval tv;
		fd_set fds;
		double tNow;
		int fdmax;

		tv.tv_sec = 1;
		tv.tv_usec = 0;

		FD_ZERO(&fds);
		fdmax = 0;
		pthread_mutex_lock(&sb->lock);
		for ( i=0; i<n_proc; i++ ) {

			int fd;

			if ( !sb->running[i] ) {
				continue;
			}

			fd = fileno(sb->result_fhs[i]);
			FD_SET(fd, &fds);
			if ( fd > fdmax ) fdmax = fd;

		}
		pthread_mutex_unlock(&sb->lock);

		r = select(fdmax+1, &fds, NULL, NULL, &tv);
		if ( r == -1 ) {
			if ( errno != EINTR ) {
				ERROR("select() failed: %s\n", strerror(errno));
			}
			continue;
		}

		if ( r == 0 ) continue; /* No progress this time.  Try again */

		pthread_mutex_lock(&sb->lock);
		for ( i=0; i<n_proc; i++ ) {

			char *nextImage;
			char results[1024];
			char *rval;
			int fd;
			int n;
			char *eptr;

			if ( !sb->running[i] ) {
				continue;
			}

			fd = fileno(sb->result_fhs[i]);
			if ( !FD_ISSET(fd, &fds) ) {
				continue;
			}

			rval = fgets(results, 1024, sb->result_fhs[i]);
			if ( rval == NULL ) {
				if ( !feof(sb->result_fhs[i]) ) {
					ERROR("fgets() failed: %s\n",
					      strerror(errno));
				}
				continue;
			}

			chomp(results);

			n = strtol(results, &eptr, 10);
			if ( eptr == results ) {
				if ( strlen(results) > 0 ) {
					ERROR("Invalid result '%s'\n", results);
				}
			} else {
				sb->n_indexable += atoi(results);
				sb->n_processed++;
			}

			/* Send next filename */
			nextImage = get_pattern(fh, &use_this_one_instead,
	                                        config_basename, prefix);

			if ( nextImage == NULL ) {
				/* No more images */
				r = write(sb->filename_pipes[i], "\n", 1);
				if ( r < 0 ) {
					ERROR("Write pipe\n");
				}
			} else {
				r = write(sb->filename_pipes[i], nextImage,
				          strlen(nextImage));
				r -= write(sb->filename_pipes[i], "\n", 1);
				if ( r < 0 ) {
					ERROR("write pipe\n");
				}
				free(nextImage);
			}

		}

		/* Update progress */
		tNow = get_monotonic_seconds();
		if ( tNow >= sb->t_last_stats+STATS_EVERY_N_SECONDS ) {

			STATUS("%i out of %i indexed so far,"
			       " %i out of %i since the last message.\n",
			       sb->n_indexable, sb->n_processed,
			       sb->n_indexable - sb->n_indexable_last_stats,
			       sb->n_processed - sb->n_processed_last_stats);

			sb->n_indexable_last_stats = sb->n_indexable;
			sb->n_processed_last_stats = sb->n_processed;
			sb->t_last_stats = tNow;

		}

		allDone = 1;
		for ( i=0; i<n_proc; i++ ) {
			if ( sb->running[i] ) allDone = 0;
		}

		pthread_mutex_unlock(&sb->lock);

	}

	fclose(fh);

	pthread_mutex_destroy(&sb->lock);

	for ( i=0; i<n_proc; i++ ) {
		int status;
		waitpid(sb->pids[i], &status, 0);
	}

	for ( i=0; i<n_proc; i++ ) {
		close(sb->filename_pipes[i]);
		fclose(sb->result_fhs[i]);
	}

	for ( i=0; i<sb->n_proc; i++ ) {
		fclose(sb->fhs[i]);
	}
	free(sb->fhs);
	free(sb->filename_pipes);
	free(sb->result_fhs);
	free(sb->pids);
	free(sb->running);

	if ( ofh != stdout ) fclose(ofh);

	STATUS("There were %i images, of which %i could be indexed.\n",
	       sb->n_processed, sb->n_indexable);

}
