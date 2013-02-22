/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2013 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2013 Thomas White <taw@physics.org>
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
#include <sys/stat.h>

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


/* Information about the indexing process for one pattern */
struct pattern_args
{
	/* "Input" */
	char *filename;

	/* "Output" */
	int n_crystals;
};


struct sb_reader
{
	pthread_mutex_t lock;
	int done;

	/* If a worker process dies unexpectedly (e.g. if it segfaults), then
	 * the pipe for its output can still stay open for a little while while
	 * its buffer empties.  The number of pipes being read from is therefore
	 * not necessarily the same as the number of worker processes. */
	int n_read;
	FILE **fhs;
	int *fds;

	/* Final output file handle */
	FILE *ofh;
};


struct sandbox
{
	pthread_mutex_t lock;

	int n_processed;
	int n_hadcrystals;
	int n_crystals;
	int n_processed_last_stats;
	int n_hadcrystals_last_stats;
	int n_crystals_last_stats;
	int t_last_stats;

	struct index_args *iargs;

	int n_proc;
	pid_t *pids;

	int *running;
	FILE **result_fhs;
	int *filename_pipes;
	int *stream_pipe_write;
	char **last_filename;

	struct sb_reader *reader;
};


/* Horrible global variable for signal handler */
int signal_pipe[2];


static void lock_sandbox(struct sandbox *sb)
{
	pthread_mutex_lock(&sb->lock);
}


static void unlock_sandbox(struct sandbox *sb)
{
	pthread_mutex_unlock(&sb->lock);
}


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
                          struct pattern_args *pargs, Stream *st,
                          int cookie)
{
	float *data_for_measurement;
	size_t data_size;
	int check;
	struct hdfile *hdfile;
	struct image image;
	int i;
	char filename[1024];

	/* Prefix to jump out of temporary folder */
	if ( pargs->filename[0] != '/' ) {
		snprintf(filename, 1023, "../../%s", pargs->filename);
	} else {
		snprintf(filename, 1023, "%s", pargs->filename);
	}

	image.features = NULL;
	image.data = NULL;
	image.flags = NULL;
	image.copyme = iargs->copyme;
	image.id = cookie;
	image.filename = pargs->filename;  /* Relative to top level */
	image.beam = iargs->beam;
	image.det = iargs->det;
	image.crystals = NULL;
	image.n_crystals = 0;

	hdfile = hdfile_open(filename);  /* Relative to temporary folder */
	if ( hdfile == NULL ) return;

	if ( iargs->element != NULL ) {

		int r;
		r = hdfile_set_image(hdfile, iargs->element);
		if ( r ) {
			ERROR("Couldn't select path '%s'\n", iargs->element);
			hdfile_close(hdfile);
			return;
		}

	} else {

		int r;
		r = hdfile_set_first_image(hdfile, "/");
		if ( r ) {
			ERROR("Couldn't select first path\n");
			hdfile_close(hdfile);
			return;
		}

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
		return;
	}

	fill_in_values(image.det, hdfile);
	fill_in_beam_parameters(image.beam, hdfile);

	image.lambda = ph_en_to_lambda(eV_to_J(image.beam->photon_energy));

	if ( (image.beam->photon_energy < 0.0) || (image.lambda > 1000) ) {
		/* Error message covers a silly value in the beam file or in
		 * the HDF5 file. */
		ERROR("Nonsensical wavelength (%e m or %e eV) value for %s.\n",
		      image.lambda, image.beam->photon_energy, image.filename);
		hdfile_close(hdfile);
		return;
	}

	if ( iargs->cmfilter ) filter_cm(&image);

	/* Take snapshot of image after CM subtraction but before
	 * the aggressive noise filter. */
	data_size = image.width * image.height * sizeof(float);
	data_for_measurement = malloc(data_size);

	if ( iargs->noisefilter ) {
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
		if ( !iargs->no_revalidate ) {
			validate_peaks(&image, iargs->min_int_snr,
				       iargs->ir_inn, iargs->ir_mid,
				       iargs->ir_out, iargs->use_saturated);
		}
		break;

		case PEAK_ZAEF:
		search_peaks(&image, iargs->threshold,
		             iargs->min_gradient, iargs->min_snr,
		             iargs->ir_inn, iargs->ir_mid, iargs->ir_out,
		             iargs->use_saturated);
		break;

	}

	/* Get rid of noise-filtered version at this point
	 * - it was strictly for the purposes of peak detection. */
	free(image.data);
	image.data = data_for_measurement;

	/* Index the pattern */
	index_pattern(&image, iargs->indm, iargs->ipriv);

	pargs->n_crystals = image.n_crystals;

	/* Default beam parameters */
	image.div = image.beam->divergence;
	image.bw = image.beam->bandwidth;

	/* Integrate each crystal's diffraction spots */
	for ( i=0; i<image.n_crystals; i++ ) {

		RefList *reflections;

		/* Set default crystal parameter(s) */
		crystal_set_profile_radius(image.crystals[i],
		                           image.beam->profile_radius);

		if ( iargs->integrate_found ) {
			reflections = select_intersections(&image,
			                                   image.crystals[i]);
		} else {
			reflections = find_intersections(&image,
			                                 image.crystals[i]);
		}

		crystal_set_reflections(image.crystals[i], reflections);

	}

	/* Integrate all the crystals at once - need all the crystals so that
	 * overlaps can be detected. */
	integrate_reflections(&image, iargs->closer,
	                              iargs->bgsub,
	                              iargs->min_int_snr,
	                              iargs->ir_inn,
	                              iargs->ir_mid,
	                              iargs->ir_out,
	                              iargs->integrate_saturated,
	                              iargs->res_cutoff);

	write_chunk(st, &image, hdfile,
	            iargs->stream_peaks, iargs->stream_refls);

	for ( i=0; i<image.n_crystals; i++ ) {
		crystal_free(image.crystals[i]);
	}

	free(image.data);
	if ( image.flags != NULL ) free(image.flags);
	image_feature_list_free(image.features);
	hdfile_close(hdfile);
}


static void run_work(const struct index_args *iargs,
                     int filename_pipe, int results_pipe, Stream *st,
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

			ERROR("Read error!\n");
			free(line);
			allDone = 1;
			continue;

		}

		chomp(line);

		if ( strlen(line) == 0 ) {

			allDone = 1;

		} else {

			pargs.filename = line;
			pargs.n_crystals = 0;

			process_image(iargs, &pargs, st, cookie);

			/* Request another image */
			c = sprintf(buf, "%i\n", pargs.n_crystals);
			w = write(results_pipe, buf, c);
			if ( w < 0 ) {
				ERROR("write P0\n");
			}

		}

		free(line);

	}

	cleanup_indexing(iargs->indm, iargs->ipriv);
	free(iargs->indm);
	free(iargs->ipriv);
	free_detector_geometry(iargs->det);
	free_beam_parameters(iargs->beam);
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

	do {

		char line[1024];
		char *rval;

		rval = fgets(line, 1024, fh);
		if ( rval == NULL ) {

			if ( feof(fh) ) {
				/* Whoops, connection lost */
				if ( chunk_started ) {
					ERROR("EOF during chunk!\n");
					fprintf(ofh, "Chunk is unfinished!\n");
					fprintf(ofh, CHUNK_END_MARKER"\n");
				} /* else normal end of output */
				return 1;
			} else {
				ERROR("fgets() failed: %s\n", strerror(errno));
			}
			break;

		}

		if ( strcmp(line, "FLUSH\n") == 0 ) break;

		fprintf(ofh, "%s", line);

		if ( strcmp(line, CHUNK_END_MARKER"\n") == 0 ) break;
		if ( strcmp(line, CHUNK_START_MARKER"\n") == 0 ) break;


	} while ( 1 );
	return 0;
}


/* Add an fd to the list of pipes to be read from */
static void add_pipe(struct sb_reader *rd, int fd)
{
	int *fds_new;
	FILE **fhs_new;
	int slot;

	pthread_mutex_lock(&rd->lock);

	fds_new = realloc(rd->fds, (rd->n_read+1)*sizeof(int));
	if ( fds_new == NULL ) {
		ERROR("Failed to allocate memory for new pipe.\n");
		return;
	}

	fhs_new = realloc(rd->fhs, (rd->n_read+1)*sizeof(FILE *));
	if ( fhs_new == NULL ) {
		ERROR("Failed to allocate memory for new FH.\n");
		return;
	}

	rd->fds = fds_new;
	rd->fhs = fhs_new;
	slot = rd->n_read;

	rd->fds[slot] = fd;

	rd->fhs[slot] = fdopen(fd, "r");
	if ( rd->fhs[slot] == NULL ) {
		ERROR("Couldn't fdopen() stream!\n");
		return;
	}

	rd->n_read++;

	pthread_mutex_unlock(&rd->lock);
}


/* Assumes that the caller is already holding rd->lock! */
static void remove_pipe(struct sb_reader *rd, int d)
{
	int i;

	for ( i=d; i<rd->n_read; i++ ) {
		if ( i < rd->n_read-1 ) {
			rd->fds[i] = rd->fds[i+1];
			rd->fhs[i] = rd->fhs[i+1];
		} /* else don't bother */
	}

	rd->n_read--;

	/* We don't bother shrinking the arrays */
}


static void *run_reader(void *rdv)
{
	struct sb_reader *rd = rdv;

	while ( 1 ) {

		int r, i;
		struct timeval tv;
		fd_set fds;
		int fdmax;

		/* Exit when:
		 *     - No fhs left open to read from
		 * AND - Main thread says "done" */
		if ( (rd->n_read == 0) && rd->done ) break;

		tv.tv_sec = 1;
		tv.tv_usec = 0;

		FD_ZERO(&fds);
		fdmax = 0;
		pthread_mutex_lock(&rd->lock);
		for ( i=0; i<rd->n_read; i++ ) {

			int fd;

			fd = rd->fds[i];

			FD_SET(fd, &fds);
			if ( fd > fdmax ) fdmax = fd;

		}
		pthread_mutex_unlock(&rd->lock);

		r = select(fdmax+1, &fds, NULL, NULL, &tv);

		if ( r == -1 ) {
			if ( errno != EINTR ) {
				ERROR("select() failed: %s\n", strerror(errno));
			} /* Otherwise no big deal */
			continue;
		}

		pthread_mutex_lock(&rd->lock);
		for ( i=0; i<rd->n_read; i++ ) {

			if ( !FD_ISSET(rd->fds[i], &fds) ) {
				continue;
			}

			/* If the chunk cannot be read, assume the connection
			 * is broken and that the process will die soon. */
			if ( pump_chunk(rd->fhs[i], rd->ofh) ) {
				/* remove_pipe() assumes that the caller is
				 * holding rd->lock ! */
				remove_pipe(rd, i);
			}

		}
		pthread_mutex_unlock(&rd->lock);

	}

	return NULL;
}


static int create_temporary_folder(signed int n)
{
	int r;
	char tmp[64];

	if ( n < 0 ) {
		snprintf(tmp, 63, "indexamajig.%i", getpid());
	} else {
		snprintf(tmp, 63, "worker.%i", n);
	}

	r = mkdir(tmp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if ( r ) {
		ERROR("Failed to create temporary folder: %s\n",
		      strerror(errno));
		return 1;
	}

	r = chdir(tmp);
	if ( r ) {
		ERROR("Failed to chdir to temporary folder: %s\n",
		      strerror(errno));
		return 1;
	}

	return 0;
}


static void start_worker_process(struct sandbox *sb, int slot,
                                 int argc, char *argv[])
{
	pid_t p;
	int filename_pipe[2];
	int result_pipe[2];
	int stream_pipe[2];

	if ( pipe(filename_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	if ( pipe(result_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	if ( pipe(stream_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	p = fork();
	if ( p == -1 ) {
		ERROR("fork() failed!\n");
		return;
	}

	if ( p == 0 ) {

		Stream *st;
		int j;
		struct sigaction sa;
		int r;

		/* First, disconnect the signal handler */
		sa.sa_flags = 0;
		sigemptyset(&sa.sa_mask);
		sa.sa_handler = SIG_DFL;
		r = sigaction(SIGCHLD, &sa, NULL);
		if ( r == -1 ) {
			ERROR("Failed to set signal handler!\n");
			return;
		}

		create_temporary_folder(slot);

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

		st = open_stream_fd_for_write(stream_pipe[1]);
		write_command(st, argc, argv);
		write_line(st, "FLUSH");
		run_work(sb->iargs, filename_pipe[0], result_pipe[1],
		         st, slot);
		close_stream(st);

		//close(filename_pipe[0]);
		close(result_pipe[1]);

		exit(0);

	}

	/* Parent process gets the 'write' end of the filename pipe
	 * and the 'read' end of the result pipe. */
	sb->pids[slot] = p;
	sb->running[slot] = 1;
	add_pipe(sb->reader, stream_pipe[0]);
	close(filename_pipe[0]);
	close(result_pipe[1]);
	close(stream_pipe[1]);
	sb->filename_pipes[slot] = filename_pipe[1];

	sb->result_fhs[slot] = fdopen(result_pipe[0], "r");
	if ( sb->result_fhs[slot] == NULL ) {
		ERROR("fdopen() failed.\n");
		return;
	}
}


static void signal_handler(int sig, siginfo_t *si, void *uc_v)
{
	write(signal_pipe[1], "\n", 1);
}


static void handle_zombie(struct sandbox *sb)
{
	int i;

	lock_sandbox(sb);
	for ( i=0; i<sb->n_proc; i++ ) {

		int status, p;

		if ( !sb->running[i] ) continue;

		p = waitpid(sb->pids[i], &status, WNOHANG);

		if ( p == -1 ) {
			ERROR("waitpid() failed.\n");
			continue;
		}

		if ( p == sb->pids[i] ) {

			sb->running[i] = 0;

			if ( WIFEXITED(status) ) {
				continue;
			}

			if ( WIFSIGNALED(status) ) {
				STATUS("Worker %i was killed by signal %i\n",
				       i, WTERMSIG(status));
				STATUS("Last filename was: %s\n",
				       sb->last_filename[i]);
				sb->n_processed++;
				start_worker_process(sb, i, 0, NULL);
			}

		}

	}
	unlock_sandbox(sb);
}


void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                    int config_basename, FILE *fh, char *use_this_one_instead,
                    FILE *ofh, int argc, char *argv[])
{
	int i;
	int allDone;
	struct sigaction sa;
	int r;
	pthread_t reader_thread;
	struct sandbox *sb;

	sb = calloc(1, sizeof(struct sandbox));
	if ( sb == NULL ) {
		ERROR("Couldn't allocate memory for sandbox.\n");
		return;
	}

	sb->reader = calloc(1, sizeof(struct sb_reader));
	if ( sb->reader == NULL ) {
		ERROR("Couldn't allocate memory for SB reader.\n");
		free(sb);
		return;
	}

	pthread_mutex_init(&sb->lock, NULL);
	pthread_mutex_init(&sb->reader->lock, NULL);

	sb->n_processed = 0;
	sb->n_hadcrystals = 0;
	sb->n_crystals = 0;
	sb->n_processed_last_stats = 0;
	sb->n_hadcrystals_last_stats = 0;
	sb->n_crystals_last_stats = 0;
	sb->t_last_stats = get_monotonic_seconds();
	sb->n_proc = n_proc;
	sb->iargs = iargs;

	sb->reader->fds = NULL;
	sb->reader->fhs = NULL;
	sb->reader->ofh = ofh;

	sb->stream_pipe_write = calloc(n_proc, sizeof(int));
	if ( sb->stream_pipe_write == NULL ) {
		ERROR("Couldn't allocate memory for pipes.\n");
		return;
	}

	lock_sandbox(sb);
	sb->filename_pipes = calloc(n_proc, sizeof(int));
	sb->result_fhs = calloc(n_proc, sizeof(FILE *));
	sb->pids = calloc(n_proc, sizeof(pid_t));
	sb->running = calloc(n_proc, sizeof(int));
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
	unlock_sandbox(sb);

	if ( pipe(signal_pipe) == -1 ) {
		ERROR("Failed to create signal pipe.\n");
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

	if ( create_temporary_folder(-1) ) return;

	/* Fork the right number of times */
	lock_sandbox(sb);
	for ( i=0; i<n_proc; i++ ) {
		start_worker_process(sb, i, argc, argv);
	}
	unlock_sandbox(sb);

	/* Start reader thread after forking, so that things are definitely
	 * "running" */
	if ( pthread_create(&reader_thread, NULL, run_reader,
	                    (void *)sb->reader) ) {
		ERROR("Failed to create reader thread.\n");
		return;
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
		lock_sandbox(sb);
		for ( i=0; i<n_proc; i++ ) {

			int fd;

			if ( sb->result_fhs[i] == NULL) continue;

			fd = fileno(sb->result_fhs[i]);
			FD_SET(fd, &fds);
			if ( fd > fdmax ) fdmax = fd;

		}
		unlock_sandbox(sb);

		FD_SET(signal_pipe[0], &fds);
		if ( signal_pipe[0] > fdmax ) fdmax = signal_pipe[0];

		r = select(fdmax+1, &fds, NULL, NULL, &tv);
		if ( r == -1 ) {
			if ( errno == EINTR ) continue;
			ERROR("select() failed: %s\n", strerror(errno));
			break;
		}

		if ( r == 0 ) continue; /* No progress this time.  Try again */

		if ( FD_ISSET(signal_pipe[0], &fds) ) {

			char d;
			read(signal_pipe[0], &d, 1);
			handle_zombie(sb);

		}

		lock_sandbox(sb);
		for ( i=0; i<n_proc; i++ ) {

			char *nextImage;
			char results[1024];
			char *rval;
			int fd;
			char *eptr;

			if ( sb->result_fhs[i] == NULL ) continue;

			fd = fileno(sb->result_fhs[i]);
			if ( !FD_ISSET(fd, &fds) ) continue;

			rval = fgets(results, 1024, sb->result_fhs[i]);
			if ( rval == NULL ) {
				if ( !feof(sb->result_fhs[i]) ) {
					ERROR("fgets() failed: %s\n",
					      strerror(errno));
				}
				sb->result_fhs[i] = NULL;
				continue;
			}

			chomp(results);

			strtol(results, &eptr, 10);
			if ( eptr == results ) {
				if ( strlen(results) > 0 ) {
					ERROR("Invalid result '%s'\n", results);
				}
			} else {

				int nc = atoi(results);
				sb->n_crystals += nc;
				if ( nc > 0 ) {
					sb->n_hadcrystals++;
				}
				sb->n_processed++;

			}

			/* Send next filename */
			nextImage = get_pattern(fh, &use_this_one_instead,
	                                        config_basename, prefix);

			free(sb->last_filename[i]);
			sb->last_filename[i] = nextImage;

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
			}

		}
		unlock_sandbox(sb);

		/* Update progress */
		lock_sandbox(sb);
		tNow = get_monotonic_seconds();
		if ( tNow >= sb->t_last_stats+STATS_EVERY_N_SECONDS ) {

			STATUS("%4i indexable out of %4i processed, "
			       "%4i crystals so far. "
			       "%4i images processed since the last message.\n",
			       sb->n_hadcrystals, sb->n_processed,
			       sb->n_crystals,
			       sb->n_processed - sb->n_processed_last_stats);

			sb->n_processed_last_stats = sb->n_processed;
			sb->n_hadcrystals_last_stats = sb->n_hadcrystals;
			sb->n_crystals_last_stats = sb->n_crystals;
			sb->t_last_stats = tNow;

		}
		unlock_sandbox(sb);

		allDone = 1;
		lock_sandbox(sb);
		for ( i=0; i<n_proc; i++ ) {
			if ( sb->running[i] ) allDone = 0;
		}

		unlock_sandbox(sb);

	}

	fclose(fh);

	/* Indicate to the reader thread that we are done */
	pthread_mutex_lock(&sb->reader->lock);
	sb->reader->done = 1;
	pthread_mutex_unlock(&sb->reader->lock);

	pthread_join(reader_thread, NULL);

	for ( i=0; i<n_proc; i++ ) {
		int status;
		waitpid(sb->pids[i], &status, 0);
	}

	for ( i=0; i<n_proc; i++ ) {
		close(sb->filename_pipes[i]);
		if ( sb->result_fhs[i] != NULL ) fclose(sb->result_fhs[i]);
	}

	free(sb->running);
	free(sb->filename_pipes);
	free(sb->result_fhs);
	free(sb->pids);

	pthread_mutex_destroy(&sb->lock);

	STATUS("Final:"
	       " %i images processed, %i had crystals, %i crystals overall.\n",
	       sb->n_processed, sb->n_hadcrystals, sb->n_crystals);

	free(sb);
}
