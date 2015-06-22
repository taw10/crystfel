/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2015 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2015 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
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
#include <pthread.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <signal.h>
#include <sys/stat.h>
#include <assert.h>
#include <sys/mman.h>
#include <semaphore.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <events.h>
#include <hdf5-file.h>
#include <detector.h>

#include "im-sandbox.h"
#include "process_image.h"


/* Write statistics at APPROXIMATELY this interval */
#define STATS_EVERY_N_SECONDS (5)

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

	/* Final output */
	Stream *stream;
};


struct sandbox
{
	int n_processed_last_stats;
	int n_hadcrystals_last_stats;
	int n_crystals_last_stats;
	int t_last_stats;

	struct index_args *iargs;

	int n_proc;
	pid_t *pids;

	int *running;
	struct filename_plus_event **last_event;
	int serial;

	struct sb_shm *shared;

	char *tmpdir;

	struct sb_reader *reader;
};


/* Horrible global variable for signal handler */
sem_t zombie_sem;


static struct filename_plus_event *get_pattern(FILE *fh, int config_basename,
                                               struct detector *det,
                                               const char *prefix)
{
	char *line = NULL;
	size_t len;
	struct filename_plus_event *fne;
	struct hdfile *hdfile;
	char filename_buf[2014];
	char event_buf[2014];

	static char *filename = NULL;
	static struct event_list *ev_list = NULL;
	static int event_index = -1;

	line = malloc(1024*sizeof(char));

	while ( event_index == -1 ) {

		int scan_check;

		do {

			/* Get the next filename */
			char *rval;

			rval = fgets(line, 1023, fh);
			if ( rval == NULL ) {
				free(line);
				free(filename);
				filename = NULL;
				return NULL;
			}

			chomp(line);

		} while ( strlen(line) == 0 );

		if ( config_basename ) {
			char *tmp;
			tmp = safe_basename(line);
			free(line);
			line = tmp;
		}

		scan_check = sscanf(line, "%s %s", filename_buf, event_buf );

		len = strlen(prefix)+strlen(filename_buf)+1;

		/* Round the length of the buffer, to keep Valgrind quiet when
		 * it gets given to write() a bit later on */
		len += 4 - (len % 4);

		if ( filename == NULL ) {
			filename = malloc(len);
		} else {
			char *new_filename;
			new_filename = realloc(filename, len*sizeof(char));
			if ( filename == NULL ) {
				return NULL;
			}
			filename = new_filename;
		}

		snprintf(filename, 1023, "%s%s", prefix, filename_buf);

		if ( det->path_dim != 0  || det->dim_dim != 0 ) {

			ev_list = initialize_event_list();

			if ( scan_check == 1 ) {

				hdfile = hdfile_open(filename);
				if ( hdfile == NULL ) {
					ERROR("Failed to open %s\n", filename);
					free(line);
					return NULL;
				}

				if ( ev_list != NULL ) {
					free_event_list(ev_list);
				}

				ev_list = fill_event_list(hdfile, det);
				if ( ev_list == NULL ) {
					ERROR("Failed to get event list.\n");
					return NULL;
				}

				if ( ev_list->num_events == 0 ) {
					event_index = -1;
				} else {
					event_index = 0;
				}

				hdfile_close(hdfile);

			} else {

				struct event *ev_to_add;

				ev_to_add = get_event_from_event_string(event_buf);
				append_event_to_event_list(ev_list, ev_to_add);
				free_event(ev_to_add);
				event_index = 0;

			}

		} else {

			event_index = 0;

		}
	}

	fne = malloc(sizeof(struct filename_plus_event));
	fne->filename = strdup(filename);

	if ( det->path_dim !=0  || det->dim_dim !=0 ) {
		fne->ev = copy_event(ev_list->events[event_index]);
		if ( event_index != ev_list->num_events-1 ) {
			event_index += 1;
		} else {
			event_index = -1;
		}
	} else {
		fne->ev = NULL;
		event_index = -1;
	}

	free(line);
	return fne;
}


static void shuffle_events(struct sb_shm *sb_shared)
{
	int i;

	for ( i=1; i<sb_shared->n_events; i++ ) {
		memcpy(sb_shared->queue[i-1], sb_shared->queue[i], MAX_EV_LEN);
	}
	sb_shared->n_events--;
}


static void run_work(const struct index_args *iargs, Stream *st,
                     int cookie, const char *tmpdir, struct sb_shm *sb_shared)
{
	int allDone = 0;

	while ( !allDone ) {

		struct pattern_args pargs;
		char filename[MAX_EV_LEN];
		char event_str[MAX_EV_LEN];
		int ser;
		struct event *ev;
		int r;

		/* Wait until an event is ready */
		sem_wait(&sb_shared->queue_sem);

		/* Get the event from the queue */
		pthread_mutex_lock(&sb_shared->queue_lock);
		if ( sb_shared->no_more ) {
			pthread_mutex_unlock(&sb_shared->queue_lock);
			allDone = 1;
			continue;
		}
		r = sscanf(sb_shared->queue[0], "%s %s %i",
		           filename, event_str, &ser);
		if ( r != 3 ) {
			STATUS("Invalid event string '%s'\n",
			       sb_shared->queue[0]);
		}
		memcpy(sb_shared->last_ev[cookie], sb_shared->queue[0],
		       MAX_EV_LEN);
		shuffle_events(sb_shared);
		pthread_mutex_unlock(&sb_shared->queue_lock);

		if ( r != 3 ) continue;

		pargs.filename_p_e = initialize_filename_plus_event();
		pargs.filename_p_e->filename = strdup(filename);

		if ( strcmp(event_str, "/") != 0 ) {

			ev = get_event_from_event_string(event_str);
			if ( ev == NULL ) {
				ERROR("Bad event string '%s'\n", event_str);
				continue;
			}
			pargs.filename_p_e->ev = ev;

		} else {

			pargs.filename_p_e->ev = NULL;

		}

		process_image(iargs, &pargs, st, cookie, tmpdir, ser,
		              sb_shared);

		free_filename_plus_event(pargs.filename_p_e);

	}

	cleanup_indexing(iargs->indm, iargs->ipriv);
	free_detector_geometry(iargs->det);
	free(iargs->hdf5_peak_path);
	free_copy_hdf5_field_list(iargs->copyme);
	cell_free(iargs->cell);
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


static ssize_t lwrite(int fd, const char *a)
{
	size_t l = strlen(a);
	return write(fd, a, l);
}


static int pump_chunk(FILE *fh, int ofd)
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
					lwrite(ofd, "Unfinished chunk!\n");
					lwrite(ofd, CHUNK_END_MARKER"\n");
				} /* else normal end of output */
				return 1;
			}

			ERROR("fgets() failed: %s\n", strerror(errno));
			return 1;

		}

		if ( strcmp(line, "FLUSH\n") == 0 ) break;
		lwrite(ofd, line);

		if ( strcmp(line, CHUNK_END_MARKER"\n") == 0 ) break;

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

	fclose(rd->fhs[d]);

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
	const int ofd = get_stream_fd(rd->stream);

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
			if ( pump_chunk(rd->fhs[i], ofd) ) {
				/* remove_pipe() assumes that the caller is
				 * holding rd->lock ! */
				remove_pipe(rd, i);
			}

		}
		pthread_mutex_unlock(&rd->lock);

	}

	return NULL;
}


static void start_worker_process(struct sandbox *sb, int slot)
{
	pid_t p;
	int stream_pipe[2];

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
		struct sigaction sa;
		int r;
		char *tmp;
		struct stat s;
		size_t ll;
		int i;

		/* First, disconnect the signal handler */
		sa.sa_flags = 0;
		sigemptyset(&sa.sa_mask);
		sa.sa_handler = SIG_DFL;
		r = sigaction(SIGCHLD, &sa, NULL);
		if ( r == -1 ) {
			ERROR("Failed to set signal handler!\n");
			return;
		}

		ll = 64 + strlen(sb->tmpdir);
		tmp = malloc(ll);
		if ( tmp == NULL ) {
			ERROR("Failed to allocate temporary dir\n");
			return;
		}

		snprintf(tmp, 63, "%s/worker.%i", sb->tmpdir, slot);

		if ( stat(tmp, &s) == -1 ) {
			if ( errno != ENOENT ) {
				ERROR("Failed to stat temporary folder.\n");
				return;
			}

			r = mkdir(tmp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if ( r ) {
				ERROR("Failed to create temporary folder: %s\n",
				strerror(errno));
				return;
			}
		}

		/* Free resources which will not be needed by worker */
		free(sb->pids);
		for ( i=0; i<sb->reader->n_read; i++ ) {
			fclose(sb->reader->fhs[i]);
		}
		free(sb->reader->fhs);
		free(sb->reader->fds);
		free(sb->reader);
		free(sb->tmpdir);
		free(sb->running);
		/* Not freed because it's not worth passing them down just for
		 * this purpose: event list file handle,
		 *               main output stream handle
		 *               original temp dir name (without indexamajig.XX)
		 *               prefix
		 */

		st = open_stream_fd_for_write(stream_pipe[1]);
		run_work(sb->iargs, st, slot, tmp, sb->shared);
		close_stream(st);

		free(tmp);
		free(sb->iargs->beam->photon_energy_from);

		munmap(sb->shared, sizeof(struct sb_shm));

		free(sb);

		exit(0);

	}

	/* Parent process gets the 'write' end of the filename pipe
	 * and the 'read' end of the result pipe. */
	sb->pids[slot] = p;
	sb->running[slot] = 1;
	add_pipe(sb->reader, stream_pipe[0]);
	close(stream_pipe[1]);
}


static void signal_handler(int sig, siginfo_t *si, void *uc_v)
{
	sem_post(&zombie_sem);
}


static void handle_zombie(struct sandbox *sb)
{
	int i;

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
				STATUS("Event ID was: %s\n",
				       sb->shared->last_ev[i]);
				start_worker_process(sb, i);
			}

		}

	}
}


static int setup_shm(struct sandbox *sb)
{
	pthread_mutexattr_t attr;

	sb->shared = mmap(NULL, sizeof(struct sb_shm), PROT_READ | PROT_WRITE,
	                  MAP_SHARED | MAP_ANON, -1, 0);

	if ( sb->shared == MAP_FAILED ) {
		ERROR("SHM setup failed: %s\n", strerror(errno));
		return 1;
	}

	if ( pthread_mutexattr_init(&attr) ) {
		ERROR("Failed to initialise mutex attr.\n");
		return 1;
	}

	if ( pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED) ) {
		ERROR("Failed to set process shared attribute.\n");
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->term_lock, &attr) ) {
		ERROR("Terminal lock setup failed.\n");
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->queue_lock, &attr) ) {
		ERROR("Queue lock setup failed.\n");
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->totals_lock, &attr) ) {
		ERROR("Totals lock setup failed.\n");
		return 1;
	}

	pthread_mutexattr_destroy(&attr);

	return 0;
}


static char *maybe_get_event_string(struct event *ev)
{
	if ( ev == NULL ) return "/";
	return get_event_string(ev);
}


/* Assumes the caller is already holding queue_lock! */
static int fill_queue(FILE *fh, int config_basename, struct detector *det,
                      const char *prefix, struct sandbox *sb)
{
	while ( sb->shared->n_events < QUEUE_SIZE ) {

		struct filename_plus_event *ne;
		char ev_string[MAX_EV_LEN];

		ne = get_pattern(fh, config_basename, det, prefix);
		if ( ne == NULL ) return 1; /* No more */

		memset(ev_string, 0, MAX_EV_LEN);
		snprintf(ev_string, MAX_EV_LEN, "%s %s %i", ne->filename,
		         maybe_get_event_string(ne->ev), sb->serial++);
		memcpy(sb->shared->queue[sb->shared->n_events++], ev_string,
		       MAX_EV_LEN);
		sem_post(&sb->shared->queue_sem);
		free_filename_plus_event(ne);

	}
	return 0;
}


void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                    int config_basename, FILE *fh,
                    Stream *stream, const char *tempdir)
{
	int i;
	struct sigaction sa;
	int r;
	pthread_t reader_thread;
	struct sandbox *sb;
	size_t ll;
	struct stat s;
	int allDone = 0;

	if ( n_proc > MAX_NUM_WORKERS ) {
		ERROR("Number of workers (%i) is too large.  Using %i\n",
		      n_proc, MAX_NUM_WORKERS);
		n_proc = MAX_NUM_WORKERS;
	}

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

	pthread_mutex_init(&sb->reader->lock, NULL);

	sb->n_processed_last_stats = 0;
	sb->n_hadcrystals_last_stats = 0;
	sb->n_crystals_last_stats = 0;
	sb->t_last_stats = get_monotonic_seconds();
	sb->n_proc = n_proc;
	sb->iargs = iargs;
	sb->serial = 1;

	sb->reader->fds = NULL;
	sb->reader->fhs = NULL;
	sb->reader->stream = stream;

	if ( setup_shm(sb) ) {
		ERROR("Failed to set up SHM.\n");
		free(sb);
		return;
	}

	sb->shared->n_processed = 0;
	sb->shared->n_hadcrystals = 0;
	sb->shared->n_crystals = 0;

	sem_init(&sb->shared->queue_sem, 1, 0);
	sem_init(&zombie_sem, 0, 0);

	sb->pids = calloc(n_proc, sizeof(pid_t));
	sb->running = calloc(n_proc, sizeof(int));
	if ( sb->pids == NULL ) {
		ERROR("Couldn't allocate memory for PIDs.\n");
		return;
	}
	if ( sb->running == NULL ) {
		ERROR("Couldn't allocate memory for process flags.\n");
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

	if ( tempdir == NULL ) {
		tempdir = "";
	}

	ll = 64+strlen(tempdir);
	sb->tmpdir = malloc(ll);
	if ( sb->tmpdir == NULL ) {
		ERROR("Failed to allocate temporary directory name\n");
		return;
	}
	snprintf(sb->tmpdir, ll, "%s/indexamajig.%i", tempdir, getpid());

	if ( stat(sb->tmpdir, &s) == -1 ) {

		int r;

		if ( errno != ENOENT ) {
			ERROR("Failed to stat temporary folder.\n");
			return;
		}

		r = mkdir(sb->tmpdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if ( r ) {
			ERROR("Failed to create temporary folder: %s\n",
			      strerror(errno));
			return;
		}

	}

	/* Fill the queue */
	pthread_mutex_lock(&sb->shared->queue_lock);
	sb->shared->n_events = 0;
	fill_queue(fh, config_basename, iargs->det, prefix, sb);
	sb->shared->no_more = 0;
	pthread_mutex_unlock(&sb->shared->queue_lock);

	/* Fork the right number of times */
	for ( i=0; i<n_proc; i++ ) {
		start_worker_process(sb, i);
	}

	/* Start reader thread after forking, so that things are definitely
	 * "running" */
	if ( pthread_create(&reader_thread, NULL, run_reader,
	                    (void *)sb->reader) ) {
		ERROR("Failed to create reader thread.\n");
		return;
	}

	do {

		int r;
		double tNow;

		sleep(5);

		/* Check for dead workers */
		if ( sem_trywait(&zombie_sem) == 0 ) handle_zombie(sb);

		/* Top up the queue if necessary */
		r = 0;
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( sb->shared->n_events < QUEUE_SIZE/2 ) {
			r = fill_queue(fh, config_basename, iargs->det, prefix,
			               sb);
		}
		pthread_mutex_unlock(&sb->shared->queue_lock);

		/* Update progress */
		tNow = get_monotonic_seconds();
		r = pthread_mutex_trylock(&sb->shared->term_lock);
		if ( r == 0 ) {

			/* Could get lock, so write status */
			int n_proc_this;
			double indexable;

			n_proc_this = sb->shared->n_processed
			              - sb->n_processed_last_stats;
			indexable = (sb->shared->n_processed == 0) ? 0 :
			            100.0 * sb->shared->n_hadcrystals
			              / sb->shared->n_processed;

			STATUS("%4i indexable out of %4i processed (%4.1f%%), "
			       "%4i crystals so far. "
			       "%4i images processed since the last message.\n",
			       sb->shared->n_hadcrystals,
			       sb->shared->n_processed, indexable,
			       sb->shared->n_crystals, n_proc_this);

			sb->n_processed_last_stats = sb->shared->n_processed;
			sb->n_hadcrystals_last_stats = sb->shared->n_hadcrystals;
			sb->n_crystals_last_stats = sb->shared->n_crystals;
			sb->t_last_stats = tNow;

			pthread_mutex_unlock(&sb->shared->term_lock);

		}

		/* Have all the events been swallowed? */
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( sb->shared->n_events == 0 ) allDone = 1;
		pthread_mutex_unlock(&sb->shared->queue_lock);

	} while ( !allDone );
	fclose(fh);

	/* Indicate to the workers that we are finished, and wake them up one
	 * last time */
	STATUS("Waiting for the last patterns to be processed...\n");
	pthread_mutex_lock(&sb->shared->queue_lock);
	sb->shared->no_more = 1;
	pthread_mutex_unlock(&sb->shared->queue_lock);
	for ( i=0; i<n_proc; i++ ) {
		sem_post(&sb->shared->queue_sem);
	}
	for ( i=0; i<n_proc; i++ ) {
		int status;
		waitpid(sb->pids[i], &status, 0);
	}

	/* Indicate to the reader thread that we are done */
	pthread_mutex_lock(&sb->reader->lock);
	sb->reader->done = 1;
	pthread_mutex_unlock(&sb->reader->lock);

	pthread_join(reader_thread, NULL);

	for ( i=0; i<sb->reader->n_read; i++ ) {
		fclose(sb->reader->fhs[i]);
	}
	free(sb->reader->fhs);
	free(sb->reader->fds);
	free(sb->running);
	free(sb->pids);
	free(sb->tmpdir);
	free(sb->reader);

	STATUS("Final: %i images processed, %i had crystals (%.1f%%),"
	       " %i crystals overall.\n",
	       sb->shared->n_processed, sb->shared->n_hadcrystals,
	       100.0 * sb->shared->n_hadcrystals / sb->shared->n_processed,
	       sb->shared->n_crystals);

	munmap(sb->shared, sizeof(struct sb_shm));
	free(sb);
}
