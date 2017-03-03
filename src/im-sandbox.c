/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2017 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2017 Thomas White <taw@physics.org>
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
#include "time-accounts.h"


struct sandbox
{
	int n_processed_last_stats;
	int t_last_stats;

	struct index_args *iargs;

	/* Worker processes */
	int n_proc;
	pid_t *pids;
	int *running;
	time_t *last_response;
	int last_ping[MAX_NUM_WORKERS];

	/* Streams to read from (NB not the same indices as the above) */
	int n_read;
	FILE **fhs;
	int *fds;

	int serial;

	struct sb_shm *shared;
	sem_t *queue_sem;

	char *tmpdir;

	/* Final output */
	Stream *stream;
};


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


static void stamp_response(struct sandbox *sb, int n)
{
	sb->last_response[n] = get_monotonic_seconds();
	sb->last_ping[n] = sb->shared->pings[n];
}


static void check_hung_workers(struct sandbox *sb)
{
	int i;
	time_t tnow = get_monotonic_seconds();
	for ( i=0; i<sb->n_proc; i++ ) {

		if ( !sb->running[i] ) continue;

		if ( sb->shared->pings[i] != sb->last_ping[i] ) {
			stamp_response(sb, i);
		}

		if ( tnow - sb->last_response[i] > 240 ) {
			STATUS("Worker %i did not respond for 240 seconds - "
			       "sending it SIGKILL.\n", i);
			kill(sb->pids[i], SIGKILL);
			stamp_response(sb, i);
		}

	}
}


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
                     int cookie, const char *tmpdir, struct sandbox *sb)
{
	int allDone = 0;
	TimeAccounts *taccs = time_accounts_init();

	while ( !allDone ) {

		struct pattern_args pargs;
		char filename[MAX_EV_LEN];
		char event_str[MAX_EV_LEN];
		int ser;
		struct event *ev;
		int r;

		/* Wait until an event is ready */
		time_accounts_set(taccs, TACC_EVENTWAIT);
		if ( sem_wait(sb->queue_sem) != 0 ) {
			ERROR("Failed to wait on queue semaphore: %s\n",
			      strerror(errno));
		}

		/* Get the event from the queue */
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( sb->shared->no_more ) {
			pthread_mutex_unlock(&sb->shared->queue_lock);
			allDone = 1;
			continue;
		}
		if ( sb->shared->n_events == 0 ) {
			ERROR("Got the semaphore, but no events in queue!\n");
			ERROR("no_more = %i\n", sb->shared->no_more);
			pthread_mutex_unlock(&sb->shared->queue_lock);
			allDone = 1;
			continue;
		}
		r = sscanf(sb->shared->queue[0], "%s %s %i",
		           filename, event_str, &ser);
		if ( r != 3 ) {
			STATUS("Invalid event string '%s'\n",
			       sb->shared->queue[0]);
		}
		memcpy(sb->shared->last_ev[cookie], sb->shared->queue[0],
		       MAX_EV_LEN);
		shuffle_events(sb->shared);
		pthread_mutex_unlock(&sb->shared->queue_lock);

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
		              sb->shared, taccs);

		free_filename_plus_event(pargs.filename_p_e);

	}

	time_accounts_set(taccs, TACC_FINALCLEANUP);
	cleanup_indexing(iargs->ipriv);
	free_detector_geometry(iargs->det);
	free(iargs->hdf5_peak_path);
	free_copy_hdf5_field_list(iargs->copyme);
	cell_free(iargs->cell);
	if ( iargs->profile ) time_accounts_print(taccs);
	time_accounts_free(taccs);
}


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
			if ( errno != EINTR ) return 1;

		}

		if ( strcmp(line, "FLUSH\n") == 0 ) break;
		lwrite(ofd, line);

		if ( strcmp(line, CHUNK_END_MARKER"\n") == 0 ) break;

	} while ( 1 );
	return 0;
}


/* Add an fd to the list of pipes to be read from */
static void add_pipe(struct sandbox *sb, int fd)
{
	int *fds_new;
	FILE **fhs_new;
	int slot;

	fds_new = realloc(sb->fds, (sb->n_read+1)*sizeof(int));
	if ( fds_new == NULL ) {
		ERROR("Failed to allocate memory for new pipe.\n");
		return;
	}

	fhs_new = realloc(sb->fhs, (sb->n_read+1)*sizeof(FILE *));
	if ( fhs_new == NULL ) {
		ERROR("Failed to allocate memory for new FH.\n");
		return;
	}

	sb->fds = fds_new;
	sb->fhs = fhs_new;
	slot = sb->n_read;

	sb->fds[slot] = fd;

	sb->fhs[slot] = fdopen(fd, "r");
	if ( sb->fhs[slot] == NULL ) {
		ERROR("Couldn't fdopen() stream!\n");
		return;
	}

	sb->n_read++;
}


static void remove_pipe(struct sandbox *sb, int d)
{
	int i;

	fclose(sb->fhs[d]);

	for ( i=d; i<sb->n_read; i++ ) {
		if ( i < sb->n_read-1 ) {
			sb->fds[i] = sb->fds[i+1];
			sb->fhs[i] = sb->fhs[i+1];
		} /* else don't bother */
	}

	sb->n_read--;

	/* We don't bother shrinking the arrays */
}


static void try_read(struct sandbox *sb, TimeAccounts *taccs)
{
	int r, i;
	struct timeval tv;
	fd_set fds;
	int fdmax;
	const int ofd = get_stream_fd(sb->stream);

	time_accounts_set(taccs, TACC_SELECT);

	tv.tv_sec = 0;
	tv.tv_usec = 500000;

	FD_ZERO(&fds);
	fdmax = 0;
	for ( i=0; i<sb->n_read; i++ ) {

		int fd;

		fd = sb->fds[i];

		FD_SET(fd, &fds);
		if ( fd > fdmax ) fdmax = fd;

	}

	r = select(fdmax+1, &fds, NULL, NULL, &tv);

	if ( r == -1 ) {
		if ( errno != EINTR ) {
			ERROR("select() failed: %s\n", strerror(errno));
		} /* Otherwise no big deal */
		return;
	}

	for ( i=0; i<sb->n_read; i++ ) {

		if ( !FD_ISSET(sb->fds[i], &fds) ) {
			continue;
		}

		/* If the chunk cannot be read, assume the connection
		 * is broken and that the process will die soon. */
		time_accounts_set(taccs, TACC_STREAMREAD);
		if ( pump_chunk(sb->fhs[i], ofd) ) {
			remove_pipe(sb, i);
		}

	}
}


static void start_worker_process(struct sandbox *sb, int slot)
{
	pid_t p;
	int stream_pipe[2];

	if ( pipe(stream_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	sb->shared->pings[slot] = 0;
	sb->last_ping[slot] = 0;

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

	        /* First, disconnect the signal handlers */
	        sa.sa_flags = 0;
	        sigemptyset(&sa.sa_mask);
	        sa.sa_handler = SIG_DFL;
	        r = sigaction(SIGCHLD, &sa, NULL);
	        if ( r == -1 ) {
			ERROR("Failed to set signal handler!\n");
			exit(1);
	        }
	        r = sigaction(SIGINT, &sa, NULL);
	        if ( r == -1 ) {
			ERROR("Failed to set signal handler!\n");
			exit(1);
	        }
	        r = sigaction(SIGQUIT, &sa, NULL);
	        if ( r == -1 ) {
			ERROR("Failed to set signal handler!\n");
			exit(1);
	        }

		ll = 64 + strlen(sb->tmpdir);
		tmp = malloc(ll);
		if ( tmp == NULL ) {
			ERROR("Failed to allocate temporary dir\n");
			exit(1);
		}

		snprintf(tmp, 63, "%s/worker.%i", sb->tmpdir, slot);

		if ( stat(tmp, &s) == -1 ) {
			if ( errno != ENOENT ) {
				ERROR("Failed to stat temporary folder.\n");
				exit(1);
			}

			r = mkdir(tmp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if ( r ) {
				ERROR("Failed to create temporary folder: %s\n",
				strerror(errno));
				exit(1);
			}
		}

		/* Free resources which will not be needed by worker */
		free(sb->pids);
		for ( i=0; i<sb->n_read; i++ ) {
			fclose(sb->fhs[i]);
		}
		free(sb->fhs);
		free(sb->fds);
		free(sb->tmpdir);
		free(sb->running);
		/* Not freed because it's not worth passing them down just for
		 * this purpose: event list file handle,
		 *               main output stream handle
		 *               original temp dir name (without indexamajig.XX)
		 *               prefix
		 */

		st = open_stream_fd_for_write(stream_pipe[1]);
		run_work(sb->iargs, st, slot, tmp, sb);
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
	stamp_response(sb, slot);
	add_pipe(sb, stream_pipe[0]);
	close(stream_pipe[1]);
}


static void handle_zombie(struct sandbox *sb, int respawn)
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

				if ( (WTERMSIG(status) == SIGINT)
				  || (WTERMSIG(status) == SIGQUIT) ) continue;

				STATUS("Worker %i was killed by signal %i\n",
				       i, WTERMSIG(status));
				STATUS("Event ID was: %s\n",
				       sb->shared->last_ev[i]);
				if ( respawn ) start_worker_process(sb, i);
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
		sem_post(sb->queue_sem);
		free_filename_plus_event(ne);

	}
	return 0;
}

volatile sig_atomic_t at_zombies = 0;
volatile sig_atomic_t at_interrupt = 0;

static void sigchld_handler(int sig, siginfo_t *si, void *uc_v)
{
	at_zombies = 1;
}


static void sigint_handler(int sig, siginfo_t *si, void *uc_v)
{
	at_interrupt = 1;
}


static void check_signals(struct sandbox *sb, const char *semname_q,
                          int respawn)
{
	if ( at_zombies ) {
		at_zombies = 0;
		handle_zombie(sb, respawn);
	}

	if ( at_interrupt ) {
		sem_unlink(semname_q);
		exit(0);
	}
}


static void try_status(struct sandbox *sb, time_t tNow)
{
	int r;
	int n_proc_this;
	double indexable;

	r = pthread_mutex_trylock(&sb->shared->term_lock);
	if ( r ) return; /* No lock -> don't bother */

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
	sb->t_last_stats = tNow;

	pthread_mutex_unlock(&sb->shared->term_lock);
}


void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                    int config_basename, FILE *fh,
                    Stream *stream, const char *tempdir)
{
	int i;
	struct sandbox *sb;
	size_t ll;
	struct stat s;
	char semname_q[64];
	struct sigaction sa;
	int r;
	int no_more = 0;
	int allDone = 0;
	TimeAccounts *taccs;

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

	sb->n_processed_last_stats = 0;
	sb->t_last_stats = get_monotonic_seconds();
	sb->n_proc = n_proc;
	sb->iargs = iargs;
	sb->serial = 1;

	sb->fds = NULL;
	sb->fhs = NULL;
	sb->stream = stream;

	if ( setup_shm(sb) ) {
		ERROR("Failed to set up SHM.\n");
		free(sb);
		return;
	}

	sb->shared->n_processed = 0;
	sb->shared->n_hadcrystals = 0;
	sb->shared->n_crystals = 0;

	/* Set up semaphore to control work queue */
	snprintf(semname_q, 64, "indexamajig-q%i", getpid());
	sb->queue_sem = sem_open(semname_q, O_CREAT | O_EXCL,
	                         S_IRUSR | S_IWUSR, 0);
	if ( sb->queue_sem == SEM_FAILED ) {
		ERROR("Failed to create semaphore: %s\n", strerror(errno));
		return;
	}

	sb->pids = calloc(n_proc, sizeof(pid_t));
	sb->running = calloc(n_proc, sizeof(int));
	sb->last_response = calloc(n_proc, sizeof(time_t));
	if ( (sb->pids == NULL) || (sb->running == NULL)
	  || (sb->last_response == NULL) )
	{
		ERROR("Couldn't allocate memory for PIDs.\n");
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

	/* Set up signal handler to take action if any children die */
	sa.sa_flags = SA_SIGINFO | SA_NOCLDSTOP | SA_RESTART;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = sigchld_handler;
	r = sigaction(SIGCHLD, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return;
	}

	/* Set up signal handler to clean up semaphore on exit */
	sa.sa_flags = SA_SIGINFO | SA_NOCLDSTOP | SA_RESTART;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = sigint_handler;
	r = sigaction(SIGINT, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return;
	}
	r = sigaction(SIGQUIT, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return;
	}

	taccs = time_accounts_init();

	do {

		time_t tNow;

		/* Check for stream output from workers */
		try_read(sb, taccs);

		/* Check for interrupt or zombies */
		time_accounts_set(taccs, TACC_SIGNALS);
		check_signals(sb, semname_q, 1);

		/* Check for hung workers */
		check_hung_workers(sb);

		/* Top up the queue if necessary */
		time_accounts_set(taccs, TACC_QUEUETOPUP);
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( !no_more && (sb->shared->n_events < QUEUE_SIZE/2) ) {
			if ( fill_queue(fh, config_basename, iargs->det,
			                prefix, sb) ) no_more = 1;
		}
		pthread_mutex_unlock(&sb->shared->queue_lock);

		/* Update progress */
		time_accounts_set(taccs, TACC_STATUS);
		tNow = get_monotonic_seconds();
		if ( tNow > sb->t_last_stats+5 ) try_status(sb, tNow);

		/* Have all the events been swallowed? */
		time_accounts_set(taccs, TACC_ENDCHECK);
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( no_more && (sb->shared->n_events == 0) ) allDone = 1;
		pthread_mutex_unlock(&sb->shared->queue_lock);

	} while ( !allDone );

	if ( iargs->profile ) time_accounts_print(taccs);

	fclose(fh);

	/* Indicate to the workers that we are finished, and wake them up one
	 * last time */
	time_accounts_set(taccs, TACC_WAKEUP);
	STATUS("Waiting for the last patterns to be processed...\n");
	pthread_mutex_lock(&sb->shared->queue_lock);
	sb->shared->no_more = 1;
	pthread_mutex_unlock(&sb->shared->queue_lock);
	for ( i=0; i<n_proc; i++ ) {
		sem_post(sb->queue_sem);
	}
	for ( i=0; i<n_proc; i++ ) {
		int status;
		time_accounts_set(taccs, TACC_WAITPID);
		while ( waitpid(sb->pids[i], &status, WNOHANG) == 0 ) {

			time_accounts_set(taccs, TACC_STREAMREAD);
			try_read(sb, taccs);

			time_accounts_set(taccs, TACC_SIGNALS);
			check_signals(sb, semname_q, 0);

			check_hung_workers(sb);

			time_accounts_set(taccs, TACC_WAITPID);
		}
		/* If this worker died and got waited by the zombie handler,
		 * waitpid() returns -1 and the loop still exits. */
	}

	if ( iargs->profile ) time_accounts_print(taccs);
	time_accounts_free(taccs);

	sem_unlink(semname_q);

	for ( i=0; i<sb->n_read; i++ ) {
		fclose(sb->fhs[i]);
	}
	free(sb->fhs);
	free(sb->fds);
	free(sb->running);
	free(sb->last_response);
	free(sb->pids);
	free(sb->tmpdir);

	STATUS("Final: %i images processed, %i had crystals (%.1f%%),"
	       " %i crystals overall.\n",
	       sb->shared->n_processed, sb->shared->n_hadcrystals,
	       100.0 * sb->shared->n_hadcrystals / sb->shared->n_processed,
	       sb->shared->n_crystals);

	munmap(sb->shared, sizeof(struct sb_shm));
	free(sb);
}
