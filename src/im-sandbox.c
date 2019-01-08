/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2019 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2019 Thomas White <taw@physics.org>
 *   2014      Valerio Mariani
 *   2011      Richard Kirian
 *   2012      Lorenzo Galli
 *   2012      Chunhong Yoon
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
#include "im-zmq.h"


struct sandbox
{
	int n_processed_last_stats;
	time_t t_last_stats;

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

	const char *tmpdir;

	/* ZMQ mode */
	int zmq;
	const char *zmq_address;

	/* Final output */
	Stream *stream;
};

struct get_pattern_ctx
{
	FILE *fh;
	int use_basename;
	struct detector *det;
	const char *prefix;
	char *filename;
	struct event_list *events;
	int event_index;
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
			STATUS("Worker %i did not respond for 4 minutes - "
			       "sending it SIGKILL.\n", i);
			kill(sb->pids[i], SIGKILL);
			stamp_response(sb, i);
		}

		if ( tnow - sb->shared->time_last_start[i] > 600 ) {
			if ( !sb->shared->warned_long_running[i] ) {
				STATUS("Worker %i has been working on one "
				       "frame for more than 10 minutes.\n", i);
				STATUS("Event ID is: %s\n",
				       sb->shared->last_ev[i]);
				STATUS("Task ID is: %s\n",
				       sb->shared->last_task[i]);
				sb->shared->warned_long_running[i] = 1;
			}
		}

	}
}


static char *read_prefixed_filename(struct get_pattern_ctx *gpctx, char **event)
{
	char* line;

	if ( event != NULL ) *event = NULL;

	line = malloc(1024);
	if ( line == NULL ) return NULL;

	do {
		if ( fgets(line, 1023, gpctx->fh) == NULL )
		{
			if ( !feof(gpctx->fh) ) {
				ERROR("Input file read error.\n");
			}
			free(line);
			return NULL;
		}
		chomp(line);

	} while ( line[0] == '\0' );

	/* Chop off event ID if requested */
	if ( event != NULL ) {
		size_t n = strlen(line);
		while ( line[n] != ' ' && n > 2 ) n--;
		if ( n != 2 ) {
			/* Event descriptor must contain "//".
			* If it doesn't, assume the filename just contains a
			* space. */
			if ( strstr(&line[n], "//") != NULL ) {
				line[n] = '\0';
				*event = strdup(&line[n+1]);
			}
		} /* else no spaces at all */
	}

	if ( gpctx->use_basename ) {
		char *tmp;
		tmp = safe_basename(line);
		free(line);
		line = tmp;
	}

	/* Add prefix */
	if ( gpctx->prefix != NULL ) {
		char *tmp;
		size_t len = strlen(line) + strlen(gpctx->prefix) + 1;
		tmp = malloc(len);
		if ( tmp == NULL ) {
			ERROR("Couldn't allocate memory for filename\n");
			return NULL;
		}
		strcpy(tmp, gpctx->prefix);
		strcat(tmp, line);
		free(line);
		line = tmp;
	}

	return line;
}


static struct filename_plus_event *get_pattern(struct get_pattern_ctx *gpctx)
{
	char *filename;
	char *evstr;

	/* If single-event geometry, just return next filename */
	if ( !multi_event_geometry(gpctx->det) )
	{
		struct filename_plus_event *fne;
		fne = malloc(sizeof(struct filename_plus_event));
		if ( fne == NULL ) return NULL;
		fne->filename = read_prefixed_filename(gpctx, NULL);
		if ( fne->filename == NULL ) {
			free(fne);
			return NULL;
		}
		fne->ev = NULL;
		return fne;
	}

	/* Ok, multi-event geometry.  Is an event available already? */
	if ( (gpctx->events != NULL)
	  && (gpctx->event_index < gpctx->events->num_events) )
	{
		struct filename_plus_event *fne;
		fne = malloc(sizeof(struct filename_plus_event));
		fne->filename = strdup(gpctx->filename);
		fne->ev = copy_event(gpctx->events->events[gpctx->event_index++]);

		return fne;
	}

	/* No events in list.  Time to top it up */
	filename = read_prefixed_filename(gpctx, &evstr);

	/* Nothing left in file -> we're done */
	if ( filename == NULL ) return NULL;

	/* Muppet check */
	if ( is_cbf_file(filename) == 1 ) {
		ERROR("Your geometry file is for a multi-event format, but "
		      "this file is in CBF format.\n");
		ERROR("Your geometry file probably needs to be changed.\n");
		return NULL;
	}

	/* Does the line from the input file contain an event ID?
	 * If so, just parse it and sent it straight back. */
	if ( evstr != NULL ) {

		/* Make an event list with only one item */
		struct event *ev = get_event_from_event_string(evstr);
		if ( ev == NULL ) {
			ERROR("Bad event descriptor: '%s'\n", evstr);
			return NULL;
		} else {
			struct filename_plus_event *fne;
			fne = malloc(sizeof(struct filename_plus_event));
			fne->filename = filename;
			fne->ev = ev;
			return fne;
		}

	} else {

		/* Enumerate all the events in the file and then send the
		 * first one back */
		struct hdfile *hdfile;
		struct filename_plus_event *fne;

		hdfile = hdfile_open(filename);
		if ( hdfile == NULL ) {
			ERROR("Failed to open %s\n", filename);
			return NULL;
		}

		if ( gpctx->events != NULL ) {
			free_event_list(gpctx->events);
		}

		gpctx->events = fill_event_list(hdfile, gpctx->det);
		if ( gpctx->events == NULL ) {
			ERROR("Failed to get event list.\n");
			return NULL;
		}

		hdfile_close(hdfile);

		/* Save filename for next time */
		if ( gpctx->filename != NULL ) {
			free(gpctx->filename);
		}
		gpctx->filename = filename;

		gpctx->event_index = 0;
		fne = malloc(sizeof(struct filename_plus_event));
		fne->filename = strdup(gpctx->filename);
		fne->ev = copy_event(gpctx->events->events[gpctx->event_index++]);
		return fne;

	}
}


static void shuffle_events(struct sb_shm *sb_shared)
{
	int i;

	for ( i=1; i<sb_shared->n_events; i++ ) {
		memcpy(sb_shared->queue[i-1], sb_shared->queue[i], MAX_EV_LEN);
	}
	sb_shared->n_events--;
}


void set_last_task(char *lt, const char *task)
{
	if ( lt == NULL ) return;
	assert(strlen(task) < MAX_TASK_LEN-1);
	strcpy(lt, task);
}


static int run_work(const struct index_args *iargs, Stream *st,
                    int cookie, const char *tmpdir, struct sandbox *sb)
{
	int allDone = 0;
	TimeAccounts *taccs;
	struct im_zmq *zmqstuff;

	/* Connect via ZMQ */
	if ( sb->zmq ) {
		zmqstuff = im_zmq_connect(sb->zmq_address);
		if ( zmqstuff == NULL ) {
			ERROR("ZMQ setup failed.\n");
			return 1;
		}
	}

	taccs = time_accounts_init();

	while ( !allDone ) {

		struct pattern_args pargs;
		char filename[MAX_EV_LEN];
		char event_str[MAX_EV_LEN];
		int ser;
		struct event *ev;
		int r;

		if ( !sb->zmq ) {

			/* Wait until an event is ready */
			time_accounts_set(taccs, TACC_EVENTWAIT);
			set_last_task(sb->shared->last_task[cookie], "wait_event");
			if ( sem_wait(sb->queue_sem) != 0 ) {
				ERROR("Failed to wait on queue semaphore: %s\n",
				      strerror(errno));
			}

			/* Get the event from the queue */
			set_last_task(sb->shared->last_task[cookie], "read_queue");
			pthread_mutex_lock(&sb->shared->queue_lock);
			if ( (sb->shared->n_events==0) && (sb->shared->no_more) ) {
				/* Queue is empty and no more coming, so exit */
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

			if ( strcmp(event_str, "(none)") != 0 ) {

				ev = get_event_from_event_string(event_str);
				if ( ev == NULL ) {
					ERROR("Bad event string '%s'\n", event_str);
					continue;
				}
				pargs.filename_p_e->ev = ev;

			} else {

				pargs.filename_p_e->ev = NULL;

			}

		} else {

			pargs.msgpack_obj = im_zmq_fetch(zmqstuff);

		}

		sb->shared->time_last_start[cookie] = get_monotonic_seconds();
		process_image(iargs, &pargs, st, cookie, tmpdir, ser,
		              sb->shared, taccs, sb->shared->last_task[cookie]);

		if ( !sb->zmq ) {
			free_filename_plus_event(pargs.filename_p_e);
		} else {
			im_zmq_clean(zmqstuff);
		}

	}

	im_zmq_shutdown(zmqstuff);

	time_accounts_set(taccs, TACC_FINALCLEANUP);
	cleanup_indexing(iargs->ipriv);
	free_detector_geometry(iargs->det);
	free(iargs->hdf5_peak_path);
	free_imagefile_field_list(iargs->copyme);
	cell_free(iargs->cell);
	if ( iargs->profile ) time_accounts_print(taccs);
	time_accounts_free(taccs);
	return 0;
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
	sb->shared->time_last_start[slot] = get_monotonic_seconds();
	sb->shared->warned_long_running[slot] = 0;

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

			r = mkdir(tmp, S_IRWXU);
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
		free(sb->running);
		/* Not freed because it's not worth passing them down just for
		 * this purpose: event list file handle,
		 *               main output stream handle
		 *               original temp dir name (without indexamajig.XX)
		 *               prefix
		 */

		st = open_stream_fd_for_write(stream_pipe[1]);
		r = run_work(sb->iargs, st, slot, tmp, sb);
		close_stream(st);

		free(tmp);
		free(sb->iargs->beam->photon_energy_from);

		munmap(sb->shared, sizeof(struct sb_shm));

		free(sb);

		exit(r);

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
				STATUS("Task ID was: %s\n",
				       sb->shared->last_task[i]);
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


/* Assumes the caller is already holding queue_lock! */
static int fill_queue(struct get_pattern_ctx *gpctx, struct sandbox *sb)
{
	if ( sb->zmq ) {
		/* Do nothing */
		return 0;
	}

	while ( sb->shared->n_events < QUEUE_SIZE ) {

		struct filename_plus_event *ne;

		ne = get_pattern(gpctx);
		if ( ne == NULL ) return 1; /* No more */

		memset(sb->shared->queue[sb->shared->n_events], 0, MAX_EV_LEN);
		snprintf(sb->shared->queue[sb->shared->n_events++], MAX_EV_LEN,
		         "%s %s %i", ne->filename, get_event_string(ne->ev),
		         sb->serial++);
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


static void try_status(struct sandbox *sb, int final)
{
	int r;
	int n_proc_this;
	time_t tNow;
	time_t time_this;
	const char *finalstr;
	char persec[64];

	tNow = get_monotonic_seconds();
	time_this = tNow - sb->t_last_stats;
	if ( !final && (time_this < 5) ) return;

	n_proc_this = sb->shared->n_processed - sb->n_processed_last_stats;

	r = pthread_mutex_trylock(&sb->shared->term_lock);
	if ( r ) return; /* No lock -> don't bother */

	finalstr = final ? "Final: " : "";
	if ( final ) {
		finalstr = "Final: ";
		persec[0] = '\0';
	} else {
		finalstr = "";
		snprintf(persec, 64, ", %.1f images/sec",
		         (double)n_proc_this/time_this);
	}
	STATUS("%s%i images processed, %i hits (%.1f%%), "
	       "%i indexable (%.1f%% of hits, %.1f%% overall), "
	       "%i crystals%s.\n",
	       finalstr, sb->shared->n_processed,
	       sb->shared->n_hits,
	       100.0 * sb->shared->n_hits / sb->shared->n_processed,
	       sb->shared->n_hadcrystals,
	       100.0 * sb->shared->n_hadcrystals / sb->shared->n_hits,
	       100.0 * sb->shared->n_hadcrystals / sb->shared->n_processed,
	       sb->shared->n_crystals, persec);

	sb->n_processed_last_stats = sb->shared->n_processed;
	sb->t_last_stats = tNow;

	pthread_mutex_unlock(&sb->shared->term_lock);
}


static void delete_temporary_folder(const char *tmpdir, int n_proc)
{
	int slot;
	size_t len;
	char *workerdir;
	char *path;

	/* List of files which it's safe to delete */
	char *files[] = {"gmon.out", "mosflm.lp", "SUMMARY", "XDS.INP",
	                 "xfel_001.img", "xfel_001.spt", "xfel.drx",
	                 "xfel.felix", "xfel.gve", "xfel.ini", "xfel.log",
	                 "IDXREF.LP", "SPOT.XDS", "xfel.newmat", "XPARM.XDS"};

	/* Number of items in the above list */
	int n_files = 15;

	if ( n_proc > 99999 ) return;  /* Paranoia */

	len = strlen(tmpdir);
	workerdir = calloc(len+32, 1);
	path = calloc(len+64, 1);
	if ( (workerdir == NULL) || (path == NULL) ) return;

	snprintf(path, 127, "%s/mosflm.lp", tmpdir);
	unlink(path);
	snprintf(path, 127, "%s/SUMMARY", tmpdir);
	unlink(path);

	for ( slot=0; slot<n_proc; slot++ ) {

		struct stat s;
		int i;

		snprintf(workerdir, 63, "%s/worker.%i", tmpdir, slot);
		if ( stat(workerdir, &s) == -1 ) continue;

		for ( i=0; i<n_files; i++ ) {
			snprintf(path, 127, "%s/%s", workerdir, files[i]);
			unlink(path);
		}

		if ( rmdir(workerdir) ) {
			ERROR("Failed to delete worker temporary folder: %s\n",
			      strerror(errno));
		}

	}

	if ( rmdir(tmpdir) ) {
		ERROR("Failed to delete temporary folder: %s\n", strerror(errno));
	}

	free(workerdir);
	free(path);
}


char *create_tempdir(const char *temp_location)
{
	char *tmpdir;
	size_t ll;
	struct stat s;

	if ( temp_location == NULL ) {
		temp_location = "";
	}

	ll = 64+strlen(temp_location);
	tmpdir = malloc(ll);
	if ( tmpdir == NULL ) {
		ERROR("Failed to allocate temporary directory name\n");
		return NULL;
	}
	snprintf(tmpdir, ll, "%s/indexamajig.%i", temp_location, getpid());

	if ( stat(tmpdir, &s) == -1 ) {

		int r;

		if ( errno != ENOENT ) {
			ERROR("Failed to stat temporary folder.\n");
			return NULL;
		}

		r = mkdir(tmpdir, S_IRWXU);
		if ( r ) {
			ERROR("Failed to create temporary folder: %s\n",
			      strerror(errno));
			return NULL;
		}

	}

	return tmpdir;
}


/* Returns the number of frames processed (not necessarily indexed).
 * If the return value is zero, something is probably wrong. */
int create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                   int config_basename, FILE *fh,
                   Stream *stream, const char *tmpdir, int serial_start,
                   const char *zmq_address)
{
	int i;
	struct sandbox *sb;
	char semname_q[64];
	struct sigaction sa;
	int r;
	int allDone = 0;
	TimeAccounts *taccs;
	struct get_pattern_ctx gpctx;

	if ( n_proc > MAX_NUM_WORKERS ) {
		ERROR("Number of workers (%i) is too large.  Using %i\n",
		      n_proc, MAX_NUM_WORKERS);
		n_proc = MAX_NUM_WORKERS;
	}

	sb = calloc(1, sizeof(struct sandbox));
	if ( sb == NULL ) {
		ERROR("Couldn't allocate memory for sandbox.\n");
		return 0;
	}

	sb->n_processed_last_stats = 0;
	sb->t_last_stats = get_monotonic_seconds();
	sb->n_proc = n_proc;
	sb->iargs = iargs;
	sb->serial = serial_start;
	sb->tmpdir = tmpdir;
	if ( zmq_address != NULL ) {
		sb->zmq = 1;
		sb->zmq_address = zmq_address;
	} else {
		sb->zmq = 0;
	}

	sb->fds = NULL;
	sb->fhs = NULL;
	sb->stream = stream;

	gpctx.fh = fh;
	gpctx.use_basename = config_basename;
	gpctx.det = iargs->det;
	gpctx.prefix = prefix;
	gpctx.filename = NULL;
	gpctx.events = NULL;
	gpctx.event_index = 0;

	if ( setup_shm(sb) ) {
		ERROR("Failed to set up SHM.\n");
		free(sb);
		return 0;
	}

	sb->shared->n_processed = 0;
	sb->shared->n_hits = 0;
	sb->shared->n_hadcrystals = 0;
	sb->shared->n_crystals = 0;

	/* Set up semaphore to control work queue */
	snprintf(semname_q, 64, "indexamajig-q%i", getpid());
	sb->queue_sem = sem_open(semname_q, O_CREAT | O_EXCL,
	                         S_IRUSR | S_IWUSR, 0);
	if ( sb->queue_sem == SEM_FAILED ) {
		ERROR("Failed to create semaphore: %s\n", strerror(errno));
		return 0;
	}

	sb->pids = calloc(n_proc, sizeof(pid_t));
	sb->running = calloc(n_proc, sizeof(int));
	sb->last_response = calloc(n_proc, sizeof(time_t));
	if ( (sb->pids == NULL) || (sb->running == NULL)
	  || (sb->last_response == NULL) )
	{
		ERROR("Couldn't allocate memory for PIDs.\n");
		return 0;
	}

	/* Fill the queue */
	pthread_mutex_lock(&sb->shared->queue_lock);
	sb->shared->n_events = 0;
	sb->shared->no_more = fill_queue(&gpctx, sb);
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
	        return 0;
	}

	/* Set up signal handler to clean up semaphore on exit */
	sa.sa_flags = SA_SIGINFO | SA_NOCLDSTOP | SA_RESTART;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = sigint_handler;
	r = sigaction(SIGINT, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return 0;
	}
	r = sigaction(SIGQUIT, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return 0;
	}

	taccs = time_accounts_init();

	if ( sb->zmq ) {
		allDone = 1;
	}

	do {

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
		if ( !sb->shared->no_more && (sb->shared->n_events < QUEUE_SIZE/2) ) {
			if ( fill_queue(&gpctx, sb) ) sb->shared->no_more = 1;
		}
		pthread_mutex_unlock(&sb->shared->queue_lock);

		/* Update progress */
		time_accounts_set(taccs, TACC_STATUS);
		try_status(sb, 0);

		/* Have all the events been swallowed? */
		time_accounts_set(taccs, TACC_ENDCHECK);
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( sb->shared->no_more && (sb->shared->n_events == 0) ) allDone = 1;
		pthread_mutex_unlock(&sb->shared->queue_lock);

	} while ( !allDone );

	fclose(fh);

	/* Indicate to the workers that we are finished, and wake them up one
	 * last time */
	time_accounts_set(taccs, TACC_WAKEUP);
	if ( sb->zmq ) {
		STATUS("Waiting for patterns to process...\n");
	} else {
		STATUS("Waiting for the last patterns to be processed...\n");
	}
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

			time_accounts_set(taccs, TACC_STATUS);
			try_status(sb, 0);

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

	try_status(sb, 1);
	r = sb->shared->n_processed;

	delete_temporary_folder(sb->tmpdir, n_proc);

	munmap(sb->shared, sizeof(struct sb_shm));
	free(sb);

	return r;
}
