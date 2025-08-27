/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2024 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2024 Thomas White <taw@physics.org>
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

#ifdef HAVE_SCHED_SETAFFINITY
#define _GNU_SOURCE
#include <sys/sysinfo.h>
#include <sched.h>
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
#include <poll.h>

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
#else
#include <sys/time.h>
#endif

#include "im-sandbox.h"
#include "im-argparse.h"
#include "process_image.h"
#include "im-zmq.h"
#include "profile.h"
#include "im-asapo.h"
#include "predict-refine.h"


typedef struct
{
	int n_read;
	int *fds;
	char **buffers;
	size_t *buffer_len;
	size_t *buffer_pos;
} PipeList;


struct sandbox
{
	int n_processed_last_stats;
	time_t t_last_stats;

	/* Processing timeout in seconds.  After this long without responding
	 * to a ping, the worker will be killed.  After 3 times this long
	 * working on one image, even with ping responses, a warning will be
	 * shown to the user. */
	int timeout;

	struct index_args *iargs;
	int argc;
	char **argv;
	const char *probed_methods;

	/* Worker processes */
	int n_proc;
	pid_t *pids;
	int *running;
	time_t *last_response;
	int last_ping[MAX_NUM_WORKERS];
	int warned_long_running[MAX_NUM_WORKERS];
	int profile;  /* Whether to do wall-clock time profiling */
	int cpu_pin;

	/* Streams to read from (NB not the same indices as the above) */
	PipeList *st_from_workers;
	PipeList *mille_from_workers;

	int serial;

	struct sb_shm *shared;
	char *shm_name;
	sem_t *queue_sem;
	char *sem_name;

	const char *tmpdir;

	/* If non-NULL, we are using ZMQ */
	struct im_zmq_params *zmq_params;

	/* If non-NULL, we are using ASAP::O */
	struct im_asapo_params *asapo_params;

	/* Final output */
	Stream *stream;
	FILE *mille_fh;
};

struct get_pattern_ctx
{
	FILE *fh;
	int use_basename;
	const DataTemplate *dtempl;
	const char *prefix;
	char *filename;
	char **events;
	int n_events;
	int event_index;
};


#ifdef HAVE_CLOCK_GETTIME

time_t get_monotonic_seconds()
{
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);
	return tp.tv_sec;
}

#else

/* Fallback version of the above.  The time according to gettimeofday() is not
 * monotonic, so measuring intervals based on it will screw up if there's a
 * timezone change (e.g. daylight savings) while the program is running. */
time_t get_monotonic_seconds()
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

	pthread_mutex_lock(&sb->shared->debug_lock);
	for ( i=0; i<sb->n_proc; i++ ) {

		if ( !sb->running[i] ) continue;

		if ( sb->shared->pings[i] != sb->last_ping[i] ) {
			stamp_response(sb, i);
		}

		if ( tnow - sb->last_response[i] > sb->timeout ) {
			STATUS("Worker %i did not respond for %i seconds - "
			       "sending it SIGKILL.\n", i, sb->timeout);
			kill(sb->pids[i], SIGKILL);
			stamp_response(sb, i);
		}

		if ( tnow - sb->shared->time_last_start[i] > sb->timeout*3 ) {
			if ( !sb->warned_long_running[i] ) {
				STATUS("Worker %i has been working on one "
				       "frame for more than %i seconds (just "
				       "for info).\n", i, sb->timeout);
				STATUS("Event ID is: %s\n",
				       sb->shared->last_ev[i]);
				STATUS("Task ID is: %s\n",
				       sb->shared->last_task[i]);
				sb->warned_long_running[i] = 1;
			}
		}

	}
	pthread_mutex_unlock(&sb->shared->debug_lock);
}


static char *read_prefixed_filename(struct get_pattern_ctx *gpctx,
                                    char **event)
{
	char* line;

	*event = NULL;

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

	/* Chop off event ID */
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


/* Return 0 for "no more" */
static int get_pattern(struct get_pattern_ctx *gpctx,
                       char **pfilename, char **pevent)
{
	char *filename;
	char *evstr;

	/* Is an event available already? */
	if ( (gpctx->events != NULL)
	  && (gpctx->event_index < gpctx->n_events) )
	{
		*pfilename = gpctx->filename;
		*pevent = gpctx->events[gpctx->event_index++];
		return 1;
	}

	do {

		/* No events in list.  Time to top it up */
		filename = read_prefixed_filename(gpctx, &evstr);

		/* Nothing left in file -> we're done */
		if ( filename == NULL ) {
			free(gpctx->filename);
			free(gpctx->events);
			return 0;
		}

		/* Does the line from the input file contain an event ID?
		 * If so, just send it straight back. */
		if ( evstr != NULL ) {
			*pfilename = filename;
			*pevent = evstr;
			free(gpctx->filename);
			gpctx->filename = filename;
			return 1;
		}

		/* We got a filename, but no event.  Attempt to expand... */
		free(gpctx->events);  /* Free the old list.
		                       * NB The actual strings were freed
		                       * by fill_queue */
		gpctx->events = image_expand_frames(gpctx->dtempl, filename,
		                                    &gpctx->n_events);
		if ( gpctx->events == NULL ) {
			ERROR("Failed to get event list from %s.\n",
			      filename);
		}

	} while ( gpctx->events == NULL );

	/* Save filename for next time */
	free(gpctx->filename);
	gpctx->filename = filename;

	gpctx->event_index = 0;
	*pfilename = gpctx->filename;
	*pevent = gpctx->events[gpctx->event_index++];
	return 1;
}


/* Like strstr(), but with a length parameter instead of nul-termination */
static const char *str_in_str(const char *haystack, size_t len, const char *needle)
{
	size_t u;
	const size_t endl = strlen(needle);

	if ( len < endl ) return NULL;

	for ( u=0; u<=len-endl; u++ ) {
		size_t v;
		int ok = 1;
		for ( v=0; v<endl; v++ ) {
			if ( haystack[u+v] != needle[v] ) {
				ok = 0;
				break;
			}
		}
		if ( ok ) {
			return &haystack[u];
		}
	}
	return NULL;
}


static size_t pump_chunk(char *buf, size_t len, struct sandbox *sb)
{
	const char *txt = (char *)buf;
	const char *endpos;
	size_t chunk_len;

	endpos = str_in_str(txt, len, STREAM_CHUNK_END_MARKER"\n");
	if ( endpos == NULL ) return 0;

	chunk_len = (endpos-txt)+strlen(STREAM_CHUNK_END_MARKER"\n");
	fwrite(buf, 1, chunk_len, stream_get_fh(sb->stream));
	fflush(stream_get_fh(sb->stream));

	return chunk_len;
}


static size_t pump_mille(char *buf, size_t len, struct sandbox *sb)
{
	int n;
	int ni;
	size_t reclen;

	if ( len < sizeof(int) ) return 0;
	ni = *(int *)buf;
	n = (ni-2)/2;

	reclen = (2+n)*sizeof(int) + (1+n)*sizeof(float);
	if ( len < reclen ) return 0;

	fwrite(buf, 1, reclen, sb->mille_fh);
	fflush(sb->mille_fh);

	return reclen;
}


static PipeList *pipe_list_new()
{
	PipeList *pd = malloc(sizeof(PipeList));
	if ( pd == NULL ) return NULL;
	pd->fds = NULL;
	pd->buffers = NULL;
	pd->buffer_len = NULL;
	pd->buffer_pos = NULL;
	pd->n_read = 0;
	return pd;
}


static void pipe_list_destroy(PipeList *pd)
{
	int i;
	for ( i=0; i<pd->n_read; i++ ) {
		free(pd->buffers[i]);
	}
	free(pd->buffers);
	free(pd->buffer_len);
	free(pd->buffer_pos);
	free(pd->fds);
	free(pd);
}


/* Add an fd to the list of pipes to be read from */
static void add_pipe(PipeList *pd, int fd)
{
	int *fds_new;
	char **buffers_new;
	size_t *buflens_new;
	size_t *bufposs_new;
	int slot;

	fds_new = realloc(pd->fds, (pd->n_read+1)*sizeof(int));
	buffers_new = realloc(pd->buffers, (pd->n_read+1)*sizeof(char *));
	buflens_new = realloc(pd->buffer_len, (pd->n_read+1)*sizeof(size_t));
	bufposs_new = realloc(pd->buffer_pos, (pd->n_read+1)*sizeof(size_t));
	if ( (fds_new == NULL) || (buffers_new == NULL)
	  || (buflens_new == NULL) || (bufposs_new == NULL) )
	{
		ERROR("Failed to allocate memory for new pipe.\n");
		return;
	}

	pd->fds = fds_new;
	pd->buffers = buffers_new;
	pd->buffer_len = buflens_new;
	pd->buffer_pos = bufposs_new;

	slot = pd->n_read;
	pd->fds[slot] = fd;
	pd->buffer_len[slot] = 64*1024;  /* Initial buffer size */
	pd->buffers[slot] = malloc(pd->buffer_len[slot]);
	if ( pd->buffers[slot] == NULL ) return;
	pd->buffer_pos[slot] = 0;

	pd->n_read++;
}


static void remove_pipe(PipeList *pd, int d)
{
	int i;

	close(pd->fds[d]);

	if ( pd->buffer_pos[d] > 0 ) {
		ERROR("WARNING: Closing pipe with %lli bytes remaining in buffer.\n",
		      (long long)pd->buffer_pos[d]);
	}

	free(pd->buffers[d]);

	for ( i=d; i<pd->n_read; i++ ) {
		if ( i < pd->n_read-1 ) {
			pd->fds[i] = pd->fds[i+1];
			pd->buffers[i] = pd->buffers[i+1];
			pd->buffer_len[i] = pd->buffer_len[i+1];
			pd->buffer_pos[i] = pd->buffer_pos[i+1];
		} /* else don't bother */
	}

	pd->n_read--;

	/* We don't bother shrinking the arrays */
}


static int find_marked(PipeList *pd)
{
	int i;
	for ( i=0; i<pd->n_read; i++ ) {
		if ( pd->buffer_len[i] == 0 ) return i;
	}
	return -1;
}


static void check_pipes(PipeList *pd, size_t(*pump)(char *, size_t len, struct sandbox *),
                        struct sandbox *sb)
{
	int r, i;
	struct pollfd *pf;

	pf = malloc(pd->n_read*sizeof(struct pollfd));
	if ( pf == NULL ) {
		ERROR("Couldn't allocate poll structure!\n");
		return;
	}

	for ( i=0; i<pd->n_read; i++ ) {
		pf[i].fd = pd->fds[i];
		pf[i].events = POLLIN;
		pf[i].revents = 0;
	}

	r = poll(pf, pd->n_read, 500);

	if ( r == -1 ) {
		if ( errno != EINTR ) {
			ERROR("poll() failed: %s\n", strerror(errno));
			return;
		}
	}

	for ( i=0; i<pd->n_read; i++ ) {

		size_t r;

		if ( pf[i].revents == 0 ) continue;

		if ( pd->buffer_len[i] == pd->buffer_pos[i] ) {
			const size_t buffer_increment = 64*1024;
			char *buf_new = realloc(pd->buffers[i],
			                        pd->buffer_len[i]+buffer_increment);
			if ( buf_new == NULL ) {
				ERROR("Failed to grow buffer\n");
				continue;
			}
			pd->buffers[i] = buf_new;
			pd->buffer_len[i] += buffer_increment;
		}

		r = read(pd->fds[i], pd->buffers[i]+pd->buffer_pos[i],
		         pd->buffer_len[i]-pd->buffer_pos[i]);

		if ( r == -1 ) {
			ERROR("read() failed\n");
			continue;
		}

		/* End of file? */
		if ( r == 0 ) {
			pd->buffer_len[i] = 0;    /* Mark for deletion */
		}

		pd->buffer_pos[i] += r;

		size_t h;
		do {
			h = pump(pd->buffers[i], pd->buffer_pos[i], sb);
			assert(h <= pd->buffer_pos[i]);
			assert(h >= 0);
			if ( h > 0 ) {
				memmove(pd->buffers[i], pd->buffers[i]+h,
				        pd->buffer_pos[i]-h);
				pd->buffer_pos[i] -= h;
			}
		} while ( h > 0 );

	}

	free(pf);

	int deleteme = find_marked(pd);
	while ( deleteme != -1 ) {
		remove_pipe(pd, deleteme);
		deleteme = find_marked(pd);
	}
}


static void try_read(struct sandbox *sb)
{
	check_pipes(sb->st_from_workers, pump_chunk, sb);
	check_pipes(sb->mille_from_workers, pump_mille, sb);
}


static void start_worker_process(struct sandbox *sb, int slot)
{
	pid_t p;
	int stream_pipe[2];
	int mille_pipe[2];
	char **nargv;
	int nargc;
	int i;
	char *tmpdir;
	char *methods_copy = NULL;
	char *worker_id;
	char *fd_stream;
	char *fd_mille;
	char buf[1024];
	const char *indexamajig = NULL;
	size_t len;
	struct stat s;

	if ( pipe(stream_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	if ( pipe(mille_pipe) == - 1 ) {
		ERROR("pipe() failed!\n");
		return;
	}

	pthread_mutex_lock(&sb->shared->queue_lock);
	sb->shared->end_of_stream[slot] = 0;
	pthread_mutex_unlock(&sb->shared->queue_lock);

	pthread_mutex_lock(&sb->shared->debug_lock);
	sb->shared->pings[slot] = 0;
	sb->last_ping[slot] = 0;
	sb->shared->time_last_start[slot] = get_monotonic_seconds();
	pthread_mutex_unlock(&sb->shared->debug_lock);

	sb->warned_long_running[slot] = 0;

	/* Set up nargv including "new" args */
	nargc = 0;
	nargv = malloc((sb->argc+16)*sizeof(char *));
	if ( nargv == NULL ) return;
	for ( i=0; i<sb->argc; i++ ) {
		nargv[nargc++] = sb->argv[i];
	}

	nargv[nargc++] = "--shm-name";
	nargv[nargc++] = sb->shm_name;

	nargv[nargc++] = "--queue-sem";
	nargv[nargc++] = sb->sem_name;

	nargv[nargc++] = "--worker-tmpdir";
	len = 64 + strlen(sb->tmpdir);
	tmpdir = malloc(len);
	if ( tmpdir == NULL ) {
		ERROR("Failed to allocate temporary dir\n");
		return;
	}
	snprintf(tmpdir, len, "%s/worker.%i", sb->tmpdir, slot);
	if ( stat(tmpdir, &s) == -1 ) {
		int r;
		if ( errno != ENOENT ) {
			ERROR("Failed to stat temporary folder.\n");
			return;
		}
		r = mkdir(tmpdir, S_IRWXU);
		if ( r ) {
			ERROR("Failed to create temporary folder: %s\n",
			strerror(errno));
			return;
		}
	}
	nargv[nargc++] = tmpdir;

	nargv[nargc++] = "--worker-id";
	worker_id = malloc(64);
	if ( worker_id == NULL ) return;
	snprintf(worker_id, 64, "%i", slot);
	nargv[nargc++] = worker_id;

	nargv[nargc++] = "--fd-stream";
	fd_stream = malloc(64);
	snprintf(fd_stream, 64, "%i", stream_pipe[1]);
	nargv[nargc++] = fd_stream;

	nargv[nargc++] = "--fd-mille";
	fd_mille = malloc(64);
	snprintf(fd_mille, 64, "%i", mille_pipe[1]);
	nargv[nargc++] = fd_mille;

	if ( sb->probed_methods != NULL ) {
		methods_copy = strdup(sb->probed_methods);
		nargv[nargc++] = "--indexing";
		nargv[nargc++] = methods_copy;
	}
	nargv[nargc++] = NULL;

	len = readlink("/proc/self/exe", buf, 1024);

	if ( (len == -1) || (len >= 1023)  ) {
		ERROR("readlink() failed: %s\n", strerror(errno));
	} else {
		buf[len] = '\0';
		if ( strstr(buf, "indexamajig") == NULL ) {
			ERROR("Didn't recognise /proc/self/exe (%s)\n", buf);
		} else {
			indexamajig = buf;
		}
	}
	if ( indexamajig == NULL ) {
		if ( strstr(sb->argv[0], "indexamajig") == NULL ) {
			ERROR("Didn't recognise argv[0] (%s)\n", sb->argv[0]);
		} else {
			indexamajig = sb->argv[0];
		}
	}
	if ( indexamajig == NULL ) {
		ERROR("Falling back on shell search path.\n");
		indexamajig = "indexamajig";
	}

	p = fork();
	if ( p == -1 ) {
		ERROR("fork() failed!\n");
		return;
	}

	if ( p == 0 ) {
		execvp(indexamajig, nargv);
		ERROR("Failed to exec!\n");
		exit(1);
	}

	free(tmpdir);
	free(methods_copy);
	free(worker_id);
	free(fd_stream);
	free(fd_mille);
	free(nargv);

	/* Parent process gets the 'write' end of the filename pipe
	 * and the 'read' end of the result pipe. */
	sb->pids[slot] = p;
	sb->running[slot] = 1;
	pthread_mutex_lock(&sb->shared->debug_lock);
	stamp_response(sb, slot);
	pthread_mutex_unlock(&sb->shared->debug_lock);
	add_pipe(sb->st_from_workers, stream_pipe[0]);
	add_pipe(sb->mille_from_workers, mille_pipe[0]);
	close(stream_pipe[1]);
	close(mille_pipe[1]);
}


static int any_running(struct sandbox *sb)
{
	int i;
	for ( i=0; i<sb->n_proc; i++ ) {
		if ( sb->running[i] ) return 1;
	}
	return 0;
}


static void handle_zombie(struct sandbox *sb, int respawn)
{
	int i;

	for ( i=0; i<sb->n_proc; i++ ) {

		int status, p;

		if ( !sb->running[i] ) continue;

		p = waitpid(sb->pids[i], &status, WNOHANG);

		if ( p == -1 ) {
			ERROR("waitpid(%i) failed: %s.\n", i, strerror(errno));
			continue;
		}

		if ( p == sb->pids[i] ) {

			sb->running[i] = 0;

			if ( WIFEXITED(status) ) {
				if ( WEXITSTATUS(status) != 0 ) {
					STATUS("Worker %i returned error code %i\n",
					       i, WEXITSTATUS(status));
					STATUS("Shutting down.\n");
					/* Error status from worker */
					pthread_mutex_lock(&sb->shared->totals_lock);
					sb->shared->should_shutdown = 1;
					pthread_mutex_unlock(&sb->shared->totals_lock);
				}
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
	char tmp[128];
	int shm_fd;

	snprintf(tmp, 127, "/indexamajig.shm.%i", getpid());
	shm_fd = shm_open(tmp, O_CREAT | O_EXCL | O_RDWR, 0600);
	if ( shm_fd == -1 ) {
		ERROR("SHM setup failed: %s\n", strerror(errno));
		return 1;
	}
	sb->shm_name = strdup(tmp);

	if ( ftruncate(shm_fd, sizeof(struct sb_shm)) == -1 ) {
		ERROR("SHM setup failed: %s\n", strerror(errno));
		free(sb->shm_name);
		return 1;
	}

	sb->shared = mmap(NULL, sizeof(struct sb_shm), PROT_READ | PROT_WRITE,
	                  MAP_SHARED, shm_fd, 0);
	if ( sb->shared == MAP_FAILED ) {
		ERROR("SHM setup failed: %s\n", strerror(errno));
		free(sb->shm_name);
		return 1;
	}

	if ( pthread_mutexattr_init(&attr) ) {
		ERROR("Failed to initialise mutex attr.\n");
		free(sb->shm_name);
		return 1;
	}

	if ( pthread_mutexattr_setpshared(&attr, PTHREAD_PROCESS_SHARED) ) {
		ERROR("Failed to set process shared attribute.\n");
		pthread_mutexattr_destroy(&attr);
		free(sb->shm_name);
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->term_lock, &attr) ) {
		ERROR("Terminal lock setup failed.\n");
		pthread_mutexattr_destroy(&attr);
		free(sb->shm_name);
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->queue_lock, &attr) ) {
		ERROR("Queue lock setup failed.\n");
		pthread_mutexattr_destroy(&attr);
		pthread_mutex_destroy(&sb->shared->term_lock);
		free(sb->shm_name);
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->debug_lock, &attr) ) {
		ERROR("Queue lock setup failed.\n");
		pthread_mutexattr_destroy(&attr);
		pthread_mutex_destroy(&sb->shared->term_lock);
		pthread_mutex_destroy(&sb->shared->queue_lock);
		free(sb->shm_name);
		return 1;
	}

	if ( pthread_mutex_init(&sb->shared->totals_lock, &attr) ) {
		ERROR("Totals lock setup failed.\n");
		pthread_mutexattr_destroy(&attr);
		pthread_mutex_destroy(&sb->shared->term_lock);
		pthread_mutex_destroy(&sb->shared->queue_lock);
		pthread_mutex_destroy(&sb->shared->debug_lock);
		free(sb->shm_name);
		return 1;
	}

	pthread_mutexattr_destroy(&attr);

	return 0;
}


/* Assumes the caller is already holding queue_lock! */
static int fill_queue(struct get_pattern_ctx *gpctx, struct sandbox *sb)
{
	while ( sb->shared->n_events < QUEUE_SIZE ) {

		char *filename;
		char *evstr;

		if ( sb->zmq_params != NULL ) {
			/* These are just semi-meaningful placeholder values to
			 * be put into the queue, instead of "(null)".
			 * A unique filename is needed so that the GUI can
			 * tell the frames apart from one another.
			 * ASAP::O, for one, will replace this with a filename
			 * that corresponds to something real. */
			filename = "ZMQdata";
			evstr = malloc(64);
			snprintf(evstr, 64, "//%i", sb->serial);
		} else if ( sb->asapo_params != NULL ) {
			filename = "ASAPOdata";
			evstr = malloc(64);
			snprintf(evstr, 64, "//%i", sb->serial);
		} else {
			if ( !get_pattern(gpctx, &filename, &evstr) ) return 1;
		}

		memset(sb->shared->queue[sb->shared->n_events], 0, MAX_EV_LEN);
		snprintf(sb->shared->queue[sb->shared->n_events++], MAX_EV_LEN,
		         "%s %s %i", filename, evstr, sb->serial++);
		sem_post(sb->queue_sem);
		free(evstr);

	}
	return 0;
}

volatile sig_atomic_t at_zombies = 0;
volatile sig_atomic_t at_interrupt = 0;
volatile sig_atomic_t at_shutdown = 0;

static void sigchld_handler(int sig, siginfo_t *si, void *uc_v)
{
	at_zombies = 1;
}


static void sigint_handler(int sig, siginfo_t *si, void *uc_v)
{
	at_interrupt = 1;
}


static void sigusr1_handler(int sig, siginfo_t *si, void *uc_v)
{
	at_shutdown = 1;
}


static void check_signals(struct sandbox *sb, int respawn)
{
	if ( at_zombies ) {
		at_zombies = 0;
		handle_zombie(sb, respawn);
	}

	if ( at_interrupt ) {
		sem_unlink(sb->sem_name);
		shm_unlink(sb->shm_name);
		exit(0);
	}

	if ( at_shutdown ) {
		at_shutdown = 0;
		STATUS("Received signal - shutting down cleanly.\n");
		pthread_mutex_lock(&sb->shared->totals_lock);
		sb->shared->should_shutdown = 1;
		pthread_mutex_unlock(&sb->shared->totals_lock);
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

	pthread_mutex_lock(&sb->shared->totals_lock);
	n_proc_this = sb->shared->n_processed - sb->n_processed_last_stats;
	pthread_mutex_unlock(&sb->shared->totals_lock);

	r = pthread_mutex_trylock(&sb->shared->term_lock);
	if ( r ) return; /* No lock -> don't bother */

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
	size_t len, pathlen, workerdirlen;
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
	workerdirlen = len+32;
	workerdir = calloc(workerdirlen, 1);
	pathlen = len+64;
	path = calloc(pathlen, 1);

	if ( (workerdir == NULL) || (path == NULL) ) return;

	snprintf(path, pathlen, "%s/mosflm.lp", tmpdir);
	unlink(path);
	snprintf(path, pathlen, "%s/SUMMARY", tmpdir);
	unlink(path);

	for ( slot=0; slot<n_proc; slot++ ) {

		struct stat s;
		int i;

		snprintf(workerdir, workerdirlen, "%s/worker.%i", tmpdir, slot);
		if ( stat(workerdir, &s) == -1 ) continue;

		for ( i=0; i<n_files; i++ ) {
			snprintf(path, pathlen, "%s/%s", workerdir, files[i]);
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


/* Call under queue_lock */
static int all_got_end_of_stream(struct sandbox *sb)
{
	int i;
	for ( i=0; i<sb->n_proc; i++ ) {
		if ( !sb->shared->end_of_stream[i] ) return 0;
	}
	return 1;
}


static int any_open(struct sandbox *sb)
{
	return (sb->st_from_workers->n_read > 0) || (sb->mille_from_workers->n_read > 0);
}


/* Returns the number of frames processed (not necessarily indexed).
 * If the return value is zero, something is probably wrong. */
int create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                   int config_basename, FILE *fh,
                   Stream *stream, const char *tmpdir, int serial_start,
                   struct im_zmq_params *zmq_params,
                   struct im_asapo_params *asapo_params,
                   int timeout, int profile, int cpu_pin,
                   int no_data_timeout, int argc, char *argv[],
                   const char *probed_methods, FILE *mille_fh)
{
	int i;
	struct sandbox *sb;
	char semname_q[64];
	struct sigaction sa;
	int r;
	int allDone = 0;
	struct get_pattern_ctx gpctx;
	time_t t_last_data;

	if ( n_proc > MAX_NUM_WORKERS ) {
		ERROR("Number of workers (%i) is too large.  Using %i\n",
		      n_proc, MAX_NUM_WORKERS);
		n_proc = MAX_NUM_WORKERS;
	}

	#ifdef HAVE_SCHED_SETAFFINITY
	int n_cpus = get_nprocs();
	if ( n_proc > n_cpus ) {
		ERROR("WARNING: Number of workers (%i) is larger than the "
		      "number of available CPUs (%i)\n", n_proc, n_cpus);
		if ( cpu_pin ) {
			ERROR("Try again with a smaller number of workers (-j) "
			      "or without --cpu-pin\n");
			return 1;
		}
	}
	#else
	if ( cpu_pin ) {
		ERROR("Option --cpu-pin not available on this system.");
		return 1;
	}
	#endif

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
	sb->profile = profile;
	sb->timeout = timeout;
	sb->argc = argc;
	sb->argv = argv;
	sb->probed_methods = probed_methods;
	sb->mille_fh = mille_fh;

	if ( zmq_params->addr != NULL ) {
		sb->zmq_params = zmq_params;
	} else {
		sb->zmq_params = NULL;
	}

	if ( asapo_params->endpoint != NULL ) {
		sb->asapo_params = asapo_params;
	} else {
		sb->asapo_params = NULL;
	}

	if ( sb->zmq_params && sb->asapo_params ) {
		ERROR("Cannot simultaneously use ZMQ and ASAP::O input.\n");
		free(sb);
		return 0;
	}

	sb->st_from_workers = pipe_list_new();
	sb->mille_from_workers = pipe_list_new();
	sb->stream = stream;

	gpctx.fh = fh;
	gpctx.use_basename = config_basename;
	gpctx.dtempl = iargs->dtempl;
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
	sb->shared->should_shutdown = 0;

	/* Set up semaphore to control work queue */
	snprintf(semname_q, 64, "indexamajig-q%i", getpid());
	sb->queue_sem = sem_open(semname_q, O_CREAT | O_EXCL,
	                         S_IRUSR | S_IWUSR, 0);
	if ( sb->queue_sem == SEM_FAILED ) {
		ERROR("Failed to create semaphore: %s\n", strerror(errno));
		return 0;
	}
	sb->sem_name = strdup(semname_q);

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

	/* Set up signal handler to shut down gracefully on request */
	sa.sa_flags = SA_SIGINFO | SA_NOCLDSTOP | SA_RESTART;
	sigemptyset(&sa.sa_mask);
	sa.sa_sigaction = sigusr1_handler;
	r = sigaction(SIGUSR1, &sa, NULL);
	if ( r == -1 ) {
	        ERROR("Failed to set signal handler!\n");
	        return 0;
	}

	t_last_data = get_monotonic_seconds();
	do {

		/* Check for stream output from workers */
		try_read(sb);

		/* Check for interrupt or zombies */
		check_signals(sb, 1);

		/* Check for hung workers */
		check_hung_workers(sb);

		/* Top up the queue if necessary */
		pthread_mutex_lock(&sb->shared->queue_lock);
		if ( !sb->shared->no_more && (sb->shared->n_events < QUEUE_SIZE/2) ) {
			if ( fill_queue(&gpctx, sb) ) sb->shared->no_more = 1;
		}
		pthread_mutex_unlock(&sb->shared->queue_lock);

		/* Update progress */
		try_status(sb, 0);

		/* Begin exit criterion checking */
		pthread_mutex_lock(&sb->shared->queue_lock);
		pthread_mutex_lock(&sb->shared->totals_lock);

		/* Case 1: Queue empty and no more coming? */
		if ( sb->shared->no_more && (sb->shared->n_events == 0) ) allDone = 1;

		/* Case 2: Worker process requested immediate shutdown */
		if ( sb->shared->should_shutdown ) {
			allDone = 1;
			sb->shared->n_events = 0;
			sb->shared->no_more = 1;
		}

		/* Case 3: No (ASAP::O) data for a long time */
		if ( get_monotonic_seconds() > t_last_data + no_data_timeout ) {
			allDone = 1;
		}
		if ( !all_got_end_of_stream(sb) ) {
			/* We are still getting data */
			t_last_data = get_monotonic_seconds();
		}

		pthread_mutex_unlock(&sb->shared->totals_lock);
		pthread_mutex_unlock(&sb->shared->queue_lock);
		/* End exit criterion checking */

	} while ( !allDone );

	if ( fh != NULL ) {
		fclose(fh);
	}

	/* Indicate to the workers that we are finished, and wake them up one
	 * last time */
	STATUS("Waiting for the last patterns to be processed...\n");
	pthread_mutex_lock(&sb->shared->queue_lock);
	sb->shared->no_more = 1;
	pthread_mutex_unlock(&sb->shared->queue_lock);
	for ( i=0; i<n_proc; i++ ) {
		sem_post(sb->queue_sem);
	}
	for ( i=0; i<n_proc; i++ ) {
		while ( any_running(sb) ) {
			try_read(sb);
			check_signals(sb, 0);
			check_hung_workers(sb);
			try_status(sb, 0);
		}
		/* If this worker died and got waited by the zombie handler,
		 * waitpid() returns -1 and the loop still exits. */
	}

	if ( any_open(sb) ) {
		STATUS("Waiting for data pipes to close...\n");
		do {
			try_read(sb);
		} while ( any_open(sb) );
	}

	sem_unlink(semname_q);
	sem_close(sb->queue_sem);

	pipe_list_destroy(sb->st_from_workers);
	pipe_list_destroy(sb->mille_from_workers);
	free(sb->running);
	free(sb->last_response);
	free(sb->pids);

	try_status(sb, 1);
	if ( sb->shared->n_processed == 0 ) r = 5;
	if ( sb->shared->should_shutdown ) r = 1;

	delete_temporary_folder(sb->tmpdir, n_proc);

	shm_unlink(sb->shm_name);
	munmap(sb->shared, sizeof(struct sb_shm));
	free(sb->shm_name);
	free(sb->sem_name);
	free(sb);

	return r;
}
