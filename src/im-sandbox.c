/*
 * im-sandbox.c
 *
 * Sandbox for indexing
 *
 * Copyright © 2012-2014 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 * Copyright © 2012 Richard Kirian
 * Copyright © 2012 Lorenzo Galli
 *
 * Authors:
 *   2010-2014 Thomas White <taw@physics.org>
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
	pthread_mutex_t lock;

	int n_processed;
	int n_hadcrystals;
	int n_crystals;
	int n_processed_last_stats;
	int n_hadcrystals_last_stats;
	int n_crystals_last_stats;
	int t_last_stats;
	int suspend_stats;

	struct index_args *iargs;

	int n_proc;
	pid_t *pids;

	int *running;
	FILE **result_fhs;
	int *filename_pipes;
	int *stream_pipe_write;
	struct filename_plus_event **last_filename;
	int serial;

	char *tmpdir;

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


struct buffer_data
{
	char *rbuffer;
	char *line;
	int fd;
	int rbufpos;
	int rbuflen;
};


static int read_fpe_data(struct buffer_data *bd)
{
	int rval;
	int no_line = 0;

	rval = read(bd->fd, bd->rbuffer+bd->rbufpos, bd->rbuflen-bd->rbufpos);
	if ( (rval == -1) || (rval == 0) ) return 1;
	bd->rbufpos += rval;
	assert(bd->rbufpos <= bd->rbuflen);

	while ( (!no_line) && (bd->rbufpos > 0) ) {

		int i;
		int line_ready = 0;
		int line_end = 0;

		/* See if there's a full line in the buffer yet */
		for ( i=0; i<bd->rbufpos; i++ ) {

			/* Is there a line in the buffer? */
			if ( bd->rbuffer[i] == '\n' ) {
				bd->rbuffer[i] = '\0';
				line_end = i;
				line_ready = 1;
				break;
			}

		}

		if ( line_ready ) {

			int new_rbuflen;

			if ( bd->line != NULL ) {
				free(bd->line);
			}

			bd->line = strdup(bd->rbuffer);

			/* Now the block's been parsed, it should be
			 * forgotten about */
			memmove(bd->rbuffer,
			        bd->rbuffer + line_end + 1,
			        bd->rbuflen - line_end - 1);

			/* Subtract the number of bytes removed */
			bd->rbufpos = bd->rbufpos - line_end - 1;
			new_rbuflen = bd->rbuflen - line_end - 1 ;
			if ( new_rbuflen == 0 ) new_rbuflen = 256;
			bd->rbuffer = realloc(bd->rbuffer,
			                      new_rbuflen*sizeof(char));
			bd->rbuflen = new_rbuflen;

			return 1;

		} else {

			if ( bd->rbufpos == bd->rbuflen ) {
				bd->rbuffer = realloc(bd->rbuffer,
				                      bd->rbuflen + 256);
				bd->rbuflen = bd->rbuflen + 256;
			}
			no_line = 1;

		}

	}

	return 0;
}


static void run_work(const struct index_args *iargs,
                     int filename_pipe, int results_pipe, Stream *st,
                     int cookie, const char *tmpdir)
{
	FILE *fh;
	int allDone = 0;
	int w;
	unsigned int opts;
	struct buffer_data bd;

	bd.rbuffer = malloc(256*sizeof(char));
	bd.rbuflen = 256;
	bd.rbufpos = 0;
	bd.line = NULL;
	bd.fd = 0;

	fh = fdopen(filename_pipe, "r");
	if ( fh == NULL ) {
		ERROR("Failed to fdopen() the filename pipe!\n");
		return;
	}

	w = write(results_pipe, "\n", 1);
	if ( w < 0 ) {
		ERROR("Failed to send request for first filename.\n");
	}

	bd.fd = fileno(fh);

	/* Set non-blocking */
	opts = fcntl(bd.fd, F_GETFL);
	fcntl(bd.fd, F_SETFL, opts | O_NONBLOCK);

	while ( !allDone ) {

		struct pattern_args pargs;
		int  c;
		int error;
		int rval;
		char buf[1024];

		error = 0;
		pargs.filename_p_e = initialize_filename_plus_event();

		rval = 0;
		do {

			fd_set fds;
			struct timeval tv;
			int sval;

			FD_ZERO(&fds);
			FD_SET(bd.fd, &fds);

			tv.tv_sec = 30;
			tv.tv_usec = 0;

			sval = select(bd.fd+1, &fds, NULL, NULL, &tv);

			if ( sval == -1 ) {

				const int err = errno;

				switch ( err ) {

					case EINTR:
					STATUS("Restarting select()\n");
					break;

					default:
					ERROR("select() failed: %s\n",
					      strerror(err));
					rval = 1;

				}

			} else if ( sval != 0 ) {
				rval = read_fpe_data(&bd);
			} else {
				ERROR("No data sent from main process..\n");
				rval = 1;
				error = 1;
			}

		} while ( !rval );

		if ( error == 1 ) {
			allDone = 1;
			continue;
		}

		chomp(bd.line);

		if ( strlen(bd.line) == 0 ) {

			allDone = 1;

		} else {

			char filename[1024];
			char event_str[1024];
			struct event* ev;
			int ser;

			sscanf(bd.line, "%s %s %i", filename, event_str, &ser);
			pargs.filename_p_e->filename = strdup(filename);

			if ( strcmp(event_str, "/") != 0 ) {

				ev = get_event_from_event_string(event_str);
				if ( ev == NULL ) {
					ERROR("Error in event recovery\n");
				}

				pargs.filename_p_e->ev = ev;

			} else {

				pargs.filename_p_e->ev = NULL;

			}

			pargs.n_crystals = 0;
			process_image(iargs, &pargs, st, cookie, tmpdir,
			              results_pipe, ser);

			/* Request another image */
			c = sprintf(buf, "%i\n", pargs.n_crystals);
			w = write(results_pipe, buf, c);
			if ( w < 0 ) {
				ERROR("write P0\n");
			}

		}

		free_filename_plus_event(pargs.filename_p_e);

	}

	free(bd.line);
	free(bd.rbuffer);

	cleanup_indexing(iargs->indm, iargs->ipriv);
	free(iargs->indm);
	free(iargs->ipriv);
	free_detector_geometry(iargs->det);
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

size_t vol = 0;


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
		char *tmp;
		struct stat s;
		size_t ll;

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
		for ( j=0; j<sb->n_proc; j++ ) {
			if ( (j != slot) && (sb->running[j]) ) {
				close(sb->stream_pipe_write[j]);
			}
		}
		for ( j=0; j<sb->n_proc; j++ ) {
			if ( (j != slot) && (sb->running[j]) ) {
				if ( sb->result_fhs[j] != NULL ) {
					fclose(sb->result_fhs[j]);
				}
				close(sb->filename_pipes[j]);
			}
		}
		free(sb->filename_pipes);
		free(sb->result_fhs);
		free(sb->pids);
		/* Also prefix, tempdir, */

		/* Child process gets the 'read' end of the filename
		 * pipe, and the 'write' end of the result pipe. */
		close(filename_pipe[1]);
		close(result_pipe[0]);

		st = open_stream_fd_for_write(stream_pipe[1]);
		run_work(sb->iargs, filename_pipe[0], result_pipe[1],
		         st, slot, tmp);
		close_stream(st);

		//close(filename_pipe[0]);
		close(result_pipe[1]);

		free(sb);

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
				STATUS("Last filename was: %s (%s)\n",
				       sb->last_filename[i]->filename,
				       get_event_string(sb->last_filename[i]->ev) );
				sb->n_processed++;
				start_worker_process(sb, i);
			}

		}

	}
	unlock_sandbox(sb);
}


void create_sandbox(struct index_args *iargs, int n_proc, char *prefix,
                    int config_basename, FILE *fh,
                    Stream *stream, const char *tempdir)
{
	int i;
	int allDone;
	struct sigaction sa;
	int r;
	pthread_t reader_thread;
	struct sandbox *sb;
	size_t ll;
	struct stat s;

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
	sb->suspend_stats = 0;
	sb->n_proc = n_proc;
	sb->iargs = iargs;
	sb->serial = 1;

	sb->reader->fds = NULL;
	sb->reader->fhs = NULL;
	sb->reader->stream = stream;

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

	if ( tempdir == NULL ) {
		tempdir = strdup("");
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

	/* Fork the right number of times */
	lock_sandbox(sb);
	for ( i=0; i<n_proc; i++ ) {
		start_worker_process(sb, i);
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

		if ( FD_ISSET(signal_pipe[0], &fds) ) {

			char d;
			read(signal_pipe[0], &d, 1);
			handle_zombie(sb);

		}

		lock_sandbox(sb);
		for ( i=0; i<n_proc; i++ ) {

			struct filename_plus_event *nextImage;
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

			if ( strcmp(results, "SUSPEND") == 0 ) {
				sb->suspend_stats++;
			} else if ( strcmp(results, "RELEASE") == 0 ) {
				if ( sb->suspend_stats > 0 ) {
					sb->suspend_stats--;
				} else {
					ERROR("RELEASE before SUSPEND.\n");
				}
			} else {

				strtol(results, &eptr, 10);
				if ( eptr == results ) {
					if ( strlen(results) > 0 ) {
						ERROR("Invalid result '%s'\n",
						      results);
					}
				} else {

					int nc = atoi(results);
					sb->n_crystals += nc;
					if ( nc > 0 ) {
						sb->n_hadcrystals++;
					}
					sb->n_processed++;

				}

			}

			/* Send next filename */
			nextImage = get_pattern(fh, config_basename,
			                        iargs->det, prefix);

			if ( sb->last_filename[i] != NULL ) {
				free_filename_plus_event(sb->last_filename[i]);
			}

			sb->last_filename[i] = nextImage;

			if ( nextImage == NULL ) {

				/* No more images */
				r = write(sb->filename_pipes[i], "\n", 1);
				if ( r < 0 ) {
					ERROR("Write pipe\n");
				}

			} else {

				char tmp[256];

				r = write(sb->filename_pipes[i],
				          nextImage->filename,
				          strlen(nextImage->filename));

				if ( r < 0 ) {
					ERROR("write pipe\n");
				}

				r = write(sb->filename_pipes[i], " ", 1);
				if ( r < 0 ) {
					ERROR("write pipe\n");
				}

				if ( nextImage->ev != NULL ) {

					r = write(sb->filename_pipes[i],
					          get_event_string(nextImage->ev),
					          strlen(get_event_string(nextImage->ev)));
					if ( r < 0 ) {
						ERROR("write pipe\n");
					}

				} else {

					r = write(sb->filename_pipes[i], "/", 1);
					if ( r < 0 ) {
						ERROR("write pipe\n");
					}

				}

				snprintf(tmp, 255, " %i", sb->serial++);
				r = write(sb->filename_pipes[i],
				          tmp, strlen(tmp));
				if ( r < 0 ) {
					ERROR("write pipe\n");
				}

				r = write(sb->filename_pipes[i], "\n", 1);
				if ( r < 0 ) {
					ERROR("write pipe\n");
				}

			}
		}

		unlock_sandbox(sb);

		/* Update progress */
		lock_sandbox(sb);
		tNow = get_monotonic_seconds();
		if ( !sb->suspend_stats
		    && (tNow >= sb->t_last_stats+STATS_EVERY_N_SECONDS) )
		{

			STATUS("%4i indexable out of %4i processed (%4.1f%%), "
			       "%4i crystals so far. "
			       "%4i images processed since the last message.\n",
			       sb->n_hadcrystals, sb->n_processed,
			       (sb->n_processed == 0 ? 0 :
			       100.0 * sb->n_hadcrystals / sb->n_processed),
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
	free(sb->tmpdir);

	pthread_mutex_destroy(&sb->lock);

	STATUS("Final:"
	       " %i images processed, %i had crystals (%.1f%%),"
	       " %i crystals overall.\n",
	       sb->n_processed, sb->n_hadcrystals,
	       100.0 * sb->n_hadcrystals / sb->n_processed, sb->n_crystals);

	free(sb);
}
