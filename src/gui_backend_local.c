/*
 * gui_backend_local.c
 *
 * GUI backend for running jobs on the local machine
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020 Thomas White <taw@physics.org>
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

#include <pty.h>
#include <glib.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gtk/gtk.h>

#include <utils.h>

#include "gui_project.h"
#include "gui_index.h"


struct local_indexing_opts
{
	int n_processes;
};


struct local_job
{
	double frac_complete;
	int n_frames;

	/* When both these are true, free the job resources */
	int indexamajig_running;
	int cancelled;

	guint indexamajig_watch;
	GPid indexamajig_pid;
	guint child_watch_source;
	guint index_readable_source;
};


static void watch_indexamajig(GPid pid, gint status, gpointer vp)
{
	struct local_job *job = vp;
	STATUS("Indexamajig exited with status %i\n", status);
	job->indexamajig_running = 0;
	g_spawn_close_pid(job->indexamajig_pid);
}


static gboolean index_readable(GIOChannel *source, GIOCondition cond,
                               void *vp)
{
	GIOStatus r;
	GError *err = NULL;
	struct local_job *job = vp;
	gchar *line;

	r = g_io_channel_read_line(source, &line, NULL, NULL, &err);
	if ( r == G_IO_STATUS_EOF ) {
		STATUS("End of output.\n");
		return FALSE;
	}
	if ( r != G_IO_STATUS_NORMAL ) {
		if ( job->indexamajig_pid != 0 ) {
			STATUS("Read error?\n");
		} else {
			STATUS("End of output (indexamajig exited)\n");
		}
		return FALSE;
	}
	chomp(line);

	if ( strncmp(line, "Final: ", 7) == 0 ) {
		job->frac_complete = 1.0;
	} else if ( strstr(line, " images processed, ") != NULL ) {
		int n_proc;
		sscanf(line, "%i ", &n_proc);
		job->frac_complete = (double)n_proc/job->n_frames;
	} else {
		STATUS("%s\n", line);
	}

	g_free(line);

	return TRUE;
}


static int write_file_list(char **filenames,
                           char **events,
                           int n_frames)
{
	FILE *fh;
	int i;

	fh = fopen("files.lst", "w");
	if ( fh == NULL ) return 1;

	for ( i=0; i<n_frames; i++ ) {
		fprintf(fh, "%s", filenames[i]);
		if ( events[i] != NULL ) {
			fprintf(fh, " %s\n", events[i]);
		} else {
			fprintf(fh, "\n");
		}
	}

	fclose(fh);

	return 0;
}


void setup_subprocess(gpointer user_data)
{
	const char *workdir = user_data;
	setsid();
	setpgid(0, 0);
	chdir(workdir);
}


static void *run_indexing(const char *job_title,
                          const char *job_notes,
                          char **filenames,
                          char **events,
                          int n_frames,
                          char *geom_filename,
                          struct peak_params *peak_search_params,
                          struct index_params *indexing_params,
                          void *opts_priv)
{
	struct local_indexing_opts *opts = opts_priv;
	GIOChannel *ioch;
	char n_thread_str[64];
	char **args;
	int i;
	int r;
	int ch_stderr;
	GError *error;
	struct local_job *job;
	struct stat s;
	char *workdir;
	const char *old_pwd;

	workdir = strdup(job_title);
	if ( workdir == NULL ) return NULL;

	if ( stat(workdir, &s) != -1 ) {
		ERROR("Working directory already exists.  "
		      "Choose a different job name.\n");
		return NULL;
	}

	if ( mkdir(workdir, S_IRWXU) ) {
		ERROR("Failed to create working directory: %s\n",
		      strerror(errno));
		return NULL;
	}

	job = malloc(sizeof(struct local_job));
	if ( job == NULL ) return NULL;

	old_pwd = getcwd(NULL, 0);
	chdir(workdir);
	if ( write_file_list(filenames, events, n_frames) ) {
		STATUS("Failed to write list\n");
		free(job);
		return NULL;
	}
	chdir(old_pwd);

	job->n_frames = n_frames;
	job->frac_complete = 0.0;

	snprintf(n_thread_str, 63, "%i", opts->n_processes);
	args = indexamajig_command_line(geom_filename,
	                                n_thread_str,
	                                peak_search_params,
	                                indexing_params);

	i = 0;
	while ( args[i] != NULL ) {
		STATUS("%s ", args[i++]);
	}
	STATUS("\n");

	r = g_spawn_async_with_pipes(NULL, args, NULL,
	                             G_SPAWN_SEARCH_PATH
	                           | G_SPAWN_DO_NOT_REAP_CHILD,
	                             setup_subprocess, workdir,
	                             &job->indexamajig_pid,
	                             NULL, NULL, &ch_stderr,
	                             &error);
	if ( r == FALSE ) {
		ERROR("Failed to run indexamajig: %s\n",
		      error->message);
		free(job);
		return NULL;
	}
	job->indexamajig_running = 1;

	job->child_watch_source = g_child_watch_add(job->indexamajig_pid,
	                                            watch_indexamajig,
	                                            job);

	ioch = g_io_channel_unix_new(ch_stderr);
	job->index_readable_source = g_io_add_watch(ioch,
	                                            G_IO_IN | G_IO_ERR | G_IO_HUP,
	                                            index_readable,
	                                            job);

	return job;
}


static int get_task_status(void *job_priv,
                           int *running,
                           float *frac_complete)
{
	struct local_job *job = job_priv;
	*frac_complete = job->frac_complete;
	*running = job->indexamajig_running;
	return 0;
}


static void cancel_task(void *job_priv)
{
	struct local_job *job = job_priv;

	if ( !job->indexamajig_running ) return;

	ERROR("Stopping indexamajig (pid %i).\n", job->indexamajig_pid);
	kill(-job->indexamajig_pid, SIGINT);
}


static void n_processes_activate_sig(GtkEntry *entry, gpointer data)
{
	struct local_indexing_opts *opts = data;
	convert_int(gtk_entry_get_text(entry), &opts->n_processes);
}


static gboolean n_processes_focus_sig(GtkEntry *entry, GdkEvent *event,
                                      gpointer data)
{
	n_processes_activate_sig(entry, data);
	return FALSE;
}


static GtkWidget *make_indexing_parameters_widget(void *opts_priv)
{
	struct local_indexing_opts *opts = opts_priv;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *entry;
	GtkWidget *label;
	char tmp[64];

	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Number of threads:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	snprintf(tmp, 63, "%i", opts->n_processes);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);

	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(n_processes_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(n_processes_focus_sig),
	                 opts);
	return vbox;
}


static struct local_indexing_opts *make_default_local_opts()
{
	struct local_indexing_opts *opts = malloc(sizeof(struct local_indexing_opts));
	if ( opts == NULL ) return NULL;

	opts->n_processes = 4;

	return opts;
}


static void write_indexing_opts(void *opts_priv, FILE *fh)
{
	struct local_indexing_opts *opts = opts_priv;

	fprintf(fh, "indexing.local.n_processes %i\n",
	        opts->n_processes);
}


static void read_indexing_opt(void *opts_priv,
                              const char *key,
                              const char *val)
{
	struct local_indexing_opts *opts = opts_priv;

	if ( strcmp(key, "indexing.local.n_processes") == 0 ) {
		if ( convert_int(val, &opts->n_processes) ) {
			ERROR("Invalid number of threads: %s\n", val);
		}
	}
}


int make_local_backend(struct crystfel_backend *be)
{
	be->name = "local";
	be->friendly_name = "Local (run on this computer)";

	be->make_indexing_parameters_widget = make_indexing_parameters_widget;
	be->run_indexing = run_indexing;
	be->cancel_task = cancel_task;
	be->task_status = get_task_status;
	be->indexing_opts_priv = make_default_local_opts();
	if ( be->indexing_opts_priv == NULL ) return 1;
	be->write_indexing_opts = write_indexing_opts;
	be->read_indexing_opt = read_indexing_opt;

	return 0;
};
