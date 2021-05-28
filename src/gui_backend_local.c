/*
 * gui_backend_local.c
 *
 * GUI backend for running jobs on the local machine
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2020-2021 Thomas White <taw@physics.org>
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

#include <glib.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <gtk/gtk.h>

#include <utils.h>

#include "crystfel_gui.h"
#include "gui_project.h"
#include "gui_index.h"
#include "gui_merge.h"
#include "gui_ambi.h"


struct local_indexing_opts
{
	int n_processes;
};


struct local_merging_opts
{
	int n_threads;
};


struct local_ambi_opts
{
	int n_threads;
};


struct local_job
{
	enum gui_job_type type;

	int n_frames;
	int niter;

	/* When both these are true, free the job resources */
	int running;
	int cancelled;

	char *stderr_filename;

	GPid pid;
	guint child_watch_source;
	GFile *workdir;
};


static void free_task(void *job_priv)
{
	struct local_job *job = job_priv;
	g_object_unref(job->workdir);
	free(job->stderr_filename);
}


static void watch_subprocess(GPid pid, gint status, gpointer vp)
{
	struct local_job *job = vp;
	STATUS("Subprocess exited with status %i\n", status);
	job->running = 0;
	g_spawn_close_pid(job->pid);
}


static int write_file_list(GFile *workdir,
                           const char *listname,
                           char **filenames,
                           char **events,
                           int n_frames)
{
	FILE *fh;
	int i;
	GFile *list_gfile;
	char *list_str;

	list_gfile = g_file_get_child(workdir, listname);
	list_str = g_file_get_path(list_gfile);
	if ( list_str == NULL ) return 1;

	fh = fopen(list_str, "w");
	free(list_str);
	if ( fh == NULL ) return 1;

	for ( i=0; i<n_frames; i++ ) {
		if ( filenames[i][0] != '/' ) {
			fprintf(fh, "../%s", filenames[i]);
		} else {
			fprintf(fh, "%s", filenames[i]);
		}
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


static struct local_job *start_local_job(char **args,
                                         const char *job_title,
                                         GFile *workdir_file,
                                         struct crystfelproject *proj,
                                         enum gui_job_type type)
{
	int i;
	int r;
	int ch_stderr;
	GError *error;
	struct local_job *job;
	char *workdir_str;
	GFile *stderr_gfile;

	workdir_str = g_file_get_path(workdir_file);
	if ( workdir_str == NULL ) return NULL;

	job = malloc(sizeof(struct local_job));
	if ( job == NULL ) return NULL;

	job->workdir = g_file_dup(workdir_file);
	job->type = type;

	STATUS("Running program: ");
	i = 0;
	while ( args[i] != NULL ) {
		STATUS("%s ", args[i++]);
	}
	STATUS("\n");

	error = NULL;
	r = g_spawn_async_with_pipes(NULL, args, NULL,
	                             G_SPAWN_SEARCH_PATH
	                           | G_SPAWN_DO_NOT_REAP_CHILD,
	                             setup_subprocess, workdir_str,
	                             &job->pid,
	                             NULL, NULL, &ch_stderr,
	                             &error);
	if ( r == FALSE ) {
		ERROR("Failed to start program: %s\n", error->message);
		free(job);
		return NULL;
	}
	job->running = 1;

	stderr_gfile = g_file_get_child(workdir_file, "stderr.log");
	job->stderr_filename = g_file_get_path(stderr_gfile);
	g_object_unref(stderr_gfile);

	job->child_watch_source = g_child_watch_add(job->pid,
	                                            watch_subprocess,
	                                            job);

	return job;
}


static int get_task_status(void *job_priv,
                           int *running,
                           float *frac_complete)
{
	int n_proc;
	struct local_job *job = job_priv;

	*running = job->running;

	switch ( job->type ) {

		case GUI_JOB_INDEXING :
		n_proc = read_number_processed(job->stderr_filename);
		*frac_complete = (double)n_proc / job->n_frames;
		break;

		case GUI_JOB_AMBIGATOR :
		*frac_complete = read_ambigator_progress(job->stderr_filename,
		                                         job->niter);
		break;

		case GUI_JOB_PROCESS_HKL :
		case GUI_JOB_PROCESS_HKL_SCALE :
		case GUI_JOB_PARTIALATOR :
		*frac_complete = read_merge_progress(job->stderr_filename,
		                                     job->type);
		break;

	}

	return 0;
}


static void cancel_task(void *job_priv)
{
	struct local_job *job = job_priv;

	if ( !job->running ) return;

	ERROR("Stopping indexamajig (pid %i).\n", job->pid);
	kill(-job->pid, SIGINT);
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


static struct local_indexing_opts *make_default_local_indexing_opts()
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


static void n_threads_activate_sig(GtkEntry *entry, gpointer data)
{
	struct local_merging_opts *opts = data;
	convert_int(gtk_entry_get_text(entry), &opts->n_threads);
}


static gboolean n_threads_focus_sig(GtkEntry *entry, GdkEvent *event,
                                    gpointer data)
{
	n_threads_activate_sig(entry, data);
	return FALSE;
}


static GtkWidget *make_merging_parameters_widget(void *opts_priv)
{
	struct local_merging_opts *opts = opts_priv;
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
	snprintf(tmp, 63, "%i", opts->n_threads);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);

	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(n_threads_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(n_threads_focus_sig),
	                 opts);
	return vbox;
}


static void *run_ambi(const char *job_title,
                      const char *job_notes,
                      struct crystfelproject *proj,
                      struct gui_indexing_result *input,
                      void *opts_priv)
{
	char n_thread_str[64];
	struct local_job *job;
	struct local_merging_opts *opts = opts_priv;
	GFile *workdir;
	GFile *sc_gfile;
	gchar *sc_filename;
	GFile *stream_gfile;
	char *stream_str;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	stream_gfile = g_file_get_child(workdir, "ambi.stream");
	stream_str = g_file_get_path(stream_gfile);
	g_object_unref(stream_gfile);

	snprintf(n_thread_str, 64, "%i", opts->n_threads);
	sc_gfile = g_file_get_child(workdir, "run_ambigator.sh");
	sc_filename = g_file_get_path(sc_gfile);
	g_object_unref(sc_gfile);
	if ( sc_filename == NULL ) return NULL;
	if ( !write_ambigator_script(sc_filename, input, n_thread_str,
	                             &proj->ambi_params, stream_str) )
	{
		char *args[3];
		args[0] = "sh";
		args[1] = sc_filename;
		args[2] = NULL;
		job = start_local_job(args, job_title, workdir,
		                      proj, GUI_JOB_AMBIGATOR);
		job->niter = proj->ambi_params.niter;
	} else {
		job = NULL;
	}

	if ( job != NULL ) {
		add_indexing_result(proj, job_title, &stream_str, 1);
	}

	g_object_unref(workdir);
	return job;
}


static void *run_merging(const char *job_title,
                         const char *job_notes,
                         struct crystfelproject *proj,
                         struct gui_indexing_result *input,
                         void *opts_priv)
{
	char n_thread_str[64];
	struct local_job *job;
	struct local_merging_opts *opts = opts_priv;
	GFile *workdir;
	GFile *sc_gfile;
	gchar *sc_filename;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	snprintf(n_thread_str, 63, "%i", opts->n_threads);
	sc_gfile = g_file_get_child(workdir, "run_merge.sh");
	sc_filename = g_file_get_path(sc_gfile);
	g_object_unref(sc_gfile);
	if ( sc_filename == NULL ) return NULL;

	if ( !write_merge_script(sc_filename, input, n_thread_str,
	                         &proj->merging_params, "crystfel.hkl") )
	{
		char *args[3];
		enum gui_job_type type;
		args[0] = "sh";
		args[1] = sc_filename;
		args[2] = NULL;
		if ( strcmp(proj->merging_params.model, "process_hkl") == 0 ) {
			if ( proj->merging_params.scale ) {
				type = GUI_JOB_PROCESS_HKL_SCALE;
			} else {
				type = GUI_JOB_PROCESS_HKL;
			}
		} else {
			type = GUI_JOB_PARTIALATOR;
		}
		job = start_local_job(args, job_title, workdir, proj, type);
	} else {
		job = NULL;
	}
	g_free(sc_filename);

	if ( job != NULL ) {

		GFile *hkl_gfile;
		char *hkl;
		char *hkl1;
		char *hkl2;

		hkl_gfile = g_file_get_child(workdir, "crystfel.hkl");
		hkl = g_file_get_path(hkl_gfile);
		g_object_unref(hkl_gfile);

		hkl_gfile = g_file_get_child(workdir, "crystfel.hkl1");
		hkl1 = g_file_get_path(hkl_gfile);
		g_object_unref(hkl_gfile);

		hkl_gfile = g_file_get_child(workdir, "crystfel.hkl2");
		hkl2 = g_file_get_path(hkl_gfile);
		g_object_unref(hkl_gfile);

		add_merge_result(proj, job_title, hkl, hkl1, hkl2);
		g_free(hkl);
		g_free(hkl1);
		g_free(hkl2);
	}

	g_object_unref(workdir);
	return job;
}


static void *run_indexing(const char *job_title,
                          const char *job_notes,
                          struct crystfelproject *proj,
                          void *opts_priv)
{
	struct local_indexing_opts *opts = opts_priv;
	struct local_job *job;
	char n_thread_str[64];
	GFile *workdir;
	GFile *sc_gfile;
	gchar *sc_filename;
	GFile *stderr_gfile;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	if ( write_file_list(workdir, "files.lst",
		             proj->filenames,
		             proj->events,
		             proj->n_frames) )
	{
		STATUS("Failed to write list\n");
		return NULL;
	}

	snprintf(n_thread_str, 63, "%i", opts->n_processes);
	sc_gfile = g_file_get_child(workdir, "run_indexamajig.sh");
	sc_filename = g_file_get_path(sc_gfile);
	g_object_unref(sc_gfile);
	if ( sc_filename == NULL ) return NULL;
	if ( !write_indexamajig_script(sc_filename,
	                               proj->geom_filename,
	                               n_thread_str,
	                               "files.lst",
	                               "crystfel.stream",
	                               NULL, 1,
	                               &proj->peak_search_params,
	                               &proj->indexing_params) )
	{
		char *args[3];
		args[0] = "sh";
		args[1] = sc_filename;
		args[2] = NULL;
		job = start_local_job(args, job_title, workdir,
		                      proj, GUI_JOB_INDEXING);
	} else {
		job = NULL;
	}

	if ( job != NULL ) {

		char *stream_fn;

		/* Indexing-specific job data */
		job->n_frames = proj->n_frames;

		stderr_gfile = g_file_get_child(workdir, "stderr.log");
		job->stderr_filename = g_file_get_path(stderr_gfile);
		g_object_unref(stderr_gfile);

		GFile *stream_gfile = g_file_get_child(job->workdir,
		                                       "crystfel.stream");
		stream_fn = g_file_get_path(stream_gfile);
		g_object_unref(stream_gfile);
		add_indexing_result(proj, job_title, &stream_fn, 1);
		g_free(stream_fn);
	}
	g_object_unref(workdir);

	return job;
}


static struct local_merging_opts *make_default_local_merging_opts()
{
	struct local_merging_opts *opts = malloc(sizeof(struct local_merging_opts));
	if ( opts == NULL ) return NULL;

	opts->n_threads = 4;

	return opts;
}


static void write_merging_opts(void *opts_priv, FILE *fh)
{
	struct local_merging_opts *opts = opts_priv;

	fprintf(fh, "merging.local.n_threads %i\n",
	        opts->n_threads);
}


static void read_merging_opt(void *opts_priv,
                             const char *key,
                             const char *val)
{
	struct local_merging_opts *opts = opts_priv;

	if ( strcmp(key, "merging.local.n_threads") == 0 ) {
		if ( convert_int(val, &opts->n_threads) ) {
			ERROR("Invalid number of threads: %s\n", val);
		}
	}
}


static GtkWidget *make_ambi_parameters_widget(void *opts_priv)
{
	struct local_ambi_opts *opts = opts_priv;
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
	snprintf(tmp, 63, "%i", opts->n_threads);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);

	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(n_threads_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(n_threads_focus_sig),
	                 opts);
	return vbox;
}


static struct local_ambi_opts *make_default_local_ambi_opts()
{
	struct local_ambi_opts *opts = malloc(sizeof(struct local_ambi_opts));
	if ( opts == NULL ) return NULL;

	opts->n_threads = 4;

	return opts;
}


static void write_ambi_opts(void *opts_priv, FILE *fh)
{
	struct local_ambi_opts *opts = opts_priv;

	fprintf(fh, "ambi.local.n_threads %i\n",
	        opts->n_threads);
}


static void read_ambi_opt(void *opts_priv,
                          const char *key,
                          const char *val)
{
	struct local_ambi_opts *opts = opts_priv;

	if ( strcmp(key, "ambi.local.n_threads") == 0 ) {
		if ( convert_int(val, &opts->n_threads) ) {
			ERROR("Invalid number of threads: %s\n", val);
		}
	}
}


int make_local_backend(struct crystfel_backend *be)
{
	be->name = "local";
	be->friendly_name = "Local (run on this computer)";

	be->cancel_task = cancel_task;
	be->free_task = free_task;
	be->task_status = get_task_status;

	be->make_indexing_parameters_widget = make_indexing_parameters_widget;
	be->run_indexing = run_indexing;
	be->indexing_opts_priv = make_default_local_indexing_opts();
	if ( be->indexing_opts_priv == NULL ) return 1;
	be->write_indexing_opts = write_indexing_opts;
	be->read_indexing_opt = read_indexing_opt;

	be->make_merging_parameters_widget = make_merging_parameters_widget;
	be->run_merging = run_merging;
	be->merging_opts_priv = make_default_local_merging_opts();
	if ( be->merging_opts_priv == NULL ) return 1;
	be->write_merging_opts = write_merging_opts;
	be->read_merging_opt = read_merging_opt;

	be->make_ambi_parameters_widget = make_ambi_parameters_widget;
	be->run_ambi = run_ambi;
	be->ambi_opts_priv = make_default_local_ambi_opts();
	if ( be->ambi_opts_priv == NULL ) return 1;
	be->write_ambi_opts = write_ambi_opts;
	be->read_ambi_opt = read_ambi_opt;

	return 0;
};
