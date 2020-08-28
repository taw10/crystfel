/*
 * gui_backend_slurm.c
 *
 * GUI backend for running jobs via SLURM
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

#include <glib.h>
#include <gtk/gtk.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <slurm/slurm.h>
#include <gio/gio.h>

#include <utils.h>

#include "gui_project.h"
#include "gui_index.h"


struct slurm_indexing_opts
{
	char *partition;
	int block_size;
	char *email_address;
};


struct slurm_job
{
	double frac_complete;
	/* FIXME: List of SLURM job numbers to track */
};


static int get_task_status(void *job_priv,
                           int *running,
                           float *frac_complete)
{
	struct slurm_job *job = job_priv;
	*frac_complete = job->frac_complete;
	*running = 1;
	return 0;
}


static void cancel_task(void *job_priv)
{
	//struct slurm_job *job = job_priv;
}


static char **create_env(uint32_t *psize)
{
	char **env;

	env = malloc(10*sizeof(char *));
	if ( env == NULL ) return NULL;

	env[0] = strdup("PATH=/path/to/indexamajig");
	*psize = 1;

	return env;
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
	struct slurm_indexing_opts *opts = opts_priv;
	struct slurm_job *job;
	job_desc_msg_t job_desc_msg;
	submit_response_msg_t *resp;
	int r;
	char *workdir;
	struct stat s;
	char **cmdline;
	char *cmdline_all;
	char *script;
	GFile *workdir_file;
	GFile *cwd_file;
	GFile *notes_file;
	char *notes_path;
	FILE *fh;

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

	cwd_file = g_file_new_for_path(".");
	workdir_file = g_file_get_child(cwd_file, workdir);
	g_object_unref(cwd_file);

	notes_file = g_file_get_child(workdir_file, "notes.txt");
	notes_path = g_file_get_path(notes_file);
	fh = fopen(notes_path, "w");
	fputs(job_notes, fh);
	fclose(fh);
	g_free(notes_path);
	g_object_unref(notes_file);

	cmdline = indexamajig_command_line(geom_filename,
	                                   "`nproc`",
	                                   peak_search_params,
	                                   indexing_params);

	cmdline_all = g_strjoinv(" ", cmdline);

	script = malloc(strlen(cmdline_all)+16);
	if ( script == NULL ) return NULL;

	strcpy(script, "#!/bin/sh\n");
	strcat(script, cmdline_all);
	g_free(cmdline_all);

	job = malloc(sizeof(struct slurm_job));
	if ( job == NULL ) return NULL;

	slurm_init_job_desc_msg(&job_desc_msg);
	job_desc_msg.user_id = getuid();
	job_desc_msg.group_id = getgid();
	job_desc_msg.mail_user = strdup(opts->email_address);
	job_desc_msg.mail_type = MAIL_JOB_FAIL;
	job_desc_msg.comment = strdup("Submitted via CrystFEL GUI");
	job_desc_msg.shared = 0;
	job_desc_msg.time_limit = 60;
	job_desc_msg.partition = strdup(opts->partition);
	job_desc_msg.min_nodes = 1;
	job_desc_msg.max_nodes = 1;
	job_desc_msg.name = strdup(job_title);
	job_desc_msg.std_err = strdup("job.err");
	job_desc_msg.std_out = strdup("job.out");
	job_desc_msg.work_dir = g_file_get_path(workdir_file);
	job_desc_msg.script = script;
	job_desc_msg.environment = create_env(&job_desc_msg.env_size);

	g_object_unref(workdir_file);

	r = slurm_submit_batch_job(&job_desc_msg, &resp);

	if ( r ) {
		ERROR("Couldn't submit job: %i\n", resp->error_code);
		free(job);
		return NULL;
	}

	STATUS("Submitted SLURM job ID %i\n", resp->job_id);
	slurm_free_submit_response_response_msg(resp);

	return job;
}


static void block_size_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	convert_int(gtk_entry_get_text(entry), &opts->block_size);
}


static gboolean block_size_focus_sig(GtkEntry *entry, GdkEvent *event,
                                     gpointer data)
{
	block_size_activate_sig(entry, data);
	return FALSE;
}


static void partition_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	opts->partition = strdup(gtk_entry_get_text(entry));
}


static gboolean partition_focus_sig(GtkEntry *entry, GdkEvent *event,
                                    gpointer data)
{
	partition_activate_sig(entry, data);
	return FALSE;
}


static void email_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	opts->email_address = strdup(gtk_entry_get_text(entry));
}


static gboolean email_focus_sig(GtkEntry *entry, GdkEvent *event,
                                gpointer data)
{
	email_activate_sig(entry, data);
	return FALSE;
}


static GtkWidget *make_indexing_parameters_widget(void *opts_priv)
{
	struct slurm_indexing_opts *opts = opts_priv;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *entry;
	GtkWidget *label;
	char tmp[64];

	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Submit job to partition:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->partition != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->partition);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry), "maxwell");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(partition_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(partition_focus_sig),
	                 opts);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Split job into blocks of");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	snprintf(tmp, 63, "%i", opts->block_size);
	entry = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(block_size_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(block_size_focus_sig),
	                 opts);
	label = gtk_label_new("frames");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Send notifications to:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->email_address != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->email_address);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "myself@example.org");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(email_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(email_focus_sig),
	                 opts);

	return vbox;
}


static struct slurm_indexing_opts *make_default_slurm_opts()
{
	struct slurm_indexing_opts *opts = malloc(sizeof(struct slurm_indexing_opts));
	if ( opts == NULL ) return NULL;

	opts->partition = NULL;
	opts->block_size = 1000;
	opts->email_address = NULL;

	return opts;
}


static void write_indexing_opts(void *opts_priv, FILE *fh)
{
	struct slurm_indexing_opts *opts = opts_priv;

	fprintf(fh, "indexing.slurm.block_size %i\n",
	        opts->block_size);

	if ( opts->partition != NULL) {
		fprintf(fh, "indexing.slurm.partition %s\n",
		        opts->partition);
	}

	if ( opts->email_address != NULL ) {
		fprintf(fh, "indexing.slurm.email_address %s\n",
		        opts->email_address);
	}
}


static void read_indexing_opt(void *opts_priv,
                              const char *key,
                              const char *val)
{
	struct slurm_indexing_opts *opts = opts_priv;

	if ( strcmp(key, "indexing.slurm.block_size") == 0 ) {
		if ( convert_int(val, &opts->block_size) ) {
			ERROR("Invalid block size: %s\n", val);
		}
	}

	if ( strcmp(key, "indexing.slurm.email_address") == 0 ) {
		opts->email_address = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.partition") == 0 ) {
		opts->partition = strdup(val);
	}
}


int make_slurm_backend(struct crystfel_backend *be)
{
	be->name = "slurm";
	be->friendly_name = "SLURM";

	be->make_indexing_parameters_widget = make_indexing_parameters_widget;
	be->run_indexing = run_indexing;
	be->write_indexing_opts = write_indexing_opts;
	be->read_indexing_opt = read_indexing_opt;
	be->cancel_task = cancel_task;
	be->task_status = get_task_status;

	be->indexing_opts_priv = make_default_slurm_opts();
	if ( be->indexing_opts_priv == NULL ) return 1;

	return 0;
};
