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
	char *path_add;
};


struct slurm_job
{
	double frac_complete;
	int n_frames;
	int n_blocks;
	uint32_t *job_ids;
	char **stderr_filenames;
};


static int read_number_processed(const char *filename)
{
	FILE *fh = fopen(filename, "r");
	int n_proc;

	/* Normal situation if SLURM job hasn't started yet */
	if ( fh == NULL ) return 0;

	do {
		char line[1024];
		if ( fgets(line, 1024, fh) == NULL ) break;

		if ( strncmp(line, "Final: ", 7) == 0 ) {
			sscanf(line, "Final: %i images processed", &n_proc);
		} else if ( strstr(line, " images processed, ") != NULL ) {
			sscanf(line, "%i ", &n_proc);
		}

	} while ( 1 );

	fclose(fh);

	return n_proc;
}


static int job_running(uint32_t job_id)
{
	job_info_msg_t *job_info;
	int running = 1;

	if ( slurm_load_job(&job_info, job_id, 0) ) {
		STATUS("Couldn't get status: %i\n",
		       slurm_strerror(slurm_get_errno()));
		running = 0;
		/* FIXME: Distinguish error cond from job complete */
	}

	switch ( job_info->job_array[0].job_state & JOB_STATE_BASE ) {

		/* Only the following states are reasons to keep on watching
		 * the job */
		case JOB_PENDING :
		case JOB_RUNNING :
		case JOB_SUSPENDED :
		running = 1;
		break;

		default :
		running = 0;
		break;
	}

	slurm_free_job_info_msg(job_info);

	return running;
}


static int get_task_status(void *job_priv,
                           int *running,
                           float *frac_complete)
{
	struct slurm_job *job = job_priv;
	int i;
	int n_proc = 0;
	int all_complete = 1;

	for ( i=0; i<job->n_blocks; i++ ) {

		n_proc += read_number_processed(job->stderr_filenames[i]);

		if ( job->job_ids[i] == 0 ) continue;

		if ( !job_running(job->job_ids[i]) ) {
			job->job_ids[i] = 0;
		} else {
			all_complete = 0;
		}
	}

	*frac_complete = (double)n_proc / job->n_frames;
	*running = 1 - all_complete;
	return 0;
}


static void cancel_task(void *job_priv)
{
	int i;
	struct slurm_job *job = job_priv;
	for ( i=0; i<job->n_blocks; i++) {
		if ( job->job_ids[i] == 0 ) continue;
		STATUS("Stopping SLURM job %i\n", job->job_ids[i]);
		if ( slurm_kill_job(job->job_ids[i], SIGINT, 0) ) {
			ERROR("Couldn't stop job: %s\n",
			      slurm_strerror(slurm_get_errno()));
		}
	}
}


static char **create_env(int *psize, char *path_add)
{
	char **env;
	const char *base_path = "PATH=/bin:/usr/bin";
	char *crystfel_path;
	size_t path_len;

	env = malloc(10*sizeof(char *));
	if ( env == NULL ) return NULL;

	crystfel_path = get_crystfel_path_str();

	path_len = 4 + strlen(base_path);

	if ( path_add != NULL ) {
		path_len += strlen(path_add);
	}

	if ( crystfel_path != NULL ) {
		path_len += strlen(crystfel_path);
	}

	env[0] = malloc(path_len);
	if ( env[0] == NULL ) return NULL;

	strcpy(env[0], base_path);
	if ( crystfel_path != NULL ) {
		strcat(env[0], ":");
		strcat(env[0], crystfel_path);
		g_free(crystfel_path);
	}
	if ( path_add != NULL ) {
		strcat(env[0], ":");
		strcat(env[0], path_add);
	}

	*psize = 1;

	return env;
}


static uint32_t submit_batch_job(const char *geom_filename,
                                 const char *file_list,
                                 const char *stream_filename,
                                 const char *email_address,
                                 const char *partition,
                                 char **env,
                                 int n_env,
                                 const char *job_name,
                                 const char *workdir,
                                 const char *stderr_file,
                                 const char *stdout_file,
                                 struct peak_params *peak_search_params,
                                 struct index_params *indexing_params)

{
	job_desc_msg_t job_desc_msg;
	submit_response_msg_t *resp;
	char **cmdline;
	char *cmdline_all;
	char *script;
	int job_id;
	int r;

	cmdline = indexamajig_command_line(geom_filename,
	                                   "`nproc`",
	                                   file_list,
	                                   stream_filename,
	                                   peak_search_params,
	                                   indexing_params);

	cmdline_all = g_strjoinv(" ", cmdline);

	script = malloc(strlen(cmdline_all)+16);
	if ( script == NULL ) return 0;

	strcpy(script, "#!/bin/sh\n");
	strcat(script, cmdline_all);
	g_free(cmdline_all);

	slurm_init_job_desc_msg(&job_desc_msg);
	job_desc_msg.user_id = getuid();
	job_desc_msg.group_id = getgid();
	job_desc_msg.mail_user = safe_strdup(email_address);
	job_desc_msg.mail_type = MAIL_JOB_FAIL;
	job_desc_msg.comment = "Submitted via CrystFEL GUI";
	job_desc_msg.shared = 0;
	job_desc_msg.time_limit = 60;
	job_desc_msg.partition = safe_strdup(partition);
	job_desc_msg.min_nodes = 1;
	job_desc_msg.max_nodes = 1;
	job_desc_msg.name = safe_strdup(job_name);
	job_desc_msg.std_err = strdup(stderr_file);
	job_desc_msg.std_out = strdup(stdout_file);
	job_desc_msg.work_dir = strdup(workdir);
	job_desc_msg.script = script;
	job_desc_msg.environment = env;
	job_desc_msg.env_size = n_env;

	r = slurm_submit_batch_job(&job_desc_msg, &resp);
	if ( r ) {
		ERROR("Couldn't submit job: %i\n", errno);
		return 0;
	}

	free(job_desc_msg.mail_user);
	free(job_desc_msg.partition);
	free(job_desc_msg.name);
	free(job_desc_msg.work_dir);
	free(job_desc_msg.std_err);
	free(job_desc_msg.std_out);

	job_id = resp->job_id;
	slurm_free_submit_response_response_msg(resp);

	return job_id;
}


static void write_partial_file_list(GFile *workdir,
                                    const char *list_filename,
                                    int j,
                                    int block_size,
                                    char **filenames,
                                    char **events,
                                    int n_frames)
{
	GFile *file;
	char *file_path;
	FILE *fh;
	int i;

	file = g_file_get_child(workdir, list_filename);
	file_path = g_file_get_path(file);

	fh = fopen(file_path, "w");
	for ( i=j*block_size;
	      (i<(j+1)*block_size) && (i<n_frames);
	      i++ )
	{
		fprintf(fh, "%s", filenames[i]);
		if ( events[i] != NULL ) {
			fprintf(fh, " %s\n", events[i]);
		} else {
			fprintf(fh, "\n");
		}
	}

	fclose(fh);
	g_free(file_path);
	g_object_unref(file);
}


static void *run_indexing(const char *job_title,
                          const char *job_notes,
                          struct crystfelproject *proj,
                          void *opts_priv)
{
	struct slurm_indexing_opts *opts = opts_priv;
	struct slurm_job *job;
	char *workdir;
	struct stat s;
	GFile *cwd_file;
	GFile *notes_file;
	GFile *workdir_file;
	char *notes_path;
	FILE *fh;
	char **env;
	int n_env;
	int i;
	int fail = 0;
	char **streams;

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

	workdir = g_file_get_path(workdir_file);

	env = create_env(&n_env, opts->path_add);

	job = malloc(sizeof(struct slurm_job));
	if ( job == NULL ) return 0;

	job->n_frames = proj->n_frames;
	job->n_blocks = proj->n_frames / opts->block_size;
	if ( proj->n_frames % opts->block_size ) job->n_blocks++;
	STATUS("Splitting job into %i blocks of max %i frames\n",
	       job->n_blocks, opts->block_size);

	job->job_ids = malloc(job->n_blocks * sizeof(uint32_t));
	if ( job->job_ids == NULL ) return NULL;

	job->stderr_filenames = malloc(job->n_blocks * sizeof(char *));
	if ( job->stderr_filenames == NULL ) return NULL;

	streams = malloc(job->n_blocks*sizeof(char *));
	if ( streams == NULL ) return NULL;

	for ( i=0; i<job->n_blocks; i++ ) {

		char job_name[128];
		char file_list[128];
		char stream_filename[128];
		char stderr_file[128];
		char stdout_file[128];
		int job_id;
		GFile *stderr_gfile;
		GFile *stream_gfile;

		snprintf(job_name, 127, "%s-%i", job_title, i);
		snprintf(file_list, 127, "files-%i.lst", i);
		snprintf(stream_filename, 127,
		         "crystfel-%i.stream", i);
		snprintf(stderr_file, 127, "stderr-%i.log", i);
		snprintf(stdout_file, 127, "stdout-%i.log", i);

		write_partial_file_list(workdir_file, file_list,
		                        i, opts->block_size,
		                        proj->filenames,
		                        proj->events,
		                        proj->n_frames);

		job_id = submit_batch_job(proj->geom_filename,
		                          file_list,
		                          stream_filename,
		                          opts->email_address,
		                          opts->partition,
		                          env,
		                          n_env,
		                          job_name,
		                          workdir,
		                          stderr_file,
		                          stdout_file,
		                          &proj->peak_search_params,
		                          &proj->indexing_params);

		if ( job_id == 0 ) {
			fail = 1;
			break;
		}

		job->job_ids[i] = job_id;

		stderr_gfile = g_file_get_child(workdir_file,
		                                stderr_file);
		job->stderr_filenames[i] = g_file_get_path(stderr_gfile);
		g_object_unref(stderr_gfile);

		stream_gfile = g_file_get_child(workdir_file,
		                                stream_filename);
		streams[i] = g_file_get_path(stream_gfile);
		g_object_unref(stream_gfile);

		STATUS("Submitted SLURM job ID %i\n", job_id);
	}

	for ( i=0; i<n_env; i++ ) free(env[i]);
	free(env);
	free(workdir);
	g_object_unref(workdir_file);

	if ( fail ) {
		free(job->job_ids);
		free(job->stderr_filenames);
		free(job);
		return NULL;
	} else {
		add_result(proj, strdup(job_title),
		           streams, job->n_blocks);
	}
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


static void pathadd_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_indexing_opts *opts = data;
	opts->path_add = strdup(gtk_entry_get_text(entry));
}


static gboolean pathadd_focus_sig(GtkEntry *entry, GdkEvent *event,
                                  gpointer data)
{
	pathadd_activate_sig(entry, data);
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

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Search path for executables:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->path_add != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->path_add);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "/path/to/indexing/programs");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(pathadd_activate_sig),
	                 opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(pathadd_focus_sig),
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
	opts->path_add = NULL;

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

	if ( opts->path_add != NULL ) {
		fprintf(fh, "indexing.slurm.path_add %s\n",
		        opts->path_add);
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

	if ( strcmp(key, "indexing.slurm.path_add") == 0 ) {
		opts->path_add = strdup(val);
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
