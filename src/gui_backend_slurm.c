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
#include <slurm/slurm.h>
#include <gio/gio.h>

#include <utils.h>

#include "gtk-util-routines.h"
#include "gui_project.h"
#include "gui_index.h"
#include "gui_merge.h"
#include "gui_ambi.h"
#include "crystfel_gui.h"


struct slurm_common_opts
{
	char *partition;
	char *email_address;
	char *account;
	char *constraint;
};


struct slurm_indexing_opts
{
	struct slurm_common_opts common;
	char *path_add;
	int block_size;
};


struct slurm_merging_opts
{
	struct slurm_common_opts common;
};


struct slurm_ambi_opts
{
	struct slurm_common_opts common;
};


struct slurm_job
{
	enum gui_job_type type;
	GFile *workdir;
	uint32_t job_id;

	/* For indexing job */
	int n_frames;
	int n_blocks;

	/* For merging/ambigator job */
	char *stderr_filename;
	int niter;
};


static int job_alive(slurm_job_info_t *job)
{
	switch ( job->job_state & JOB_STATE_BASE ) {

		/* Only the following states are reasons to keep on watching
		 * the job */
		case JOB_PENDING :
		case JOB_RUNNING :
		case JOB_SUSPENDED :
		return 1;

		default :
		return 0;
	}
}


static int job_running(uint32_t job_id)
{
	job_info_msg_t *job_info;
	int running = 1;

	if ( slurm_load_job(&job_info, job_id, 0) ) {
		STATUS("Couldn't get status: %i\n",
		       slurm_strerror(slurm_get_errno()));
		return 0;
	}

	running = job_alive(&job_info->job_array[0]);
	slurm_free_job_info_msg(job_info);
	return running;
}


static double indexing_progress(struct slurm_job *job, int *running)
{
	job_info_msg_t *array_job_info;
	int i;
	int n_running;
	int lowest_alive_task;

	if ( slurm_load_job(&array_job_info, job->job_id, 0) ) {
		STATUS("Couldn't get status: %i\n",
		       slurm_strerror(slurm_get_errno()));
		*running = 0;
		return 0.0;
	}

	n_running = 0;
	lowest_alive_task = job->n_blocks;
	for ( i=0; i<array_job_info->record_count; i++ ) {

		slurm_job_info_t *job_info = &array_job_info->job_array[i];

		/* Ignore the array job itself */
		if ( job_info->array_job_id == 0 ) continue;

		if ( job_alive(job_info) ) {
			if ( job_info->array_task_id < lowest_alive_task ) {
				lowest_alive_task = job_info->array_task_id;
			}
			n_running++;
		}
	}
	slurm_free_job_info_msg(array_job_info);

	*running = (n_running > 0);

	/* If there are lots of blocks, just count running jobs instead of
	 * reading loads of log files */
	if ( job->n_blocks > 15 ) {

		/* Didn't find any alive jobs at all?
		 * Then we've either just started or just finished. */
		if ( lowest_alive_task == job->n_blocks ) {
			if ( n_running > 0 ) {
				return 0.0;
			} else {
				return 1.0;
			}
		} else {
			return (double)lowest_alive_task / job->n_blocks;
		}

	} else {

		int i;
		int n_proc = 0;

		for ( i=0; i<job->n_blocks; i++ ) {

			char tmp[128];
			GFile *stderr_gfile;
			char *stderr_filename;

			snprintf(tmp, 127, "stderr-%i.log", i);
			stderr_gfile = g_file_get_child(job->workdir, tmp);
			stderr_filename = g_file_get_path(stderr_gfile);
			g_object_unref(stderr_gfile);

			n_proc += read_number_processed(stderr_filename);
			g_free(stderr_filename);

		}

		return (double)n_proc / job->n_frames;
	}
}


static int get_task_status(void *job_priv,
                           int *running,
                           float *frac_complete)
{
	struct slurm_job *job = job_priv;

	switch ( job->type ) {

		case GUI_JOB_INDEXING :
		*frac_complete = indexing_progress(job, running);
		break;

		case GUI_JOB_AMBIGATOR :
		*frac_complete = read_ambigator_progress(job->stderr_filename,
		                                         job->niter);
		*running = job_running(job->job_id);
		break;

		case GUI_JOB_PROCESS_HKL :
		case GUI_JOB_PROCESS_HKL_SCALE :
		case GUI_JOB_PARTIALATOR :
		*frac_complete = read_merge_progress(job->stderr_filename,
		                                     job->type);
		*running = job_running(job->job_id);
		break;

	}

	return 0;
}


static void cancel_task(void *job_priv)
{
	char jobid[128];
	const gchar *args[3];
	GError *error = NULL;
	GSubprocess *sp;
	struct slurm_job *job = job_priv;

	snprintf(jobid, 127, "%i", job->job_id);
	args[0] = "scancel";
	args[1] = jobid;
	args[2] = NULL;
	sp = g_subprocess_newv(args, G_SUBPROCESS_FLAGS_NONE, &error);
	if ( sp == NULL ) {
		ERROR("Failed to invoke scancel: %s\n", error->message);
		g_error_free(error);
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


static void partition_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->partition = safe_strdup(get_text_or_null(entry));
}


static gboolean partition_focus_sig(GtkEntry *entry, GdkEvent *event,
                                    gpointer data)
{
	partition_activate_sig(entry, data);
	return FALSE;
}


static void email_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->email_address = safe_strdup(get_text_or_null(entry));
}


static gboolean email_focus_sig(GtkEntry *entry, GdkEvent *event,
                                gpointer data)
{
	email_activate_sig(entry, data);
	return FALSE;
}


static void account_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->account = safe_strdup(get_text_or_null(entry));
}


static gboolean account_focus_sig(GtkEntry *entry, GdkEvent *event,
                                  gpointer data)
{
	account_activate_sig(entry, data);
	return FALSE;
}


static void constraint_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->constraint = safe_strdup(get_text_or_null(entry));
}


static gboolean constraint_focus_sig(GtkEntry *entry, GdkEvent *event,
                                     gpointer data)
{
	constraint_activate_sig(entry, data);
	return FALSE;
}


static void add_common_opts(GtkWidget *vbox,
                            struct slurm_common_opts *opts)
{
	GtkWidget *hbox;
	GtkWidget *entry;
	GtkWidget *label;

	/* Partition */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Submit job to partition:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->partition != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->partition);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry), "maxwell");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(partition_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(partition_focus_sig), opts);

	/* Email address */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Send notifications to:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->email_address != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->email_address);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "myself@example.org");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(email_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(email_focus_sig), opts);

	/* Account */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Charge resource use to:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->account != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->account);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "SLURM account");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(account_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(account_focus_sig), opts);

	/* Constraint */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Required node features:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->constraint != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->constraint);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "SLURM constraint");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(constraint_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(constraint_focus_sig), opts);
}


static void write_common_opts(FILE *fh,
                              struct slurm_common_opts *opts,
                              const char *prefix)
{
	if ( opts->partition != NULL) {
		fprintf(fh, "%s.slurm.partition %s\n",
		        prefix, opts->partition);
	}

	if ( opts->email_address != NULL ) {
		fprintf(fh, "%s.slurm.email_address %s\n",
		        prefix, opts->email_address);
	}

	if ( opts->account != NULL ) {
		fprintf(fh, "%s.slurm.account %s\n",
		        prefix, opts->account);
	}

	if ( opts->constraint != NULL ) {
		fprintf(fh, "%s.slurm.constraint %s\n",
		        prefix, opts->constraint);
	}
}


static struct slurm_job *start_slurm_job(enum gui_job_type type,
                                         const char *script_filename,
                                         const char *jobname,
                                         const char *array_inx,
                                         GFile *workdir,
                                         const char *stdout_filename,
                                         const char *stderr_filename,
                                         struct slurm_common_opts *opts)
{
	char **env;
	int n_env;
	char *script;
	struct slurm_job *job;
	job_desc_msg_t job_desc_msg;
	submit_response_msg_t *resp;
	GFile *stderr_gfile;
	int r;

	script = load_entire_file(script_filename);
	if ( script == NULL ) return NULL;

	job = malloc(sizeof(struct slurm_job));
	if ( job == NULL ) {
		free(script);
		return NULL;
	}

	job->type = type;

	env = create_env(&n_env, NULL);

	slurm_init_job_desc_msg(&job_desc_msg);
	job_desc_msg.user_id = getuid();
	job_desc_msg.group_id = getgid();
	job_desc_msg.mail_user = safe_strdup(opts->email_address);
	job_desc_msg.mail_type = MAIL_JOB_FAIL;
	job_desc_msg.comment = "Submitted via CrystFEL GUI";
	job_desc_msg.shared = 0;
	job_desc_msg.time_limit = 60;
	job_desc_msg.partition = safe_strdup(opts->partition);
	job_desc_msg.min_nodes = 1;
	job_desc_msg.max_nodes = 1;
	job_desc_msg.name = safe_strdup(jobname);
	job_desc_msg.std_err = strdup(stderr_filename);
	job_desc_msg.std_out = strdup(stdout_filename);
	job_desc_msg.work_dir = g_file_get_path(workdir);
	job_desc_msg.script = script;
	job_desc_msg.environment = env;
	job_desc_msg.env_size = n_env;
	job_desc_msg.features = safe_strdup(opts->constraint);
	job_desc_msg.account = safe_strdup(opts->account);
	job_desc_msg.array_inx = safe_strdup(array_inx);

	r = slurm_submit_batch_job(&job_desc_msg, &resp);
	free(job_desc_msg.mail_user);
	free(job_desc_msg.partition);
	free(job_desc_msg.name);
	free(job_desc_msg.work_dir);
	free(job_desc_msg.std_err);
	free(job_desc_msg.std_out);
	free(job_desc_msg.features);
	free(job_desc_msg.account);
	free(job_desc_msg.script);
	if ( r ) {
		ERROR("Couldn't submit job: %i\n", errno);
		free(job);
		return NULL;
	}

	STATUS("Submitted SLURM job ID %i\n", resp->job_id);

	job->job_id = resp->job_id;
	slurm_free_submit_response_response_msg(resp);

	stderr_gfile = g_file_get_child(workdir, stderr_filename);
	job->stderr_filename = g_file_get_path(stderr_gfile);
	g_object_unref(stderr_gfile);

	job->workdir = g_file_dup(workdir);

	return job;
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
	int i;
	char **streams;
	GFile *workdir;
	GFile *sc_gfile;
	char *sc_filename;
	int n_blocks;
	char array_inx[128];
	char serial_offs[128];

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	n_blocks = proj->n_frames / opts->block_size;
	if ( proj->n_frames % opts->block_size ) n_blocks++;
	STATUS("Splitting job into %i blocks of max %i frames\n",
	       n_blocks, opts->block_size);

	streams = malloc(n_blocks*sizeof(char *));
	if ( streams == NULL ) return NULL;

	for ( i=0; i<n_blocks; i++ ) {

		char file_list[128];
		char stream_filename[128];
		GFile *stream_gfile;

		/* Create (sub-)list of files */
		snprintf(file_list, 127, "files-%i.lst", i);
		write_partial_file_list(workdir,
		                        file_list,
		                        i,
		                        opts->block_size,
		                        proj->filenames,
		                        proj->events,
		                        proj->n_frames);

		/* Work out the stream filename */
		snprintf(stream_filename, 127, "crystfel-%i.stream", i);
		stream_gfile = g_file_get_child(workdir,
		                                stream_filename);
		streams[i] = g_file_get_path(stream_gfile);
		g_object_unref(stream_gfile);
	}

	sc_gfile = g_file_get_child(workdir, "run_indexamajig.sh");
	sc_filename = g_file_get_path(sc_gfile);
	g_object_unref(sc_gfile);
	if ( sc_filename == NULL ) return NULL;

	snprintf(array_inx, 127, "0-%i", n_blocks-1);
	snprintf(serial_offs, 127, "$((${SLURM_ARRAY_TASK_ID}*%i+1))",
	         opts->block_size);

	if ( !write_indexamajig_script(sc_filename,
	                               proj->geom_filename,
	                               "`nproc`",
	                               "files-${SLURM_ARRAY_TASK_ID}.lst",
	                               "crystfel-${SLURM_ARRAY_TASK_ID}.stream",
	                               serial_offs, 0,
	                               &proj->peak_search_params,
	                               &proj->indexing_params) )
	{
		job = start_slurm_job(GUI_JOB_INDEXING,
		                      sc_filename,
		                      job_title,
		                      array_inx,
		                      workdir,
		                      "stdout-%a.log",
		                      "stderr-%a.log",
		                      &opts->common);
	} else {
		job = NULL;
	}
	g_free(sc_filename);

	if ( job != NULL ) {
		job->n_frames = proj->n_frames;
		job->n_blocks = n_blocks;
		add_indexing_result(proj, strdup(job_title), streams, n_blocks);
	}

	g_object_unref(workdir);

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

	add_common_opts(vbox, &opts->common);

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


static void set_default_common_opts(struct slurm_common_opts *opts)
{
	opts->partition = NULL;
	opts->email_address = NULL;
	opts->account = NULL;
	opts->constraint = NULL;
}


static struct slurm_indexing_opts *make_default_slurm_indexing_opts()
{
	struct slurm_indexing_opts *opts = malloc(sizeof(struct slurm_indexing_opts));
	if ( opts == NULL ) return NULL;

	set_default_common_opts(&opts->common);
	opts->block_size = 1000;
	opts->path_add = NULL;

	return opts;
}


static void write_indexing_opts(void *opts_priv, FILE *fh)
{
	struct slurm_indexing_opts *opts = opts_priv;

	write_common_opts(fh, &opts->common, "indexing");

	fprintf(fh, "indexing.slurm.block_size %i\n",
	        opts->block_size);

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
		opts->common.email_address = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.partition") == 0 ) {
		opts->common.partition = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.path_add") == 0 ) {
		opts->path_add = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.account") == 0 ) {
		opts->common.account = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.constraint") == 0 ) {
		opts->common.constraint = strdup(val);
	}
}


static void *run_ambi(const char *job_title,
                      const char *job_notes,
                      struct crystfelproject *proj,
                      struct gui_indexing_result *input,
                      void *opts_priv)
{
	struct slurm_job *job;
	struct slurm_ambi_opts *opts = opts_priv;
	GFile *workdir;
	GFile *sc_gfile;
	char *sc_filename;
	GFile *stream_gfile;
	char *stream_str;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	stream_gfile = g_file_get_child(workdir, "ambi.stream");
	stream_str = g_file_get_path(stream_gfile);
	g_object_unref(stream_gfile);

	sc_gfile = g_file_get_child(workdir, "run_ambigator.sh");
	sc_filename = g_file_get_path(sc_gfile);
	g_object_unref(sc_gfile);
	if ( sc_filename == NULL ) return NULL;

	if ( !write_ambigator_script(sc_filename, input, "`nproc`",
	                             &proj->ambi_params, stream_str) )
	{
		job = start_slurm_job(GUI_JOB_AMBIGATOR,
		                      sc_filename, job_title, NULL, workdir,
		                      "stdout.log", "stderr.log", &opts->common);
		job->niter = proj->ambi_params.niter;
	} else {
		job = NULL;
	}
	g_free(sc_filename);

	if ( job != NULL ) {
		char **streams = malloc(sizeof(char *));
		if ( streams != NULL ) {
			streams[0] = stream_str;
			add_indexing_result(proj, strdup(job_title), streams, 1);
		}
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
	struct slurm_job *job;
	struct slurm_merging_opts *opts = opts_priv;
	GFile *workdir;
	GFile *sc_gfile;
	char *sc_filename;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	sc_gfile = g_file_get_child(workdir, "run_merge.sh");
	sc_filename = g_file_get_path(sc_gfile);
	g_object_unref(sc_gfile);
	if ( sc_filename == NULL ) return NULL;

	if ( !write_merge_script(sc_filename, input, "`nproc`",
	                         &proj->merging_params, "crystfel.hkl") )
	{
		enum gui_job_type type;
		if ( strcmp(proj->merging_params.model, "process_hkl") == 0 ) {
			if ( proj->merging_params.scale ) {
				type = GUI_JOB_PROCESS_HKL_SCALE;
			} else {
				type = GUI_JOB_PROCESS_HKL;
			}
		} else {
			type = GUI_JOB_PARTIALATOR;
		}
		job = start_slurm_job(type, sc_filename, job_title, NULL,
		                      workdir, "stdout.log", "stderr.log",
		                      &opts->common);
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

		add_merge_result(proj, strdup(job_title), hkl, hkl1, hkl2);
	}

	g_object_unref(workdir);
	return job;
}


static struct slurm_merging_opts *make_default_slurm_merging_opts()
{
	struct slurm_merging_opts *opts = malloc(sizeof(struct slurm_merging_opts));
	if ( opts == NULL ) return NULL;
	set_default_common_opts(&opts->common);
	return opts;
}


static void write_merging_opts(void *opts_priv, FILE *fh)
{
	struct slurm_merging_opts *opts = opts_priv;
	write_common_opts(fh, &opts->common, "merging");
}


static void read_merging_opt(void *opts_priv,
                             const char *key,
                             const char *val)
{
	struct slurm_merging_opts *opts = opts_priv;

	if ( strcmp(key, "merging.slurm.email_address") == 0 ) {
		opts->common.email_address = strdup(val);
	}

	if ( strcmp(key, "merging.slurm.partition") == 0 ) {
		opts->common.partition = strdup(val);
	}

	if ( strcmp(key, "merging.slurm.account") == 0 ) {
		opts->common.account = strdup(val);
	}

	if ( strcmp(key, "merging.slurm.constraint") == 0 ) {
		opts->common.constraint = strdup(val);
	}

}


static GtkWidget *make_merging_parameters_widget(void *opts_priv)
{
	GtkWidget *vbox;
	struct slurm_merging_opts *opts = opts_priv;
	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	add_common_opts(vbox, &opts->common);
	return vbox;
}


static struct slurm_ambi_opts *make_default_slurm_ambi_opts()
{
	struct slurm_ambi_opts *opts = malloc(sizeof(struct slurm_ambi_opts));
	if ( opts == NULL ) return NULL;
	set_default_common_opts(&opts->common);
	return opts;
}


static void write_ambi_opts(void *opts_priv, FILE *fh)
{
	struct slurm_ambi_opts *opts = opts_priv;
	write_common_opts(fh, &opts->common, "ambi");
}


static void read_ambi_opt(void *opts_priv,
                          const char *key,
                          const char *val)
{
	struct slurm_ambi_opts *opts = opts_priv;

	if ( strcmp(key, "ambi.slurm.email_address") == 0 ) {
		opts->common.email_address = strdup(val);
	}

	if ( strcmp(key, "ambi.slurm.partition") == 0 ) {
		opts->common.partition = strdup(val);
	}

	if ( strcmp(key, "ambi.slurm.account") == 0 ) {
		opts->common.account = strdup(val);
	}

	if ( strcmp(key, "ambi.slurm.constraint") == 0 ) {
		opts->common.constraint = strdup(val);
	}
}


static GtkWidget *make_ambi_parameters_widget(void *opts_priv)
{
	GtkWidget *vbox;
	struct slurm_ambi_opts *opts = opts_priv;
	vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	add_common_opts(vbox, &opts->common);
	return vbox;
}


int make_slurm_backend(struct crystfel_backend *be)
{
	be->name = "slurm";
	be->friendly_name = "SLURM";

	be->cancel_task = cancel_task;
	be->task_status = get_task_status;

	be->make_indexing_parameters_widget = make_indexing_parameters_widget;
	be->run_indexing = run_indexing;
	be->indexing_opts_priv = make_default_slurm_indexing_opts();
	if ( be->indexing_opts_priv == NULL ) return 1;
	be->write_indexing_opts = write_indexing_opts;
	be->read_indexing_opt = read_indexing_opt;

	be->make_merging_parameters_widget = make_merging_parameters_widget;
	be->run_merging = run_merging;
	be->merging_opts_priv = make_default_slurm_merging_opts();
	if ( be->merging_opts_priv == NULL ) return 1;
	be->write_merging_opts = write_merging_opts;
	be->read_merging_opt = read_merging_opt;

	be->make_ambi_parameters_widget = make_ambi_parameters_widget;
	be->run_ambi = run_ambi;
	be->ambi_opts_priv = make_default_slurm_ambi_opts();
	if ( be->ambi_opts_priv == NULL ) return 1;
	be->write_ambi_opts = write_ambi_opts;
	be->read_ambi_opt = read_ambi_opt;

	return 0;
};
