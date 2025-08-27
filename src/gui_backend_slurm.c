/*
 * gui_backend_slurm.c
 *
 * GUI backend for running jobs via SLURM
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
#include <gtk/gtk.h>
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
	char *reservation;
	char *qos;
	char *partition;
	char *email_address;
	char *account;
	char *constraint;
	int time_limit;
	int exclusive;
};


struct slurm_indexing_opts
{
	struct slurm_common_opts common;
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


static const char *get_str_val(const char *line, const char *key)
{
	const char *pos = strstr(line, key);
	if ( pos == NULL ) return NULL;

	const char *eq = strchr(pos, '=');
	if ( eq == NULL ) return NULL;

	const char *sp = strchr(pos, ' ');
	if ( sp == NULL ) return NULL;

	return strndup(eq+1, sp-eq-1);
}


static char *g_bytes_to_terminated_array(GBytes *bytes)
{
	gpointer arr;
	gsize size;
	char *buf;

	arr = g_bytes_unref_to_data(bytes, &size);

	buf = malloc(size+1);
	if ( buf == NULL ) return NULL;

	memcpy(buf, arr, size);
	buf[size] = '\0';

	g_free(arr);

	return buf;
}


static int get_job_status(int job_id, int *n_alive, int *n_running)
{
	const gchar *args[6];
	GError *error = NULL;
	GSubprocess *sp;
	char job_id_str[64];
	char *line;
	char *nl;
	GBytes *stdout_buf;
	GBytes *stderr_buf;
	char *buf;
	char *buf_stderr;

	snprintf(job_id_str, 63, "%i", job_id);
	args[0] = "scontrol";
	args[1] = "-o";
	args[2] = "show";
	args[3] = "job";
	args[4] = job_id_str;
	args[5] = NULL;

	sp = g_subprocess_newv(args, G_SUBPROCESS_FLAGS_STDOUT_PIPE
	                           | G_SUBPROCESS_FLAGS_STDERR_PIPE, &error);
	if ( sp == NULL ) {
		ERROR("Failed to invoke scontrol: %s\n", error->message);
		g_error_free(error);
		return 1;
	}

	if ( !g_subprocess_communicate(sp, NULL, NULL,
	                               &stdout_buf, &stderr_buf, &error) )
	{
		ERROR("Error communicating with scontrol: %s\n", error->message);
		g_error_free(error);
		return 1;
	}

	buf = g_bytes_to_terminated_array(stdout_buf);
	buf_stderr = g_bytes_to_terminated_array(stderr_buf);

	if ( buf_stderr[0] != '\0' ) {
		ERROR("scontrol error: %s\n", buf_stderr);
		/* ... but carry on */
	}
	free(buf_stderr);

	*n_alive = 0;
	*n_running = 0;

	/* Parse output */
	line = &buf[0];
	nl = strchr(line, '\n');
	while ( nl != NULL ) {

		int p1, p2;

		nl[0] = '\0';

		const char *state = get_str_val(line, "JobState");
		const char *array_task_str = get_str_val(line, "ArrayTaskId");

		if ((strcmp(state, "PENDING") == 0)
		 || (strcmp(state, "SUSPENDED") == 0))
		{
			(*n_alive)++;
		}

		if ((strcmp(state, "RUNNING") == 0)
		 || (strcmp(state, "COMPLETING") == 0))
		{
			(*n_running)++;
		}

		if ( (array_task_str != NULL)
		  && (sscanf(array_task_str, "%i-%i", &p1, &p2) == 2) )
		{
			/* This is a "job array leader" */
			if ((strcmp(state, "PENDING") == 0)
			 || (strcmp(state, "SUSPENDED") == 0)) {
				(*n_alive) += p2-p1;
			}
		}

		/* We are not interested in: FAILED, COMPLETED, CANCELLED */

		line = nl+1;
		nl = strchr(line, '\n');
	}

	free(buf);

	return 0;
}


static double indexing_progress(struct slurm_job *job, int n_alive, int n_running)
{
	/* If there are lots of blocks, just count running jobs instead of
	 * reading loads of log files */
	if ( job->n_blocks > 15 ) {

		return (job->n_blocks - n_alive - 0.5*n_running) / job->n_blocks;

	} else {

		int ijob;
		int n_proc = 0;

		for ( ijob=0; ijob<job->n_blocks; ijob++ ) {

			char tmp[128];
			char *stderr_filename;

			snprintf(tmp, 127, "stderr-%i.log", ijob);
			stderr_filename = relative_to_cwd(job->workdir, tmp);

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
	int n_running, n_alive;

	if ( get_job_status(job->job_id, &n_alive, &n_running) ) {
		ERROR("Failed to get task status: %i\n", job->job_id);
		return 1;
	}

	switch ( job->type ) {

		case GUI_JOB_INDEXING :
		*frac_complete = indexing_progress(job, n_alive, n_running);
		*running = (n_alive+n_running > 0);
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


static void free_task(void *job_priv)
{
	struct slurm_job *job = job_priv;
	g_object_unref(job->workdir);
	free(job->stderr_filename);
}


static void cancel_task(void *job_priv)
{
	char jobid[128];
	const gchar *args[3];
	GError *error = NULL;
	GSubprocess *sp;
	struct slurm_job *job = job_priv;

	snprintf(jobid, 127, "%u", job->job_id);
	args[0] = "scancel";
	args[1] = jobid;
	args[2] = NULL;
	sp = g_subprocess_newv(args, G_SUBPROCESS_FLAGS_NONE, &error);
	if ( sp == NULL ) {
		ERROR("Failed to invoke scancel: %s\n", error->message);
		g_error_free(error);
	}
}


static void reservation_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->reservation = safe_strdup(get_text_or_null(entry));
}


static gboolean reservation_focus_sig(GtkEntry *entry, GdkEvent *event,
                                      gpointer data)
{
	reservation_activate_sig(entry, data);
	return FALSE;
}


static void qos_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->qos = safe_strdup(get_text_or_null(entry));
}


static gboolean qos_focus_sig(GtkEntry *entry, GdkEvent *event,
                              gpointer data)
{
	qos_activate_sig(entry, data);
	return FALSE;
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


static void timelimit_activate_sig(GtkEntry *entry, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->time_limit = get_uint(GTK_WIDGET(entry));
}


static gboolean timelimit_focus_sig(GtkEntry *entry, GdkEvent *event,
                                    gpointer data)
{
	timelimit_activate_sig(entry, data);
	return FALSE;
}


static void exclusive_toggle_sig(GtkToggleButton *toggle, gpointer data)
{
	struct slurm_common_opts *opts = data;
	opts->exclusive = gtk_toggle_button_get_active(toggle);
}


static void add_common_opts(GtkWidget *vbox,
                            struct slurm_common_opts *opts)
{
	GtkWidget *hbox;
	GtkWidget *entry;
	GtkWidget *label;
	GtkWidget *toggle;
	char tmp[64];

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

	/* Time limit */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Job time limit (minutes):");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	snprintf(tmp, 63, "%i", opts->time_limit);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "time limit in minutes");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(timelimit_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(timelimit_focus_sig), opts);

	/* Reservation */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Reservation:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->reservation != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->reservation);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "SLURM reservation");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(reservation_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(reservation_focus_sig), opts);

	/* QoS */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	label = gtk_label_new("Quality of service:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0);
	entry = gtk_entry_new();
	if ( opts->qos != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(entry), opts->qos);
	}
	gtk_entry_set_placeholder_text(GTK_ENTRY(entry),
	                               "SLURM QoS");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(entry), "activate",
	                 G_CALLBACK(qos_activate_sig), opts);
	g_signal_connect(G_OBJECT(entry), "focus-out-event",
	                 G_CALLBACK(qos_focus_sig), opts);

	/* Exclusive */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 0);
	toggle = gtk_check_button_new_with_label("Request exclusive use of compute node(s)");
	set_active(toggle, opts->exclusive);
	gtk_widget_set_tooltip_text(toggle, "--exclusive");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(toggle), FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(toggle), "toggled",
	                 G_CALLBACK(exclusive_toggle_sig), opts);
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

	if ( opts->reservation != NULL ) {
		fprintf(fh, "%s.slurm.reservation %s\n",
		        prefix, opts->reservation);
	}

	if ( opts->qos != NULL ) {
		fprintf(fh, "%s.slurm.qos %s\n",
		        prefix, opts->qos);
	}

	fprintf(fh, "%s.slurm.time_limit %i\n",
	        prefix, opts->time_limit);

	fprintf(fh, "%s.slurm.exclusive %i\n",
	        prefix, opts->exclusive);
}


static int empty(const char *str)
{
	if ( str == NULL ) return 1;
	if ( str[0] == '\0' ) return 1;
	return 0;
}


static char *add_bits(char *old, const char *new1, const char *new2)
{
	size_t len;
	char *nn;

	if ( old == NULL ) return NULL;

	len = strlen(new1) + strlen(new2) + strlen(old) + 12;
	nn = malloc(len);
	if ( nn == NULL ) {
		ERROR("Failed to expand string for #SBATCH bits.\n");
		return NULL;
	}

	strcpy(nn, old);
	strcat(nn, "\n#SBATCH ");
	strcat(nn, new1);
	if ( strlen(new2) > 0 ) {
		strcat(nn, " ");
		strcat(nn, new2);
	}
	free(old);
	return nn;
}


static char *sbatch_bits(struct slurm_common_opts *opts,
                         const char *jobname,
                         const char *array_inx,
                         const char *stdout_filename,
                         const char *stderr_filename)
{
	char time_limit[64];
	char *str = strdup("");

	if ( !empty(array_inx) ) {
		str = add_bits(str, "--array", array_inx);
	}
	str = add_bits(str, "--job-name", jobname);
	str = add_bits(str, "--output", stdout_filename);
	str = add_bits(str, "--error", stderr_filename);
	snprintf(time_limit, 63, "%i", opts->time_limit);
	str = add_bits(str, "--time", time_limit);
	if ( !empty(opts->email_address) ) {
		str = add_bits(str, "--mail-user", opts->email_address);
	}
	if ( !empty(opts->partition) ) {
		str = add_bits(str, "--partition", opts->partition);
	}
	if ( !empty(opts->constraint) ) {
		str = add_bits(str, "--constraint", opts->constraint);
	}
	if ( !empty(opts->account) ) {
		str = add_bits(str, "--account", opts->account);
	}
	if ( !empty(opts->reservation) ) {
		str = add_bits(str, "--reservation", opts->reservation);
	}
	if ( !empty(opts->qos) ) {
		str = add_bits(str, "--qos", opts->qos);
	}
	if ( opts->exclusive ) {
		str = add_bits(str, "--exclusive", "");
	}
	str = add_bits(str, "--nodes", "1");
	str = add_bits(str, "--mail-type", "FAIL\n\n");

	return str;
}


static struct slurm_job *start_slurm_job(enum gui_job_type type,
                                         const char *script_filename,
                                         const char *jobname,
                                         GFile *workdir,
                                         const char *stderr_filename)
{
	const gchar *args[5];
	GError *error = NULL;
	GSubprocess *sp;
	char buf[256];
	gsize bytes_read;

	args[0] = "sbatch";
	args[1] = "--comment";
	args[2] = "Submitted via CrystFEL GUI";
	args[3] = script_filename;
	args[4] = NULL;

	sp = g_subprocess_newv(args, G_SUBPROCESS_FLAGS_STDOUT_PIPE
	                           | G_SUBPROCESS_FLAGS_STDERR_MERGE, &error);
	if ( sp == NULL ) {
		ERROR("Failed to invoke sbatch: %s\n", error->message);
		g_error_free(error);
		return NULL;
	}

	g_subprocess_wait(sp, NULL, &error);

	bytes_read = g_input_stream_read(g_subprocess_get_stdout_pipe(sp),
	                                 buf, 256, NULL, &error);
	buf[bytes_read] = '\0';
	chomp(buf);

	if ( strncmp(buf, "Submitted batch job ", 20) == 0 ) {

		struct slurm_job *job;
		int job_id;

		if ( convert_int(buf+20, &job_id) ) {
			ERROR("Didn't get batch job ID from '%s'\n", buf);
			return NULL;
		}

		job = malloc(sizeof(struct slurm_job));
		if ( job == NULL ) return NULL;

		job->type = type;
		job->job_id = job_id;
		job->stderr_filename = strdup(stderr_filename);
		job->workdir = g_file_dup(workdir);

		STATUS("Submitted batch job ID %i\n", job_id);

		return job;

	} else {
		ERROR("Didn't understand sbatch reply: '%s'\n", buf);
		return NULL;
	}
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
	if ( fh == NULL ) {
		ERROR("Failed to write %s\n", file_path);
		g_free(file_path);
		g_object_unref(file);
		return;
	}

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
                          void *opts_priv,
                          double wavelength_estimate,
                          double clen_estimate)
{
	struct slurm_indexing_opts *opts = opts_priv;
	struct slurm_job *job;
	int i;
	char **streams;
	GFile *workdir;
	int n_blocks;
	char array_inx[128];
	char serial_offs[128];
	gchar *sc_rel_filename;
	gchar *stdout_rel_filename;
	gchar *stderr_rel_filename;
	gchar *files_rel_filename;
	gchar *stream_rel_filename;
	gchar *harvest_rel_filename;
	gchar *mille_rel_filename;
	char *slurm_prologue;
	GFile *ggeom;
	GFile *ggeomcopy;
	GError *error;

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
		streams[i] = relative_to_cwd(workdir, stream_filename);
	}

	snprintf(array_inx, 127, "0-%i", n_blocks-1);
	snprintf(serial_offs, 127, "$((${SLURM_ARRAY_TASK_ID}*%i+1))",
	         opts->block_size);

	sc_rel_filename = relative_to_cwd(workdir, "run_indexamajig.sh");
	files_rel_filename = relative_to_cwd(workdir,
	                                     "files-${SLURM_ARRAY_TASK_ID}.lst");
	stream_rel_filename = relative_to_cwd(workdir,
	                                      "crystfel-${SLURM_ARRAY_TASK_ID}.stream");
	stdout_rel_filename = relative_to_cwd(workdir, "stdout-%a.log");
	stderr_rel_filename = relative_to_cwd(workdir, "stderr-%a.log");
	harvest_rel_filename = relative_to_cwd(workdir, "parameters.json");
	mille_rel_filename = relative_to_cwd(workdir, "mille-data-${SLURM_ARRAY_TASK_ID}");

	slurm_prologue = sbatch_bits(&opts->common, job_title, array_inx,
	                             stdout_rel_filename, stderr_rel_filename);

	/* Copy geometry file into working directory
	 * Used for geometry refinement, not indexing! */
	ggeom = g_file_new_for_path(proj->geom_filename);
	ggeomcopy = g_file_get_child(workdir, "detector.geom");
	error = NULL;
	g_file_copy(ggeom, ggeomcopy, G_FILE_COPY_BACKUP | G_FILE_COPY_ALL_METADATA,
	            NULL, NULL, NULL, &error);
	g_object_unref(ggeom);
	g_object_unref(ggeomcopy);

	if ( !write_indexamajig_script(sc_rel_filename,
	                               proj->geom_filename,
	                               "$(nproc)",
	                               files_rel_filename,
	                               stream_rel_filename,
	                               NULL, NULL,
	                               harvest_rel_filename,
	                               mille_rel_filename,
	                               serial_offs,
	                               &proj->peak_search_params,
	                               &proj->indexing_params,
	                               wavelength_estimate,
	                               clen_estimate,
	                               slurm_prologue) )
	{
		/* The stderr filename isn't used by indexing_progress() in the
		 * Slurm backend - it knows where to find the files. */
		job = start_slurm_job(GUI_JOB_INDEXING,
		                      sc_rel_filename,
		                      job_title, workdir, "notused");
	} else {
		job = NULL;
	}

	free(slurm_prologue);

	if ( job != NULL ) {
		job->n_frames = proj->n_frames;
		job->n_blocks = n_blocks;
		add_indexing_result(proj, job_title, streams, n_blocks);
	}

	for ( i=0; i<n_blocks; i++ ) {
		free(streams[i]);
	}
	free(streams);

	free(sc_rel_filename);
	free(files_rel_filename);
	free(stream_rel_filename);
	free(stdout_rel_filename);
	free(stderr_rel_filename);
	free(harvest_rel_filename);
	free(mille_rel_filename);
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

	return vbox;
}


static void set_default_common_opts(struct slurm_common_opts *opts)
{
	opts->partition = NULL;
	opts->email_address = NULL;
	opts->account = NULL;
	opts->constraint = NULL;
	opts->time_limit = 60;
	opts->qos = NULL;
	opts->exclusive = 1;
	opts->reservation = NULL;
}


static struct slurm_indexing_opts *make_default_slurm_indexing_opts()
{
	struct slurm_indexing_opts *opts = malloc(sizeof(struct slurm_indexing_opts));
	if ( opts == NULL ) return NULL;

	set_default_common_opts(&opts->common);
	opts->block_size = 1000;

	return opts;
}


static void write_indexing_opts(void *opts_priv, FILE *fh)
{
	struct slurm_indexing_opts *opts = opts_priv;

	write_common_opts(fh, &opts->common, "indexing");

	fprintf(fh, "indexing.slurm.block_size %i\n",
	        opts->block_size);
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
		STATUS("The 'search path for executables' input for the SLURM "
		       "backend is no longer used.\n");
		STATUS("Add the correct locations to PATH before starting the "
		       "CrystFEL GUI.  This will be propagated to your SLURM "
		       "jobs.\n");
	}

	if ( strcmp(key, "indexing.slurm.account") == 0 ) {
		opts->common.account = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.constraint") == 0 ) {
		opts->common.constraint = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.reservation") == 0 ) {
		opts->common.reservation = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.qos") == 0 ) {
		opts->common.qos = strdup(val);
	}

	if ( strcmp(key, "indexing.slurm.time_limit") == 0 ) {
		opts->common.time_limit = atoi(val);
	}

	if ( strcmp(key, "indexing.slurm.exclusive") == 0 ) {
		opts->common.exclusive = atoi(val);
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
	char *sc_rel_filename;
	char *stream_rel_filename;
	char *stdout_rel_filename;
	char *stderr_rel_filename;
	char *fg_rel_filename;
	char *intermediate_rel_filename;
	char *harvest_rel_filename;
	char *slurm_prologue;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	stream_rel_filename = relative_to_cwd(workdir, "ambi.stream");
	sc_rel_filename = relative_to_cwd(workdir, "run_ambigator.sh");
	stdout_rel_filename = relative_to_cwd(workdir, "stdout.log");
	stderr_rel_filename = relative_to_cwd(workdir, "stderr.log");
	fg_rel_filename = relative_to_cwd(workdir, "fg.dat");
	intermediate_rel_filename = relative_to_cwd(workdir, "ambigator-input.stream");
	harvest_rel_filename = relative_to_cwd(workdir, "parameters.json");

	slurm_prologue = sbatch_bits(&opts->common, job_title, NULL,
	                             stdout_rel_filename, stderr_rel_filename);

	if ( !write_ambigator_script(sc_rel_filename, input, "$(nproc)",
	                             &proj->ambi_params, stream_rel_filename,
	                             stdout_rel_filename, stderr_rel_filename,
	                             fg_rel_filename,
	                             intermediate_rel_filename,
	                             harvest_rel_filename,
	                             slurm_prologue) )
	{
		job = start_slurm_job(GUI_JOB_AMBIGATOR,
		                      sc_rel_filename, job_title, NULL,
		                      stderr_rel_filename);
		job->niter = proj->ambi_params.niter;
	} else {
		job = NULL;
	}

	free(slurm_prologue);
	g_free(sc_rel_filename);
	g_free(stdout_rel_filename);
	g_free(stderr_rel_filename);
	g_free(harvest_rel_filename);

	if ( job != NULL ) {
		add_indexing_result(proj, job_title, &stream_rel_filename, 1);
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
	char *sc_rel_filename;
	char *output_rel_filename;
	char *stdout_rel_filename;
	char *stderr_rel_filename;
	char *harvest_rel_filename;
	char *log_folder_rel;
	char *slurm_prologue;

	workdir = make_job_folder(job_title, job_notes);
	if ( workdir == NULL ) return NULL;

	sc_rel_filename = relative_to_cwd(workdir, "run_merge.sh");
	output_rel_filename = relative_to_cwd(workdir, "crystfel.hkl");
	stdout_rel_filename = relative_to_cwd(workdir, "stdout.log");
	stderr_rel_filename = relative_to_cwd(workdir, "stderr.log");
	harvest_rel_filename = relative_to_cwd(workdir, "parameters.json");
	log_folder_rel = relative_to_cwd(workdir, "pr-logs");

	slurm_prologue = sbatch_bits(&opts->common, job_title, NULL,
	                             stdout_rel_filename, stderr_rel_filename);

	if ( !write_merge_script(sc_rel_filename, input, "$(nproc)",
	                         &proj->merging_params, output_rel_filename,
	                         stdout_rel_filename, stderr_rel_filename,
	                         harvest_rel_filename,
	                         log_folder_rel,
	                         slurm_prologue) )
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
		job = start_slurm_job(type, sc_rel_filename, job_title, NULL,
		                      stderr_rel_filename);
	} else {
		job = NULL;
	}

	free(slurm_prologue);

	if ( job != NULL ) {

		char *hkl1;
		char *hkl2;

		hkl1 = relative_to_cwd(workdir, "crystfel.hkl1");
		hkl2 = relative_to_cwd(workdir, "crystfel.hkl2");

		add_merge_result(proj, job_title, input->name,
		                 output_rel_filename, hkl1, hkl2);
		g_free(hkl1);
		g_free(hkl2);
	}

	g_object_unref(workdir);
	g_free(sc_rel_filename);
	g_free(stdout_rel_filename);
	g_free(stderr_rel_filename);
	g_free(output_rel_filename);
	g_free(harvest_rel_filename);
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

	if ( strcmp(key, "merging.slurm.reservation") == 0 ) {
		opts->common.reservation = strdup(val);
	}

	if ( strcmp(key, "merging.slurm.qos") == 0 ) {
		opts->common.qos = strdup(val);
	}

	if ( strcmp(key, "merging.slurm.time_limit") == 0 ) {
		opts->common.time_limit = atoi(val);
	}

	if ( strcmp(key, "merging.slurm.exclusive") == 0 ) {
		opts->common.exclusive = atoi(val);
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

	if ( strcmp(key, "ambi.slurm.reservation") == 0 ) {
		opts->common.reservation = strdup(val);
	}

	if ( strcmp(key, "ambi.slurm.qos") == 0 ) {
		opts->common.qos = strdup(val);
	}

	if ( strcmp(key, "ambi.slurm.time_limit") == 0 ) {
		opts->common.time_limit = atoi(val);
	}

	if ( strcmp(key, "ambi.slurm.exclusive") == 0 ) {
		opts->common.exclusive = atoi(val);
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
	be->free_task = free_task;
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
