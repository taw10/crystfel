/*
 * gui_merge.c
 *
 * Peak search parts of GUI
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <gtk/gtk.h>
#include <assert.h>

#include <utils.h>

#include "gui_project.h"
#include "crystfel_gui.h"
#include "crystfelmergeopts.h"
#include "gtk-util-routines.h"


struct new_merging_job_params {
	struct crystfelproject *proj;
	GtkWidget *backend_combo;
	GtkWidget *backend_opts_widget;
	GtkWidget *backend_opts_box;
	GtkWidget *job_title_entry;
	GtkWidget *job_notes_text;
	GtkWidget *model_combo;
};


static void free_new_merging_job_params(gpointer njp, GClosure *closure)
{
	free(njp);
}


static void get_merging_opts(struct merging_params *opts,
                             CrystFELMergeOpts *mo)
{
	free(opts->model);
	opts->model = strdup(crystfel_merge_opts_get_model(mo));
	free(opts->symmetry);
	opts->symmetry = strdup(crystfel_merge_opts_get_symmetry(mo));
	opts->scale = crystfel_merge_opts_get_scale(mo);
	opts->bscale = crystfel_merge_opts_get_bscale(mo);
	opts->postref = crystfel_merge_opts_get_postref(mo);
	opts->niter = crystfel_merge_opts_get_niter(mo);
	free(opts->polarisation);
	opts->polarisation = strdup(crystfel_merge_opts_get_polarisation(mo));
	opts->deltacchalf = crystfel_merge_opts_get_deltacchalf(mo);
	opts->min_measurements = crystfel_merge_opts_get_min_measurements(mo);
	opts->max_adu = crystfel_merge_opts_get_max_adu(mo);
	free(opts->custom_split);
	opts->custom_split = safe_strdup(crystfel_merge_opts_get_custom_split(mo));
	free(opts->twin_sym);
	opts->twin_sym = safe_strdup(crystfel_merge_opts_get_twin_sym(mo));
	opts->min_res = crystfel_merge_opts_get_min_res(mo);
	opts->push_res = crystfel_merge_opts_get_push_res(mo);
}


static int run_merging(struct crystfelproject *proj,
                       int backend_idx,
                       const char *job_title,
                       const char *job_notes)
{
	struct crystfel_backend *be;
	void *job_priv;
	const gchar *results_name;
	struct gui_result *input;

	/* Which result to merge? */
	results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->results_combo));
	input = find_result_by_name(proj, results_name);
	if ( input == NULL ) {
		ERROR("Please select a result first\n");
		return 1;
	}

	be = &proj->backends[backend_idx];
	job_priv = be->run_merging(job_title, job_notes, proj, input,
	                           be->merging_opts_priv);

	if ( job_priv != NULL ) {
		char name[256];
		snprintf(name, 255, "Merging data (%s)", job_title);
		add_running_task(proj, name, be, job_priv);
		return 0;
	} else {
		return 1;
	}
}


static void merging_response_sig(GtkWidget *dialog, gint resp,
                                 struct new_merging_job_params *njp)
{
	if ( resp == GTK_RESPONSE_OK ) {

		int backend_idx;
		const char *job_title;
		char *job_notes;

		get_merging_opts(&njp->proj->merging_params,
		                 CRYSTFEL_MERGE_OPTS(njp->proj->merging_opts));

		backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(njp->backend_combo));
		if ( backend_idx < 0 ) return;

		job_title = gtk_entry_get_text(GTK_ENTRY(njp->job_title_entry));
		job_notes = get_all_text(GTK_TEXT_VIEW(njp->job_notes_text));

		if ( job_title[0] == '\0' ) {
			ERROR("You must provide a job name.\n");
			return;
		}

		free(njp->proj->merging_new_job_title);
		njp->proj->merging_new_job_title = strdup(job_title);

		if ( run_merging(njp->proj, backend_idx,
		                 job_title, job_notes) == 0 )
		{
			gtk_widget_destroy(dialog);
			njp->proj->merging_opts = NULL;
		}

		free(job_notes);

	} else {
		gtk_widget_destroy(dialog);
		njp->proj->merging_opts = NULL;
	}
}


static void set_merging_opts(struct merging_params *opts,
                             CrystFELMergeOpts *mo)
{
	crystfel_merge_opts_set_model(mo, opts->model);
	crystfel_merge_opts_set_symmetry(mo, opts->symmetry);
	crystfel_merge_opts_set_scale(mo, opts->scale);
	crystfel_merge_opts_set_bscale(mo, opts->bscale);
	crystfel_merge_opts_set_postref(mo, opts->postref);
	crystfel_merge_opts_set_niter(mo, opts->niter);
	crystfel_merge_opts_set_polarisation(mo, opts->polarisation);
	crystfel_merge_opts_set_deltacchalf(mo, opts->deltacchalf);
	crystfel_merge_opts_set_min_measurements(mo, opts->min_measurements);
	crystfel_merge_opts_set_max_adu(mo, opts->max_adu);
	crystfel_merge_opts_set_custom_split(mo, opts->custom_split);
	crystfel_merge_opts_set_twin_sym(mo, opts->twin_sym);
	crystfel_merge_opts_set_min_res(mo, opts->min_res);
	crystfel_merge_opts_set_push_res(mo, opts->push_res);
}


static GtkWidget *make_merging_job_opts(struct crystfelproject *proj,
                                        struct new_merging_job_params *njp)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *scroll;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Job name:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	njp->job_title_entry = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->job_title_entry),
	                   TRUE, TRUE, 2.0);
	if ( proj->merging_new_job_title != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(njp->job_title_entry),
		                   proj->merging_new_job_title);
	}

	label = gtk_label_new("This name will be used for a working subfolder");
	gtk_label_set_markup(GTK_LABEL(label),
	                     "<i>This name will be used for a working subfolder</i>");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	gtk_entry_set_placeholder_text(GTK_ENTRY(njp->job_title_entry),
	                               "merge-trial-1");

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   TRUE, TRUE, 0);
	label = gtk_label_new("Notes:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	njp->job_notes_text = gtk_text_view_new();
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_add(GTK_CONTAINER(scroll), njp->job_notes_text);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroll),
	                                    GTK_SHADOW_ETCHED_IN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(scroll),
	                   TRUE, TRUE, 2.0);

	label = gtk_label_new("The notes above will be placed in the job's folder as 'notes.txt'");
	gtk_label_set_markup(GTK_LABEL(label),
	                     "<i>The notes above will be placed in the job's folder as 'notes.txt'</i>");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	return box;
}


static void merging_backend_changed_sig(GtkWidget *combo,
                                        struct new_merging_job_params *njp)
{
	int backend_idx;
	struct crystfel_backend *be;

	backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	if ( backend_idx < 0 ) return;
	njp->proj->merging_backend_selected = backend_idx;

	be = &njp->proj->backends[backend_idx];

	if ( njp->backend_opts_widget != NULL ) {
		gtk_widget_destroy(njp->backend_opts_widget);
	}

	njp->backend_opts_widget = be->make_merging_parameters_widget(be->merging_opts_priv);

	gtk_box_pack_start(GTK_BOX(njp->backend_opts_box),
	                   GTK_WIDGET(njp->backend_opts_widget),
	                   FALSE, FALSE, 0);
	gtk_widget_show_all(njp->backend_opts_widget);
}


static GtkWidget *make_merging_backend_opts(struct new_merging_job_params *njp)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	int i;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Batch system:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	njp->backend_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->backend_combo),
	                   FALSE, FALSE, 0);

	for ( i=0; i<njp->proj->n_backends; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(njp->backend_combo),
		                          njp->proj->backends[i].name,
		                          njp->proj->backends[i].friendly_name);
	}

	njp->backend_opts_box = gtk_box_new(GTK_ORIENTATION_VERTICAL,
	                                    0);
	gtk_box_pack_start(GTK_BOX(box),
	                   GTK_WIDGET(njp->backend_opts_box),
	                   FALSE, FALSE, 0);
	njp->backend_opts_widget = NULL;

	/* njp->backend_opts{_box} must exist before the following */
	g_signal_connect(G_OBJECT(njp->backend_combo), "changed",
	                 G_CALLBACK(merging_backend_changed_sig), njp);
	gtk_combo_box_set_active(GTK_COMBO_BOX(njp->backend_combo),
	                         njp->proj->merging_backend_selected);

	return box;
}


gint merge_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *backend_page;
	GtkWidget *job_page;
	struct new_merging_job_params *njp;

	njp = malloc(sizeof(struct new_merging_job_params));
	if ( njp == NULL ) return FALSE;

	njp->proj = proj;

	dialog = gtk_dialog_new_with_buttons("Merge",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect_data(G_OBJECT(dialog), "response",
	                      G_CALLBACK(merging_response_sig),
	                      njp, free_new_merging_job_params, 0);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	proj->merging_opts = crystfel_merge_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(proj->merging_opts),
	                   FALSE, FALSE, 8.0);

	set_merging_opts(&proj->merging_params, CRYSTFEL_MERGE_OPTS(proj->merging_opts));

	job_page = make_merging_job_opts(proj, njp);
	gtk_notebook_prepend_page(GTK_NOTEBOOK(proj->merging_opts),
	                          job_page,
	                          gtk_label_new("Job name/notes"));

	backend_page = make_merging_backend_opts(njp);
	gtk_notebook_append_page(GTK_NOTEBOOK(proj->merging_opts),
	                         backend_page,
	                         gtk_label_new("Cluster/batch system"));

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}


static GSList *append_arg_str(GSList *args,
                              const char *label,
                              const char *val)
{
	size_t len;
	char *str;

	len = strlen(label)+strlen(val)+4;
	str = malloc(len);
	if ( str == NULL ) return args;
	snprintf(str, 63, "--%s=%s", label, val);

	return g_slist_append(args, str);
}


static GSList *append_arg_int(GSList *args,
                              const char *label,
                              int val)
{
	char *str = malloc(64);
	if ( str == NULL ) return args;
	snprintf(str, 63, "--%s=%i", label, val);
	return g_slist_append(args, str);
}


static GSList *append_arg_float(GSList *args,
                                const char *label,
                                float val)
{
	char *str = malloc(64);
	if ( str == NULL ) return args;
	snprintf(str, 63, "--%s=%f", label, val);
	return g_slist_append(args, str);
}


static GSList *process_hkl_command_line(struct gui_result *input,
                                        struct merging_params *params)
{
	GSList *args = NULL;
	char *exe_path;

	exe_path = get_crystfel_exe("process_hkl");
	if ( exe_path == NULL ) return NULL;
	args = g_slist_append(args, exe_path);

	/* FIXME: For each stream */
	args = append_arg_str(args, "input", input->streams[0]);

	args = append_arg_str(args, "symmetry", params->symmetry);

	if ( params->scale ) {
		args = g_slist_append(args, strdup("--scale"));
	}

	args = append_arg_str(args, "polarisation",
	                      params->polarisation);

	args = append_arg_int(args, "min-measurements",
	                      params->min_measurements);

	args = append_arg_float(args, "max-adu", params->max_adu);

	args = append_arg_float(args, "min-res", params->min_res);

	args = append_arg_float(args, "push-res", params->push_res);

	return args;
}


static GSList *partialator_command_line(const char *n_thread_str,
                                        struct gui_result *input,
                                        struct merging_params *params)
{
	GSList *args = NULL;
	char *exe_path;

	exe_path = get_crystfel_exe("partialator");
	if ( exe_path == NULL ) return NULL;
	args = g_slist_append(args, exe_path);

	/* FIXME: For each stream */
	args = append_arg_str(args, "input", input->streams[0]);

	args = append_arg_str(args, "symmetry", params->symmetry);

	if ( params->twin_sym != NULL ) {
		args = g_slist_append(args, "-w");
		args = g_slist_append(args, strdup(params->twin_sym));
	}

	if ( params->custom_split != NULL ) {
		args = append_arg_str(args, "custom-split",
		                      params->custom_split);
	}

	if ( !params->scale ) {
		args = g_slist_append(args, "--no-scale");
	}

	if ( !params->bscale ) {
		args = g_slist_append(args, "--no-bscale");
	}

	if ( !params->postref ) {
		args = g_slist_append(args, "--no-pr");
	}

	if ( !params->deltacchalf ) {
		args = g_slist_append(args, "no-deltacchalf");
	}

	args = append_arg_int(args, "iterations", params->niter);

	args = append_arg_str(args, "polarisation",
	                      params->polarisation);

	args = append_arg_int(args, "min-measurements",
	                      params->min_measurements);

	args = append_arg_float(args, "max-adu", params->max_adu);

	args = append_arg_float(args, "min-res", params->min_res);

	args = append_arg_float(args, "push-res", params->push_res);

	return args;
}


char **merging_command_line(const char *n_thread_str,
                            struct gui_result *input,
                            struct merging_params *params)
{
	GSList *args;
	char **arg_strings;
	GSList *args2;
	int i, n;

	if ( strcmp(params->model, "process_hkl") == 0 ) {
		args = process_hkl_command_line(input, params);
	} else {
		args = partialator_command_line(n_thread_str,
		                                input,
		                                params);
	}

	if ( args == NULL ) return NULL;

	n = g_slist_length(args);
	arg_strings = malloc((n+1)*sizeof(char *));
	if ( arg_strings == NULL ) return NULL;

	args2 = args;
	for ( i=0; i<n; i++ ) {
		arg_strings[i] = args2->data;
		args2 = args2->next;
	}
	arg_strings[n] = NULL;
	g_slist_free(args);

	return arg_strings;
}
