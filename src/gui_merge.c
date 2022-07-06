/*
 * gui_merge.c
 *
 * Merging via CrystFEL GUI
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
	struct gui_job_notes_page *notes_page;
	GtkWidget *backend_combo;
	GtkWidget *backend_opts_widget;
	GtkWidget *backend_opts_box;
	GtkWidget *job_title_entry;
	GtkWidget *model_combo;
	GtkWidget *input_combo;
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
	opts->pr_logs = crystfel_merge_opts_get_pr_logs(mo);
	opts->twin_sym = safe_strdup(crystfel_merge_opts_get_twin_sym(mo));
	opts->min_res = crystfel_merge_opts_get_min_res(mo);
	opts->push_res = crystfel_merge_opts_get_push_res(mo);
}


static int run_merging(struct crystfelproject *proj,
                       const char *results_name,
                       int backend_idx,
                       const char *job_title,
                       const char *job_notes)
{
	struct crystfel_backend *be;
	void *job_priv;
	struct gui_indexing_result *input;

	/* Which result to merge? */
	input = find_indexing_result_by_name(proj, results_name);
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
		const char *results_name;

		get_merging_opts(&njp->proj->merging_params,
		                 CRYSTFEL_MERGE_OPTS(njp->proj->merging_opts));

		backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(njp->backend_combo));
		if ( backend_idx < 0 ) return;

		job_title = gtk_entry_get_text(GTK_ENTRY(njp->job_title_entry));
		job_notes = get_all_text(GTK_TEXT_VIEW(njp->notes_page->textview));

		if ( job_title[0] == '\0' ) {
			ERROR("You must provide a job name.\n");
			return;
		}

		free(njp->proj->merging_new_job_title);
		njp->proj->merging_new_job_title = strdup(job_title);

		results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(njp->input_combo));
		if ( results_name == NULL ) {
			ERROR("Please select the input\n");
			return;
		}
		if ( run_merging(njp->proj, results_name,
		                 backend_idx, job_title, job_notes) == 0 )
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
	crystfel_merge_opts_set_pr_logs(mo, opts->pr_logs);
	crystfel_merge_opts_set_twin_sym(mo, opts->twin_sym);
	crystfel_merge_opts_set_min_res(mo, opts->min_res);
	crystfel_merge_opts_set_push_res(mo, opts->push_res);
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
	GtkWidget *label;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *backend_page;
	int i;
	char *new_title;
	struct new_merging_job_params *njp;

	if ( proj->merging_opts != NULL ) return FALSE;

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

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Job/output name:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	njp->job_title_entry = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->job_title_entry),
	                   TRUE, TRUE, 2.0);
	gtk_entry_set_placeholder_text(GTK_ENTRY(njp->job_title_entry),
	                               "merge-trial-1");
	new_title = make_new_job_title(proj->merging_new_job_title);
	if ( new_title != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(njp->job_title_entry), new_title);
		free(new_title);
	}

	label = gtk_label_new("Input:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	njp->input_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->input_combo),
	                   FALSE, FALSE, 4.0);
	for ( i=0; i<proj->n_results; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(njp->input_combo),
		                          proj->results[i].name,
		                          proj->results[i].name);
	}
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(njp->input_combo),
	                            selected_result(proj));


	proj->merging_opts = crystfel_merge_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(proj->merging_opts),
	                   FALSE, FALSE, 8.0);

	set_merging_opts(&proj->merging_params, CRYSTFEL_MERGE_OPTS(proj->merging_opts));

	backend_page = make_merging_backend_opts(njp);
	gtk_notebook_append_page(GTK_NOTEBOOK(proj->merging_opts),
	                          backend_page,
	                          gtk_label_new("Cluster/batch system"));

	njp->notes_page = add_job_notes_page(proj->merging_opts);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}


static int write_partialator_script(const char *filename,
                                    struct gui_indexing_result *input,
                                    const char *n_thread_str,
                                    struct merging_params *params,
                                    const char *out_hkl,
                                    const char *stdout_filename,
                                    const char *stderr_filename,
                                    const char *harvest_filename,
                                    const char *log_folder)
{
	FILE *fh;
	int i;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	fprintf(fh, "#!/bin/sh\n");

	fprintf(fh, "partialator \\\n");

	for ( i=0; i<input->n_streams; i++ ) {
		fprintf(fh, "\"%s\" \\\n", input->streams[i]);
	}

	fprintf(fh, " --model=%s", params->model);

	fprintf(fh, " -j %s", n_thread_str);
	fprintf(fh, " -o \"%s\"", out_hkl);
	fprintf(fh, " -y %s", params->symmetry);

	fprintf(fh, " --polarisation=%s", params->polarisation);
	fprintf(fh, " --min-measurements=%i", params->min_measurements);
	fprintf(fh, " --max-adu=%f", params->max_adu);
	fprintf(fh, " --min-res=%f", params->min_res);
	fprintf(fh, " --push-res=%f", params->push_res);

	if ( params->twin_sym != NULL ) {
		fprintf(fh, " -w %s", params->twin_sym);
	}

	if ( params->custom_split != NULL ) {
		fprintf(fh, " --custom-split=\"%s\"", params->custom_split);
	}

	if ( !params->scale ) {
		fprintf(fh, " --no-scale");
	}

	if ( !params->bscale ) {
		fprintf(fh, " --no-Bscale");
	}

	if ( !params->postref ) {
		fprintf(fh, " --no-pr");
	}

	if ( !params->deltacchalf ) {
		fprintf(fh, " --no-deltacchalf");
	}

	if ( !params->pr_logs ) {
		fprintf(fh, " --no-logs");
	}

	fprintf(fh, " --iterations=%i", params->niter);
	fprintf(fh, " --harvest-file=%s", harvest_filename);
	fprintf(fh, " --log-folder=%s", log_folder);

	fprintf(fh, " >%s 2>%s\n", stdout_filename, stderr_filename);

	fclose(fh);
	return 0;
}


static void add_process_hkl(FILE *fh,
                            struct gui_indexing_result *input,
                            struct merging_params *params,
                            const char *out_hkl,
                            const char *stdout_filename,
                            const char *stderr_filename,
                            const char *extra_arg,
                            const char *out_suffix)
{
	int i;

	fprintf(fh, "process_hkl \\\n");

	for ( i=0; i<input->n_streams; i++ ) {
		fprintf(fh, " \"%s\" \\\n", input->streams[i]);
	}

	fprintf(fh, " -o \"%s%s\"", out_hkl, out_suffix);
	fprintf(fh, " -y %s", params->symmetry);

	if ( params->scale ) {
		fprintf(fh, " --scale");
	}

	fprintf(fh, " --polarisation=%s", params->polarisation);
	fprintf(fh, " --min-measurements=%i", params->min_measurements);
	fprintf(fh, " --max-adu=%f", params->max_adu);
	fprintf(fh, " --min-res=%f", params->min_res);
	fprintf(fh, " --push-res=%f", params->push_res);
	fprintf(fh, " %s >>%s 2>>%s\n",
	        extra_arg, stdout_filename, stderr_filename);
}


static int write_process_hkl_script(const char *filename,
                                    struct gui_indexing_result *input,
                                    struct merging_params *params,
                                    const char *out_hkl,
                                    const char *stdout_filename,
                                    const char *stderr_filename)
{
	FILE *fh;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	fprintf(fh, "#!/bin/sh\n");

	add_process_hkl(fh, input, params, out_hkl,
	                stdout_filename, stderr_filename, "", "");
	add_process_hkl(fh, input, params, out_hkl,
	                stdout_filename, stderr_filename, "--even-only", "1");
	add_process_hkl(fh, input, params, out_hkl,
	                stdout_filename, stderr_filename, "--odd-only", "2");

	fclose(fh);
	return 0;
}


int write_merge_script(const char *filename,
                       struct gui_indexing_result *input,
                       const char *n_thread_str,
                       struct merging_params *params,
                       const char *out_hkl,
                       const char *stdout_filename,
                       const char *stderr_filename,
                       const char *harvest_filename,
                       const char *log_folder)
{
	if ( strcmp(params->model, "process_hkl") == 0 ) {
		return write_process_hkl_script(filename, input,
		                                params, out_hkl,
		                                stdout_filename,
		                                stderr_filename);
	} else {
		return write_partialator_script(filename, input, n_thread_str,
		                                params, out_hkl,
		                                stdout_filename,
		                                stderr_filename,
		                                harvest_filename,
		                                log_folder);
	}
}


double read_merge_progress(char *logfile_str, enum gui_job_type type)
{
	FILE *fh;
	double frac_complete = 0.0;

	fh = fopen(logfile_str, "r");
	if ( fh == NULL ) return 0.0;

	do {
		char line[1024];

		if ( fgets(line, 1024, fh) == NULL ) break;

		if ( type == GUI_JOB_PROCESS_HKL ) {
			frac_complete += 1.0/3.0;
		} else if ( type == GUI_JOB_PROCESS_HKL_SCALE ) {
			frac_complete += 1.0/6.0;
		} else {
			int cycle, max_cycles;
			assert(type == GUI_JOB_PARTIALATOR);
			if ( strcmp(line, "Initial partiality calculation...\n") == 0 ) {
				frac_complete = 0.1;
			}
			if ( sscanf(line, "Scaling and refinement cycle %d of %d\n",
			            &cycle, &max_cycles) == 2 )
			{
				frac_complete = 0.1 + 0.8*(double)cycle/max_cycles;
			}
			if ( strcmp(line, "Final merge...\n") == 0 ) {
				frac_complete = 0.9;
			}
			if ( strncmp(line, "Writing two-way split", 20) == 0 ) {
				frac_complete = 1.0;
			}

		}
	} while ( 1 );

	fclose(fh);

	return frac_complete;
}
