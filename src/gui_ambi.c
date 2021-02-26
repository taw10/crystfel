/*
 * gui_ambi.c
 *
 * Resolve indexing ambiguities via CrystFEL GUI
 *
 * Copyright © 2021 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2021 Thomas White <taw@physics.org>
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

#include "crystfelsymmetryselector.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"


struct ambi_window
{
	struct crystfelproject *proj;
	struct gui_job_notes_page *notes_page;
	GtkWidget *jobname;
	GtkWidget *dataset;
	GtkWidget *limit_res;
	GtkWidget *min_res;
	GtkWidget *max_res;
	GtkWidget *niter;
	GtkWidget *use_ncorr;
	GtkWidget *ncorr;
	GtkWidget *sym;
	GtkWidget *source_sym;
	GtkWidget *operator;
	GtkWidget *use_source_sym;
	GtkWidget *use_operator;
	GtkWidget *backend_combo;
	GtkWidget *backend_opts_widget;
	GtkWidget *backend_opts_box;
};


static int run_ambi(struct crystfelproject *proj,
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
	job_priv = be->run_ambi(job_title, job_notes, proj, input,
	                        be->ambi_opts_priv);

	if ( job_priv != NULL ) {
		char name[256];
		snprintf(name, 255, "Resolving indexing ambiguity (%s)", job_title);
		add_running_task(proj, name, be, job_priv);
		return 0;
	} else {
		return 1;
	}
}



static char *get_sym(GtkWidget *w)
{
	return crystfel_symmetry_selector_get_group_symbol(CRYSTFEL_SYMMETRY_SELECTOR(w));
}


static char *get_str(GtkWidget *w)
{
	const char *text = gtk_entry_get_text(GTK_ENTRY(w));
	if ( text == NULL ) return NULL;
	if ( text[0] == '\0' ) return NULL;
	return strdup(text);
}


static void ambi_response_sig(GtkWidget *dialog, gint resp,
                              struct ambi_window *win)
{
	int r = 0;

	if ( resp == GTK_RESPONSE_ACCEPT ) {

		int backend_idx;
		const char *job_title;
		char *job_notes;
		const char *results_name;

		win->proj->ambi_params.use_res = get_bool(win->limit_res);
		win->proj->ambi_params.res_min = get_float(win->min_res);
		win->proj->ambi_params.res_max = get_float(win->max_res);
		win->proj->ambi_params.niter = get_uint(win->niter);
		win->proj->ambi_params.use_ncorr = get_bool(win->use_ncorr);
		win->proj->ambi_params.ncorr = get_uint(win->ncorr);
		win->proj->ambi_params.sym = get_sym(win->sym);
		win->proj->ambi_params.source_sym = get_sym(win->source_sym);
		win->proj->ambi_params.operator = get_str(win->operator);
		win->proj->ambi_params.use_operator = get_bool(win->use_operator);

		/* "Minimum resolution" should be the bigger number */
		if ( win->proj->ambi_params.res_min < win->proj->ambi_params.res_max ) {
			double tmp = win->proj->ambi_params.res_min;
			win->proj->ambi_params.res_min = win->proj->ambi_params.res_max;
			win->proj->ambi_params.res_max = tmp;
		}

		backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(win->backend_combo));
		if ( backend_idx < 0 ) return;

		job_title = gtk_entry_get_text(GTK_ENTRY(win->jobname));
		job_notes = get_all_text(GTK_TEXT_VIEW(win->notes_page->textview));

		if ( job_title[0] == '\0' ) {
			ERROR("You must provide a job name.\n");
			return;
		}

		free(win->proj->ambi_new_job_title);
		win->proj->ambi_new_job_title = strdup(job_title);

		results_name = gtk_combo_box_get_active_id(GTK_COMBO_BOX(win->dataset));
		if ( results_name == NULL ) {
			ERROR("Please select the input\n");
			return;
		}
		if ( run_ambi(win->proj, results_name,
		              backend_idx, job_title, job_notes) == 0 )
		{
			gtk_widget_destroy(dialog);
			win->proj->ambi_opts = NULL;
		}

		free(job_notes);

	}

	if ( !r ) {
		gtk_widget_destroy(dialog);
		win->proj->ambi_opts = NULL;
	}
}


static GtkWidget *make_ambigator_options(struct ambi_window *win)
{
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	char tmp[64];
	struct crystfelproject *proj = win->proj;

	vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), 4);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Target ('real') point group:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->sym = crystfel_symmetry_selector_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->sym),
	                   FALSE, FALSE, 4.0);
	crystfel_symmetry_selector_set_group_symbol(CRYSTFEL_SYMMETRY_SELECTOR(win->sym),
	                                            proj->ambi_params.sym);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->use_source_sym = gtk_radio_button_new_with_label(NULL, "Source ('twinned') point group:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->use_source_sym),
	                   FALSE, FALSE, 4.0);
	win->source_sym = crystfel_symmetry_selector_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->source_sym),
	                   FALSE, FALSE, 4.0);
	crystfel_symmetry_selector_set_group_symbol(CRYSTFEL_SYMMETRY_SELECTOR(win->source_sym),
	                                            proj->ambi_params.source_sym);
	g_signal_connect(G_OBJECT(win->use_source_sym), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->source_sym);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->use_operator = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(win->use_source_sym),
	                                                                "Ambiguity operation:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->use_operator),
	                   FALSE, FALSE, 4.0);
	win->operator = gtk_entry_new();
	if ( proj->ambi_params.operator != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(win->operator), proj->ambi_params.operator);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->operator),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(win->use_operator), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->operator);

	if ( proj->ambi_params.use_operator ) {
		gtk_widget_set_sensitive(win->source_sym, FALSE);
		set_active(win->use_operator, TRUE);
	} else {
		gtk_widget_set_sensitive(win->operator, FALSE);
		set_active(win->use_source_sym, TRUE);
	}

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->limit_res = gtk_check_button_new_with_label("Restrict resolution range:");
	set_active(win->limit_res, proj->ambi_params.use_res);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->limit_res),
	                   FALSE, FALSE, 4.0);
	win->min_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->min_res), 6);
	snprintf(tmp, 64, "%.2f", proj->ambi_params.res_min);
	gtk_entry_set_text(GTK_ENTRY(win->min_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->min_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->max_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->max_res), 6);
	snprintf(tmp, 64, "%.2f", proj->ambi_params.res_max);
	gtk_entry_set_text(GTK_ENTRY(win->max_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->max_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Å");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(win->limit_res), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->min_res);
	g_signal_connect(G_OBJECT(win->limit_res), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->max_res);
	gtk_widget_set_sensitive(win->min_res, proj->ambi_params.use_res);
	gtk_widget_set_sensitive(win->max_res, proj->ambi_params.use_res);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->use_ncorr = gtk_check_button_new_with_label("Limit number of correlations to:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->use_ncorr),
	                   FALSE, FALSE, 4.0);
	win->ncorr = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->ncorr), 6);
	snprintf(tmp, 64, "%i", proj->ambi_params.ncorr);
	gtk_entry_set_text(GTK_ENTRY(win->ncorr), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->ncorr),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(win->use_ncorr), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->ncorr);
	set_active(win->use_ncorr, proj->ambi_params.use_ncorr);
	gtk_widget_set_sensitive(win->ncorr, proj->ambi_params.use_ncorr);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Number of iterations:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->niter = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->niter), 2);
	snprintf(tmp, 64, "%i", proj->ambi_params.niter);
	gtk_entry_set_text(GTK_ENTRY(win->niter), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->niter),
	                   FALSE, FALSE, 4.0);

	gtk_widget_show_all(vbox);
	return vbox;
}


static void ambi_backend_changed_sig(GtkWidget *combo,
                                     struct ambi_window *win)
{
	int backend_idx;
	struct crystfel_backend *be;

	backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	if ( backend_idx < 0 ) return;
	win->proj->ambi_backend_selected = backend_idx;

	be = &win->proj->backends[backend_idx];

	if ( win->backend_opts_widget != NULL ) {
		gtk_widget_destroy(win->backend_opts_widget);
	}

	win->backend_opts_widget = be->make_ambi_parameters_widget(be->ambi_opts_priv);

	gtk_box_pack_start(GTK_BOX(win->backend_opts_box),
	                   GTK_WIDGET(win->backend_opts_widget),
	                   FALSE, FALSE, 0);
	gtk_widget_show_all(win->backend_opts_widget);
}


static GtkWidget *make_ambi_backend_opts(struct ambi_window *win)
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

	win->backend_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->backend_combo),
	                   FALSE, FALSE, 0);

	for ( i=0; i<win->proj->n_backends; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->backend_combo),
		                          win->proj->backends[i].name,
		                          win->proj->backends[i].friendly_name);
	}

	win->backend_opts_box = gtk_box_new(GTK_ORIENTATION_VERTICAL,
	                                    0);
	gtk_box_pack_start(GTK_BOX(box),
	                   GTK_WIDGET(win->backend_opts_box),
	                   FALSE, FALSE, 0);
	win->backend_opts_widget = NULL;

	/* win->backend_opts{_box} must exist before the following */
	g_signal_connect(G_OBJECT(win->backend_combo), "changed",
	                 G_CALLBACK(ambi_backend_changed_sig), win);
	gtk_combo_box_set_active(GTK_COMBO_BOX(win->backend_combo),
	                         win->proj->ambi_backend_selected);

	gtk_widget_show_all(box);

	return box;
}


gint ambi_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *notebook;
	GtkWidget *backend_page;
	struct ambi_window *win;
	char *new_title;
	int i;

	if ( proj->ambi_opts != NULL ) return FALSE;

	win = malloc(sizeof(struct ambi_window));
	if ( win == NULL ) return 0;

	win->proj = proj;

	dialog = gtk_dialog_new_with_buttons("Resolve indexing ambiguity",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Close", GTK_RESPONSE_CLOSE,
	                                     "Run", GTK_RESPONSE_ACCEPT,
	                                     NULL);
	proj->ambi_opts = dialog;

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(ambi_response_sig),
	                 win);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), 4);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Job/output name:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->jobname = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->jobname), 16);
	gtk_entry_set_placeholder_text(GTK_ENTRY(win->jobname),
	                               "ambi-trial-1");
	new_title = make_new_job_title(proj->ambi_new_job_title);
	if ( new_title != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(win->jobname), new_title);
		free(new_title);
	}
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->jobname),
	                   TRUE, TRUE, 4.0);

	label = gtk_label_new("Input:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->dataset = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->dataset),
	                   FALSE, FALSE, 4.0);
	for ( i=0; i<proj->n_results; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->dataset),
		                          proj->results[i].name,
		                          proj->results[i].name);
	}

	notebook = gtk_notebook_new();
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_LEFT);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(notebook),
	                   FALSE, FALSE, 4.0);

	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
	                         make_ambigator_options(win),
	                         gtk_label_new("Indexing ambiguity"));

	backend_page = make_ambi_backend_opts(win);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
	                          backend_page,
	                          gtk_label_new("Cluster/batch system"));

	win->notes_page = add_job_notes_page(notebook);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(dialog);

	return FALSE;
}


int write_ambigator_script(const char *filename,
                           struct gui_indexing_result *input,
                           const char *n_thread_str,
                           struct ambi_params *params,
                           const char *out_stream)
{
	FILE *fh;
	char *exe_path;
	int i;

	fh = fopen(filename, "w");
	if ( fh == NULL ) return 1;

	fprintf(fh, "#!/bin/sh\n");

	fprintf(fh, "cat \\\n");
	for ( i=0; i<input->n_streams; i++ ) {
		fprintf(fh, "%s \\\n", input->streams[i]);
	}
	fprintf(fh, " > ambigator-input.stream\n");

	exe_path = get_crystfel_exe("ambigator");
	if ( exe_path == NULL ) return 1;
	fprintf(fh, "%s ambigator-input.stream \\\n", exe_path);

	fprintf(fh, " -j %s", n_thread_str);
	fprintf(fh, " -o %s", out_stream);
	fprintf(fh, " -y %s", params->sym);
	if ( params->use_operator ) {
		fprintf(fh, " --operator=%s", params->operator);
	} else {
		fprintf(fh, " -w %s", params->source_sym);
	}

	if ( params->use_res ) {
		fprintf(fh, " --lowres=%f", params->res_min);
		fprintf(fh, " --highres=%f", params->res_max);
	}

	if ( params->use_ncorr ) {
		fprintf(fh, " --ncorr=%i", params->ncorr);
	}

	fprintf(fh, " --iterations=%i", params->niter);
	fprintf(fh, " --fg-graph=fg.dat");
	fprintf(fh, " >stdout.log 2>stderr.log\n");

	fclose(fh);
	return 0;
}


double read_ambigator_progress(char *logfile_str, int niter)
{
	FILE *fh;
	double iter_inc;
	double frac_complete = 0.0;

	iter_inc = 0.8/niter;

	fh = fopen(logfile_str, "r");
	if ( fh == NULL ) return 0.0;

	do {
		char line[1024];
		int junk;

		if ( fgets(line, 1024, fh) == NULL ) break;

		if ( strncmp(line, "Mean number of correlations per crystal:", 40) == 0 ) {
			frac_complete = 0.1;
		}
		if ( strncmp(line, "Mean f,g =", 10) == 0 ) {
			frac_complete += iter_inc;
		}
		if ( sscanf(line, "%d assignments are different from "
		                  "their starting values\n", &junk) == 1 )
		{
			frac_complete = 1.0;
		}

	} while ( 1 );

	fclose(fh);

	printf("got %f\n", frac_complete);
	return frac_complete;
}
