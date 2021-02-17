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
};


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
		win->proj->ambi_use_res = get_bool(win->limit_res);
		win->proj->ambi_res_min = get_float(win->min_res);
		win->proj->ambi_res_max = get_float(win->max_res);
		win->proj->ambi_niter = get_uint(win->niter);
		win->proj->ambi_use_ncorr = get_bool(win->use_ncorr);
		win->proj->ambi_ncorr = get_uint(win->ncorr);
		win->proj->ambi_sym = get_sym(win->sym);
		win->proj->ambi_source_sym = get_sym(win->source_sym);
		win->proj->ambi_operator = get_str(win->operator);
	}

	if ( !r ) gtk_widget_destroy(dialog);
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
	                                            proj->ambi_sym);

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
	                                            proj->ambi_source_sym);
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
	gtk_entry_set_text(GTK_ENTRY(win->operator), proj->ambi_operator);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->operator),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(win->use_operator), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->operator);

	if ( proj->ambi_source_sym != NULL ) {
		gtk_widget_set_sensitive(win->operator, FALSE);
		set_active(win->use_source_sym, TRUE);
	} else {
		gtk_widget_set_sensitive(win->source_sym, FALSE);
		set_active(win->use_operator, TRUE);
	}

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->limit_res = gtk_check_button_new_with_label("Restrict resolution range:");
	set_active(win->limit_res, proj->ambi_use_res);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->limit_res),
	                   FALSE, FALSE, 4.0);
	win->min_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->min_res), 6);
	snprintf(tmp, 64, "%.2f", proj->ambi_res_min);
	gtk_entry_set_text(GTK_ENTRY(win->min_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->min_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->max_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->max_res), 6);
	snprintf(tmp, 64, "%.2f", proj->ambi_res_max);
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
	gtk_widget_set_sensitive(win->min_res, proj->ambi_use_res);
	gtk_widget_set_sensitive(win->max_res, proj->ambi_use_res);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->use_ncorr = gtk_check_button_new_with_label("Limit number of correlations to:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->use_ncorr),
	                   FALSE, FALSE, 4.0);
	win->ncorr = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->ncorr), 6);
	snprintf(tmp, 64, "%i", proj->ambi_ncorr);
	gtk_entry_set_text(GTK_ENTRY(win->ncorr), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->ncorr),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(win->use_ncorr), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->ncorr);
	set_active(win->use_ncorr, proj->ambi_use_ncorr);
	gtk_widget_set_sensitive(win->ncorr, proj->ambi_use_ncorr);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Number of iterations:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->niter = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->niter), 2);
	snprintf(tmp, 64, "%i", proj->ambi_niter);
	gtk_entry_set_text(GTK_ENTRY(win->niter), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->niter),
	                   FALSE, FALSE, 4.0);

	return vbox;
}


static GtkWidget *make_ambi_backend_opts(struct ambi_window *win)
{
	return gtk_vbox_new(FALSE, 0.0);
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
	int i;

	win = malloc(sizeof(struct ambi_window));
	if ( win == NULL ) return 0;

	win->proj = proj;

	dialog = gtk_dialog_new_with_buttons("Resolve indexing ambiguity",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Close", GTK_RESPONSE_CLOSE,
	                                     "Run", GTK_RESPONSE_ACCEPT,
	                                     NULL);

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
