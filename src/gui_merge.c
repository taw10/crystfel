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

#include "gui_project.h"
#include "crystfel_gui.h"
#include "crystfelmergeopts.h"


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


static void merging_response_sig(GtkWidget *dialog, gint resp,
                                 struct new_merging_job_params *njp)
{
	if ( resp == GTK_RESPONSE_OK ) {
		STATUS("Doing it!\n");
	}
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
	GtkWidget *notebook;
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

	notebook = crystfel_merge_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(notebook),
	                   FALSE, FALSE, 8.0);

	job_page = make_merging_job_opts(proj, njp);
	gtk_notebook_prepend_page(GTK_NOTEBOOK(notebook),
	                          job_page,
	                          gtk_label_new("Job name/notes"));

	backend_page = make_merging_backend_opts(njp);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),
	                         backend_page,
	                         gtk_label_new("Cluster/batch system"));

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}
