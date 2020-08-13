/*
 * gui_index.c
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
#include <getopt.h>
#include <string.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms-compat.h>
#include <assert.h>

#include <datatemplate.h>
#include <peaks.h>

#include "crystfel_gui.h"
#include "crystfelimageview.h"
#include "crystfelindexingopts.h"

void cell_explorer_sig(struct crystfelproject *proj)
{
	GSubprocess *sp;
	GError *error = NULL;

	sp = g_subprocess_new(G_SUBPROCESS_FLAGS_NONE, &error,
	                      "cell_explorer", "test.stream", NULL);
	if ( sp == NULL ) {
		ERROR("Failed to start cell_explorer: %s\n",
		      error->message);
		g_error_free(error);
		return;
	}

	STATUS("Starting cell_explorer...\n");
}


static void index_all_response_sig(GtkWidget *dialog, gint resp,
                                   struct crystfelproject *proj)
{
	if ( resp == GTK_RESPONSE_OK ) {
		STATUS("OK!\n");
	}

	gtk_widget_destroy(dialog);
}


static GtkWidget *make_backend_opts(struct crystfelproject *proj)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *combo;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Batch system:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo),
	                   FALSE, FALSE, 0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "local",
	                "Local (run on this computer)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "slurm",
	                "SLURM");

	return box;
}


gint index_all_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *indexing_opts;

	dialog = gtk_dialog_new_with_buttons("Index all frames",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(index_all_response_sig), proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	indexing_opts = crystfel_indexing_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(indexing_opts),
	                   FALSE, FALSE, 8.0);

	gtk_notebook_append_page(GTK_NOTEBOOK(indexing_opts),
	                         make_backend_opts(proj),
	                         gtk_label_new("Cluster/batch system"));

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}


static void index_one_response_sig(GtkWidget *dialog, gint resp,
                                   struct crystfelproject *proj)
{
	if ( resp == GTK_RESPONSE_OK ) {
		STATUS("OK!\n");
	}

	gtk_widget_destroy(dialog);
}


gint index_one_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *indexing_opts;

	dialog = gtk_dialog_new_with_buttons("Index one frame",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(index_one_response_sig), proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	indexing_opts = crystfel_indexing_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(indexing_opts),
	                   FALSE, FALSE, 8.0);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}
