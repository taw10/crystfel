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

static void unitcell_response_sig(GtkWidget *dialog, gint resp,
                                  struct crystfelproject *proj)
{
	const char *algo;

	algo = gtk_combo_box_get_active_id(GTK_COMBO_BOX(proj->unitcell_combo));

	gtk_widget_destroy(dialog);
	if ( resp != GTK_RESPONSE_OK ) {
		proj->unitcell_combo = NULL;
		return;
	}

	proj->backend->run_unitcell(proj, algo);
	proj->unitcell_combo = NULL;
}


gint unitcell_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *combo;

	if ( proj->unitcell_combo != NULL ) return FALSE;

	dialog = gtk_dialog_new_with_buttons("Determine unit cell",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(unitcell_response_sig), proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);
	proj->peak_vbox = vbox;

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox), FALSE, FALSE, 8.0);
	label = gtk_label_new("Unit cell determination algorithm");
	gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 2.0);
	combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo), TRUE, TRUE, 2.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
	                          "mosflm-nocell-nolatt",
	                          "MOSFLM");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
	                          "dirax",
	                          "DirAx");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
	                          "asdf-nocell",
	                          "ASDF");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
	                          "xds-nocell-nolatt",
	                          "xds");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
	                          "xgandalf-nocell",
	                          "XGANDALF");
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
	proj->unitcell_combo = combo;

	gtk_widget_show_all(dialog);

	return FALSE;
}
