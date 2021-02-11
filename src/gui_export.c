/*
 * gui_export.c
 *
 * Export data from CrystFEL GUI
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

#include <utils.h>
#include <reflist-utils.h>
#include <cell-utils.h>

#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"


struct export_window
{
	struct crystfelproject *proj;
	GtkWidget *cell_chooser;
	GtkWidget *limit_res;
	GtkWidget *min_res;
	GtkWidget *max_res;
	GtkWidget *dataset;
	GtkWidget *format;
};


static void export_response_sig(GtkWidget *dialog, gint resp,
                                struct export_window *f)
{
	if ( resp != GTK_RESPONSE_APPLY ) {
		gtk_widget_destroy(dialog);
		return;
	}
}


gint export_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	char tmp[64];
	struct export_window *win;

	win = malloc(sizeof(struct export_window));
	if ( win == NULL ) return 0;

	win->proj = proj;

	dialog = gtk_file_chooser_dialog_new("Export data",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_FILE_CHOOSER_ACTION_SAVE,
	                                     GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                     GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
	                                     NULL);
	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog),
	                                               TRUE);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(export_response_sig),
	                 win);

	vbox = gtk_vbox_new(FALSE, 0.0);
	gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(dialog),
	                                  GTK_WIDGET(vbox));
	gtk_container_set_border_width(GTK_CONTAINER(vbox), 4);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Results to export:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->dataset = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->dataset),
	                   FALSE, FALSE, 4.0);

	label = gtk_label_new("Format");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->format = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->format),
	                   FALSE, FALSE, 4.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "mtz",
	                          "MTZ, plain");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "mtz-bij",
	                          "MTZ, Bijvoet pairs together");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "xscale",
	                          "XSCALE");

	label = gtk_label_new("Unit cell file:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->cell_chooser = gtk_file_chooser_button_new("Unit cell file",
	                                              GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(win->cell_chooser),
	                                TRUE);
	/* Use the "FoM" cell file because there should only be one
	 * point of truth for the "final" cell parameters.  Eventually, I hope
	 * to determine this automatically. */
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(win->cell_chooser),
	                              proj->fom_cell_filename);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->cell_chooser),
	                   FALSE, FALSE, 4.0);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	win->limit_res = gtk_check_button_new_with_label("Restrict resolution range:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->limit_res),
	                   FALSE, FALSE, 4.0);
	win->min_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->min_res), 4);
	snprintf(tmp, 64, "%.2f", proj->export_res_min);
	gtk_entry_set_text(GTK_ENTRY(win->min_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->min_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->max_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(win->max_res), 4);
	snprintf(tmp, 64, "%.2f", proj->export_res_max);
	gtk_entry_set_text(GTK_ENTRY(win->max_res), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->max_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Å");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(dialog);

	return FALSE;
}
