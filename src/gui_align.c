/*
 * gui_align.c
 *
 * Align detector via CrystFEL GUI
 *
 * Copyright © 2024 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2024 Thomas White <taw@physics.org>
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

#include "version.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"


struct align_window
{
	struct crystfelproject *proj;
	GtkWidget *window;
	GtkWidget *input_combo;
	GtkWidget *out_of_plane;
	GtkWidget *level;
};


static int run_align(const char *input_name, int level, int out_of_plane,
                     const char *out_geom)
{
	STATUS("Mock detector alignment from run '%s', level %i, %s  -----> %s\n",
	       input_name, level, out_of_plane ? "out of plane" : "in plane",
	       out_geom);

	return 0;
}


static void align_response_sig(GtkWidget *dialog, gint resp,
                               struct align_window *win)
{
	int r = 0;

	if ( resp == GTK_RESPONSE_ACCEPT ) {

		int level;
		const char *input_name;
		int out_of_plane;
		gchar *filename;

		level = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(win->level));

		input_name = get_combo_id(win->input_combo);
		if ( input_name == NULL ) {
			ERROR("Please select the input\n");
			r = 1;
		}

		out_of_plane = get_bool(win->out_of_plane);
		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

		r = run_align(input_name, level, out_of_plane, filename);

		g_free(filename);
	}

	if ( !r ) gtk_widget_destroy(dialog);
}


gint align_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *hbox;
	GtkWidget *label;
	struct align_window *win;
	int i;

	win = malloc(sizeof(struct align_window));
	if ( win == NULL ) return 0;

	win->proj = proj;

	win->window = gtk_file_chooser_dialog_new("Align detector",
	                                          GTK_WINDOW(proj->window),
	                                          GTK_FILE_CHOOSER_ACTION_SAVE,
	                                          GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
	                                          GTK_STOCK_SAVE, GTK_RESPONSE_ACCEPT,
	                                          NULL);
	gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(win->window),
	                                               TRUE);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_file_chooser_set_extra_widget(GTK_FILE_CHOOSER(win->window),
	                                  GTK_WIDGET(hbox));
	gtk_container_set_border_width(GTK_CONTAINER(hbox), 4);

	label = gtk_label_new("Refine using indexing result:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->input_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->input_combo),
	                   FALSE, FALSE, 4.0);
	for ( i=0; i<proj->n_results; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->input_combo),
		                          proj->results[i].name,
		                          proj->results[i].name);
	}
	gtk_combo_box_set_active_id(GTK_COMBO_BOX(win->input_combo),
	                            selected_result(proj));

	label = gtk_label_new("Hierarchy level:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 16.0);
	win->level = gtk_spin_button_new_with_range(0.0, 9.0, 1.0);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->level),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(win->level, "--level");

	win->out_of_plane = gtk_check_button_new_with_label("Include out-of-plane positions and tilts");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->out_of_plane),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(win->out_of_plane, "--out-of-plane");

	g_signal_connect(G_OBJECT(win->window), "response",
	                 G_CALLBACK(align_response_sig), win);

	gtk_dialog_set_default_response(GTK_DIALOG(win->window),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(win->window);

	return FALSE;
}
