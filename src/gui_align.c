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
	GtkWidget *out_of_plane;
	GtkWidget *level;
};


static int run_align(struct align_window *win)
{
	STATUS("Mock detector alignment!\n");
	return 0;
}


static void align_response_sig(GtkWidget *dialog, gint resp,
                               struct align_window *win)
{
	int r = 0;

	if ( resp == GTK_RESPONSE_ACCEPT ) {
		r = run_align(win);
	}

	if ( !r ) gtk_widget_destroy(dialog);
}


gint align_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *content_area;
	struct align_window *win;

	win = malloc(sizeof(struct align_window));
	if ( win == NULL ) return 0;

	win->proj = proj;

	win->window = gtk_dialog_new_with_buttons("Align detector",
	                                          GTK_WINDOW(proj->window),
	                                          GTK_DIALOG_DESTROY_WITH_PARENT,
	                                          "Cancel", GTK_RESPONSE_CANCEL,
	                                          "Run", GTK_RESPONSE_OK,
	                                          NULL);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(win->window));
	gtk_box_pack_start(GTK_BOX(content_area), GTK_WIDGET(vbox), TRUE, TRUE, 0.0);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0.0);
	label = gtk_label_new("Refinement level:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 16.0);
	win->level = gtk_spin_button_new_with_range(0.0, 9.0, 1.0);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->level),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(win->level, "--level");

	win->out_of_plane = gtk_check_button_new_with_label("Refine out-of-plane positions and tilts");
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(win->out_of_plane),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(win->out_of_plane, "--out-of-plane");

	g_signal_connect(G_OBJECT(win->window), "response",
	                 G_CALLBACK(align_response_sig), win);

	gtk_dialog_set_default_response(GTK_DIALOG(win->window),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(win->window);

	return FALSE;
}
