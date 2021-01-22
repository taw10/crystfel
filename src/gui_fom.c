/*
 * gui_fom.c
 *
 * Figures of merit via CrystFEL GUI
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
#include "gtk-util-routines.h"


static void fom_response_sig(GtkWidget *dialog, gint resp,
                             struct crystfelproject *proj)
{
	gtk_widget_destroy(dialog);
}


static void add_fom(GtkWidget *menu,
                    const char *text,
                    const char *markup)
{
	GtkWidget *label;
	GtkWidget *item;

	item = gtk_check_menu_item_new();

	label = gtk_label_new(text);
	if ( markup != NULL ) {
		gtk_label_set_markup(GTK_LABEL(label), markup);
	}
	gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);

	gtk_container_add(GTK_CONTAINER(item), label);
	gtk_widget_show_all(item);

	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
}


static void add_separator(GtkWidget *menu)
{
	GtkWidget *item;
	item = gtk_separator_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	gtk_widget_show_all(item);
}


static GtkWidget *make_fom_menu()
{
	GtkWidget *menu;
	GtkWidget *item;

	menu = gtk_menu_new();

	item = gtk_tearoff_menu_item_new();
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	gtk_widget_show_all(item);

	add_fom(menu, "I/σ(I), completeness etc", NULL);
	add_fom(menu, "Rsplit", "R<sub>split</sub>");
	add_fom(menu, "CC", "CC<sub>½</sub>");
	add_fom(menu, "CC*", "CC<sup>*</sup>");
	add_separator(menu);
	add_fom(menu, "CCano", "CC<sub>ano</sub>");
	add_fom(menu, "Rano", "R<sub>ano</sub>");
	add_fom(menu, "Rano ÷ Rsplit", "R<sub>ano</sub> ÷ R<sub>split</sub>");
	add_fom(menu, "RMS anomalous correlation ratio", NULL);
	add_separator(menu);
	add_fom(menu, "Fraction of differences within 1σ(I)", NULL);
	add_fom(menu, "Fraction of differences within 2σ(I)", NULL);

	return menu;
}


gint fom_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *combo;
	GtkWidget *button;
	GtkWidget *entry;
	GtkWidget *check;
	GtkWidget *da;
	int i;
	char tmp[64];

	dialog = gtk_dialog_new_with_buttons("Calculate figures of merit",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Close", GTK_RESPONSE_CLOSE,
	                                     NULL);

	g_signal_connect(G_OBJECT(dialog), "response",
	                 G_CALLBACK(fom_response_sig),
	                 proj);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Merged results:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo),
	                   FALSE, FALSE, 4.0);
	for ( i=0; i<proj->n_merge_results; i++) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo),
		                          proj->merge_results[i].name,
		                          proj->merge_results[i].name);
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(combo),
	                         proj->n_merge_results-1);

	label = gtk_label_new("Figures of merit to show:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	button = gtk_menu_button_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);
	gtk_menu_button_set_popup(GTK_MENU_BUTTON(button),
	                          make_fom_menu());

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Resolution range:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
	snprintf(tmp, 64, "%.2f", proj->fom_res_min);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
	snprintf(tmp, 64, "%.2f", proj->fom_res_max);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Å.  Number of bins:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 4.0);
	snprintf(tmp, 64, "%i", proj->fom_nbins);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);

	hbox = gtk_hbox_new(FALSE, 0.0);
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	check = gtk_check_button_new_with_label("Discard reflections with I/sigI less than");
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check),
	                             proj->fom_use_min_snr);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(check),
	                   FALSE, FALSE, 4.0);
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 4.0);
	snprintf(tmp, 64, "%.2f", proj->fom_min_snr);
	gtk_entry_set_text(GTK_ENTRY(entry), tmp);

	da = gtk_drawing_area_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(da),
	                   FALSE, FALSE, 4.0);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(dialog);

	return FALSE;
}
