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
#include <symmetry.h>

#include "version.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gtk-util-routines.h"


struct export_window
{
	struct crystfelproject *proj;
	GtkWidget *window;
	GtkWidget *cell_chooser;
	GtkWidget *limit_res;
	GtkWidget *min_res;
	GtkWidget *max_res;
	GtkWidget *dataset;
	GtkWidget *format;
};


static int export_to_xds(struct gui_merge_result *result,
                         const char *filename, UnitCell *cell,
                         double min_res, double max_res)
{
	RefList *reflist;
	char *sym_str;
	SymOpList *sym;
	int r;

	reflist = read_reflections_2(result->hkl, &sym_str);
	if ( reflist == NULL ) return 1;
	if ( sym_str == NULL ) return 1;

	sym = get_pointgroup(sym_str);
	if ( sym == NULL ) return 1;

	r = write_to_xds(reflist, sym, cell, min_res, max_res, filename);

	free_symoplist(sym);
	free(sym_str);
	reflist_free(reflist);

	return r;
}


static int export_to_mtz(struct gui_merge_result *result,
                         const char *filename, UnitCell *cell,
                         double min_res, double max_res,
                         int bij, const char *spg)
{
	RefList *reflist;
	char *sym_str;
	SymOpList *sym;
	int r;
	char *crystal_name;
	char *dir_name;
	char *dir_basename;

	reflist = read_reflections_2(result->hkl, &sym_str);
	if ( reflist == NULL ) return 1;
	if ( sym_str == NULL ) return 1;

	sym = get_pointgroup(sym_str);
	if ( sym == NULL ) return 1;

	dir_name = getcwd(NULL, 0);
	dir_basename = safe_basename(dir_name);
	crystal_name = result->indexing_result_name;
	if ( crystal_name == NULL ) {
		crystal_name = "unknown";
	}
	r = write_to_mtz(reflist, sym, cell, min_res, max_res, filename,
	                 result->name, crystal_name, dir_basename, bij,
			 spg);

	free(dir_name);
	free(dir_basename);
	free_symoplist(sym);
	free(sym_str);
	reflist_free(reflist);

	return r;
}


static void reminder(GtkWidget *window, const char *message)
{
	GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(window),
	                                           GTK_DIALOG_DESTROY_WITH_PARENT,
	                                           GTK_MESSAGE_INFO,
	                                           GTK_BUTTONS_OK,
	                                           NULL);
	gtk_message_dialog_set_markup(GTK_MESSAGE_DIALOG(dialog), message);
	gtk_window_set_title(GTK_WINDOW(dialog), "Reminder");
	gtk_dialog_run(GTK_DIALOG(dialog));
	gtk_widget_destroy(dialog);
}


static int export_data(struct export_window *win, char *filename)
{
	gchar *cell_filename;
	const char *dataset;
	const char *format;
	struct gui_merge_result *result;
	UnitCell *cell;
	int r = 0;
	double min_res = 0;
	double max_res = +INFINITY;

	dataset = gtk_combo_box_get_active_id(GTK_COMBO_BOX(win->dataset));
	if ( dataset == NULL ) {
		error_box(win->proj, "Please select the dataset to export.\n");
		return 1;
	}

	format = gtk_combo_box_get_active_id(GTK_COMBO_BOX(win->format));
	if ( format == NULL ) {
		error_box(win->proj, "Please select the data format to use.\n");
		return 1;
	}

	cell_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(win->cell_chooser));
	if ( cell_filename == NULL ) {
		error_box(win->proj, "Please choose the unit cell file.\n");
		return 1;
	}

	cell = load_cell_from_file(cell_filename);
	if ( cell == NULL ) {
		ERROR("Failed to load unit cell file %s\n", cell_filename);
		return 1;
	}

	if ( get_bool(win->limit_res) ) {
		min_res = 1e10/get_float(win->min_res);
		max_res = 1e10/get_float(win->max_res);
	}

	result = find_merge_result_by_name(win->proj, dataset);
	if ( result == NULL ) {
		ERROR("Couldn't find merged dataset '%s'\n", dataset);
		return 1;
	}

	STATUS("Exporting dataset %s to %s, in format %s, using unit cell %s,"
	       "%f to %f m^-1\n", dataset, filename, format, cell_filename,
	       min_res, max_res);

	if ( strcmp(format, "mtz") == 0 ) {
		r = export_to_mtz(result, filename, cell, min_res, max_res, 0, NULL);
		reminder(win->window,
		         "CrystFEL <b>does not know</b> the space "
		         "group of your structure.\n\n"
		         "The space group written into the MTZ header is just "
		         "a representation of the processing done within "
		         "CrystFEL.\n\n"
		         "You may need to change the space group in the MTZ "
		         "header, to allow structure solution programs to find "
		         "the correct space group.\n\n"
		         "<a href=\"https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/pointgroup.rst\">"
		         "Click here for some background information</a>");
	} else if ( strcmp(format, "mtz-bij") == 0 ) {
		r = export_to_mtz(result, filename, cell, min_res, max_res, 1, NULL);
		reminder(win->window,
		         "CrystFEL <b>does not know</b> the space "
		         "group of your structure.\n\n"
		         "The space group written into the MTZ header is just "
		         "a representation of the processing done within "
		         "CrystFEL.\n\n"
		         "You may need to change the space group in the MTZ "
		         "header, to allow structure solution programs to find "
		         "the correct space group.\n\n"
		         "<a href=\"https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/pointgroup.rst\">"
		         "Click here for some background information</a>");
	} else if ( strcmp(format, "xds") == 0 ) {
		r = export_to_xds(result, filename, cell, min_res, max_res);
		reminder(win->window,
		         "CrystFEL <b>does not know</b> the space "
		         "group of your structure.\n\n"
		         "The space group written into the XDS file "
		         "(header SPACE_GROUP_NUMBER=...) is just a "
		         "representation of the processing done within "
		         "CrystFEL.\n\n"
		         "You may need to change the space group in the XDS "
		         "header, to allow structure solution programs to find "
		         "the correct space group.\n\n"
		         "<a href=\"https://gitlab.desy.de/thomas.white/crystfel/-/blob/master/doc/articles/pointgroup.rst\">"
		         "Click here for some background information</a>");
	} else {
		ERROR("Unrecognised export format '%s'\n", format);
		return 1;
	}

	if ( r ) {
		ERROR("Export failed\n");
	}

	g_free(cell_filename);

	return 0;
}


static void export_response_sig(GtkWidget *dialog, gint resp,
                                struct export_window *win)
{
	int r = 0;

	if ( resp == GTK_RESPONSE_ACCEPT ) {
		gchar *filename;
		filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
		r = export_data(win, filename);
		g_free(filename);
	}

	if ( !r ) gtk_widget_destroy(dialog);
}


static void format_changed_sig(GtkWidget *format, struct export_window *win)
{
	const char *new_format;
	gchar *cur_name;
	new_format = gtk_combo_box_get_active_id(GTK_COMBO_BOX(format));
	cur_name = gtk_file_chooser_get_current_name(GTK_FILE_CHOOSER(win->window));
	strip_extension(cur_name);
	if ( cur_name[0] == '\0' ) {
		if ( strncmp(new_format, "mtz", 3) == 0 ) {
			gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(win->window),
			                                  "data_from_crystfel.mtz");
		} else if ( strcmp(new_format, "xds") == 0 ) {
			gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(win->window),
			                                  "XDS_ASCII.HKL");
		}
	} else {
		gchar *new_name = malloc(strlen(cur_name)+5);
		if ( new_name != NULL ) {
			strcpy(new_name, cur_name);
			if ( strncmp(new_format, "mtz", 3) == 0 ) {
				strcat(new_name, ".mtz");
			} else if ( strcmp(new_format, "xds") == 0 ) {
				strcat(new_name, ".HKL");
			}
			gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(win->window),
			                                  new_name);
			free(new_name);
		}
	}
	g_free(cur_name);
}


/* Dataset changed, update unit cell if appropriate */
static void ds_changed_sig(GtkWidget *widget, struct export_window *win)
{
	const char *dataset;
	struct gui_merge_result *result;
	char *fn;

	dataset = gtk_combo_box_get_active_id(GTK_COMBO_BOX(win->dataset));
	if ( dataset == NULL ) return;

	result = find_merge_result_by_name(win->proj, dataset);
	if ( result == NULL ) return;

	fn = cell_file_for_result(result);
	if ( fn != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(win->cell_chooser), fn);
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
	int i;

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
	win->window = dialog;

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
	for ( i=0; i<proj->n_merge_results; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->dataset),
		                          proj->merge_results[i].name,
		                          proj->merge_results[i].name);
	}
	g_signal_connect(G_OBJECT(win->dataset), "changed",
	                 G_CALLBACK(ds_changed_sig), win);

	label = gtk_label_new("Format");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->format = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(win->format),
	                   FALSE, FALSE, 4.0);
	if ( libcrystfel_can_write_mtz() ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "mtz",
		                          "MTZ");
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "mtz-bij",
		                          "MTZ, Bijvoet pairs together");
	}
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(win->format), "xds",
	                          "XDS ASCII");
	gtk_combo_box_set_active(GTK_COMBO_BOX(win->format), 0);
	g_signal_connect(G_OBJECT(win->format), "changed",
	                 G_CALLBACK(format_changed_sig), win);
	format_changed_sig(win->format, win);

	label = gtk_label_new("Unit cell file:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	win->cell_chooser = gtk_file_chooser_button_new("Unit cell file",
	                                              GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(win->cell_chooser),
	                                TRUE);
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
	g_signal_connect(G_OBJECT(win->limit_res), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->min_res);
	g_signal_connect(G_OBJECT(win->limit_res), "toggled",
	                 G_CALLBACK(i_maybe_disable), win->max_res);
	gtk_widget_set_sensitive(win->min_res, FALSE);
	gtk_widget_set_sensitive(win->max_res, FALSE);

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_CLOSE);
	gtk_widget_show_all(dialog);

	return FALSE;
}
