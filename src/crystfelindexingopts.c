/*
 * crystfelindexingopts.h
 *
 * A GTK widget to set indexing options
 *
 * Copyright © 2020-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020-2021 Thomas White <taw@physics.org>
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gtk/gtk.h>
#include <glib-object.h>
#include <errno.h>
#include <math.h>

#include <integration.h>
#include <index.h>

#include "crystfelindexingopts.h"
#include "gtk-util-routines.h"


G_DEFINE_TYPE(CrystFELIndexingOpts,
              crystfel_indexing_opts,
              GTK_TYPE_NOTEBOOK)


static void crystfel_indexing_opts_class_init(CrystFELIndexingOptsClass *klass)
{
}


static void crystfel_indexing_opts_init(CrystFELIndexingOpts *io)
{
}


static void i_disable_if_not(GtkWidget *toggle, GtkWidget *widget)
{
	i_maybe_disable(toggle, widget);
	g_signal_connect(G_OBJECT(toggle),
	                 "toggled",
	                 G_CALLBACK(i_maybe_disable),
	                 widget);
}


static GtkWidget *add_tol(GtkGrid *grid, const char *spec_t,
                          const char *unit_t, gint left, gint top)
{
	GtkWidget *spec;
	GtkWidget *entry;
	GtkWidget *unit;

	spec = gtk_label_new(spec_t);
	g_object_set(G_OBJECT(spec), "margin-left", 12, NULL);
	gtk_grid_attach(grid, spec, left, top, 1, 1);

	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 5);
	gtk_grid_attach(grid, entry, left+1, top, 1, 1);

	unit = gtk_label_new(unit_t);
	g_object_set(G_OBJECT(unit), "margin-right", 12, NULL);
	gtk_grid_attach(grid, unit, left+2, top, 1, 1);

	return entry;
}


static GtkWidget *make_tolerances(CrystFELIndexingOpts *io)
{
	GtkWidget *grid;

	grid = gtk_grid_new();
	gtk_grid_set_row_spacing(GTK_GRID(grid), 4);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 4);
	gtk_container_set_border_width(GTK_CONTAINER(grid), 6);

	io->tols[0] = add_tol(GTK_GRID(grid), "a", "%", 0, 0);
	io->tols[1] = add_tol(GTK_GRID(grid), "b", "%", 4, 0);
	io->tols[2] = add_tol(GTK_GRID(grid), "c", "%", 8, 0);
	io->tols[3] = add_tol(GTK_GRID(grid), "α", "°", 0, 1);
	io->tols[4] = add_tol(GTK_GRID(grid), "β", "°", 4, 1);
	io->tols[5] = add_tol(GTK_GRID(grid), "ɣ", "°", 8, 1);

	return grid;
}


static void add_method(GtkListStore *store, const char *name,
                       const char *friendly_name,
                       gboolean enabled, gboolean prior_cell,
                       gboolean prior_latt)
{
	GtkTreeIter iter;
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter,
	                   0, enabled,
	                   1, friendly_name,
	                   2, prior_cell,
	                   3, prior_latt,
	                   4, name,
	                   -1);
}


static void toggle_column(GtkListStore *store,
                          gchar *path_str, int column)
{
	GtkTreePath *path;
	GtkTreeIter iter;
	gboolean val;

	path = gtk_tree_path_new_from_string(path_str);
	gtk_tree_model_get_iter(GTK_TREE_MODEL(store),
	                        &iter, path);
	gtk_tree_path_free(path);

	gtk_tree_model_get(GTK_TREE_MODEL(store),
	                   &iter, column, &val, -1);

	gtk_list_store_set(store, &iter, column, !val, -1);
}


static void indm_toggled(GtkCellRendererToggle *cr,
                         gchar *path_str,
                         CrystFELIndexingOpts *io)
{
	toggle_column(io->indm_store, path_str, 0);
}


static void prior_cell_toggled(GtkCellRendererToggle *cr,
                               gchar *path_str,
                               CrystFELIndexingOpts *io)
{
	toggle_column(io->indm_store, path_str, 2);
}


static void prior_latt_toggled(GtkCellRendererToggle *cr,
                               gchar *path_str,
                               CrystFELIndexingOpts *io)
{
	toggle_column(io->indm_store, path_str, 3);
}


static GtkWidget *make_indexing_methods(CrystFELIndexingOpts *io)
{
	GtkWidget *treeview;
	GtkCellRenderer *renderer;
	GtkTreeViewColumn *column;

	io->indm_store = gtk_list_store_new(5,
	                           G_TYPE_BOOLEAN,  /* Algo on */
	                           G_TYPE_STRING,   /* Friendly name */
	                           G_TYPE_BOOLEAN,  /* Prior cell */
	                           G_TYPE_BOOLEAN,  /* Prior latt */
	                           G_TYPE_STRING);  /* Real name */

	add_method(io->indm_store, "dirax", "DirAx", TRUE, FALSE, FALSE);
	add_method(io->indm_store, "mosflm", "MOSFLM", TRUE, TRUE, TRUE);
	add_method(io->indm_store, "xds", "XDS", TRUE, TRUE, TRUE);
	add_method(io->indm_store, "xgandalf", "XGANDALF", TRUE, TRUE, FALSE);
	add_method(io->indm_store, "pinkIndexer", "PinkIndexer", FALSE, TRUE, FALSE);
	add_method(io->indm_store, "taketwo", "TakeTwo", TRUE, TRUE, FALSE);
	add_method(io->indm_store, "asdf", "ASDF", TRUE, TRUE, FALSE);
	add_method(io->indm_store, "felix", "Felix", FALSE, TRUE, FALSE);

	treeview = gtk_tree_view_new_with_model(GTK_TREE_MODEL(io->indm_store));

	renderer = gtk_cell_renderer_toggle_new();
	column = gtk_tree_view_column_new_with_attributes(NULL,
	                                                  renderer,
	                                                  "active", 0,
	                                                  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(treeview), column);
	g_signal_connect(G_OBJECT(renderer), "toggled",
	                 G_CALLBACK(indm_toggled), io);

	renderer = gtk_cell_renderer_text_new();
	column = gtk_tree_view_column_new_with_attributes("Method",
	                                                  renderer,
	                                                  "text", 1,
	                                                  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(treeview), column);

	renderer = gtk_cell_renderer_toggle_new();
	column = gtk_tree_view_column_new_with_attributes("Prior unit cell",
	                                                  renderer,
	                                                  "active", 2,
	                                                  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(treeview), column);
	g_signal_connect(G_OBJECT(renderer), "toggled",
	                 G_CALLBACK(prior_cell_toggled), io);

	renderer = gtk_cell_renderer_toggle_new();
	column = gtk_tree_view_column_new_with_attributes("Prior lattice type",
	                                                  renderer,
	                                                  "active", 3,
	                                                  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(treeview), column);
	g_signal_connect(G_OBJECT(renderer), "toggled",
	                 G_CALLBACK(prior_latt_toggled), io);

	return treeview;
}


static void cell_file_clear_sig(GtkButton *buton,
                                CrystFELIndexingOpts *io)
{
	io->cell_file = NULL;
	gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(io->cell_chooser),
	                              "(none)");
}


static void auto_indm_toggle_sig(GtkToggleButton *togglebutton,
                                 CrystFELIndexingOpts *io)
{
	gtk_widget_set_sensitive(GTK_WIDGET(io->indm_chooser),
	                         !gtk_toggle_button_get_active(togglebutton));
}


static void cell_file_set_sig(GtkFileChooserButton *widget,
                              CrystFELIndexingOpts *io)
{
	gchar *filename;
	filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(widget));
	g_free(io->cell_file);
	io->cell_file = filename;
}


static GtkWidget *indexing_parameters(CrystFELIndexingOpts *io)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *expander;
	GtkWidget *frame;
	GtkWidget *tolerances;
	GtkWidget *button;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	/* Cell file chooser */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Unit cell file:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	io->cell_chooser = gtk_file_chooser_button_new("Unit cell file",
	                                               GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_file_chooser_set_local_only(GTK_FILE_CHOOSER(io->cell_chooser),
	                                TRUE);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->cell_chooser),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(io->cell_chooser), "file-set",
	                 G_CALLBACK(cell_file_set_sig), io);
	button = gtk_button_new_from_icon_name("edit-clear",
	                                       GTK_ICON_SIZE_BUTTON);
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(cell_file_clear_sig), io);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(hbox, "-p / --pdb (or CrystFEL format)");

	/* Indexing method selector */
	io->auto_indm = gtk_check_button_new_with_label("Automatically choose the indexing methods");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->auto_indm),
	                   FALSE, FALSE, 4.0);
	g_signal_connect(G_OBJECT(io->auto_indm), "toggled",
	                 G_CALLBACK(auto_indm_toggle_sig), io);
	expander = gtk_expander_new("Select indexing methods and prior information");
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(expander), GTK_WIDGET(frame));
	io->indm_chooser = make_indexing_methods(io);
	gtk_container_add(GTK_CONTAINER(frame),
	                  GTK_WIDGET(io->indm_chooser));
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(expander),
	                   FALSE, FALSE, 4.0);
	gtk_container_set_border_width(GTK_CONTAINER(frame), 6);

	/* --multi */
	io->multi = gtk_check_button_new_with_label("Attempt to find multiple lattices per frame");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->multi),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->multi, "--multi");

	/* --no-refine (NB inverse) */
	io->refine = gtk_check_button_new_with_label("Refine the indexing solution");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->refine),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->refine, "--no-refine");

	/* --no-retry (NB inverse) */
	io->retry = gtk_check_button_new_with_label("Retry indexing if unsuccessful");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->retry),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->retry, "--no-retry");

	/* --no-check-peaks (NB inverse) */
	io->check_peaks = gtk_check_button_new_with_label("Check indexing solutions match peaks");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->check_peaks),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->check_peaks, "--no-check-peaks");

	/* --no-check-cell (NB inverse) and --tolerance */
	io->check_cell = gtk_check_button_new_with_label("Check indexing solutions against reference cell");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->check_cell),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->check_cell, "--no-check-cell");
	expander = gtk_expander_new("Unit cell tolerances");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(expander),
	                   FALSE, FALSE, 4.0);
	tolerances = make_tolerances(io);
	gtk_container_add(GTK_CONTAINER(expander), tolerances);
	i_disable_if_not(io->check_cell, tolerances);
	gtk_widget_set_tooltip_text(expander, "--tolerance");

	/* --min-peaks (NB add one) */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	io->enable_hitfind = gtk_check_button_new_with_label("Skip frames with fewer than");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->enable_hitfind),
	                   FALSE, FALSE, 0.0);
	io->ignore_fewer_peaks = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->ignore_fewer_peaks), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->ignore_fewer_peaks),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("peaks");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	i_disable_if_not(io->enable_hitfind, io->ignore_fewer_peaks);
	gtk_widget_set_tooltip_text(hbox, "--min-peaks");

	gtk_widget_show_all(box);

	return box;
}


static GtkWidget *integration_parameters(CrystFELIndexingOpts *io)
{
	GtkWidget *box;
	GtkWidget *label;
	GtkWidget *hbox;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	/* --integration=method */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Integration method:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	io->integration_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->integration_combo),
	                   FALSE, FALSE, 4.0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(io->integration_combo), "none",
	                "No integration (only spot prediction)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(io->integration_combo), "rings",
	                "Ring summation");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(io->integration_combo), "prof2d",
	                "Two dimensional profile fitting");
	gtk_widget_set_tooltip_text(hbox, "--integration={none,rings,prof2d}");

	/* -cen */
	io->centering = gtk_check_button_new_with_label("Center integration boxes on observed reflections");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->centering),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->centering, "--integration=xxx-cen");

	/* --overpredict */
	io->overpredict = gtk_check_button_new_with_label("Over-predict reflections (for post-refinement)");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->overpredict),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->overpredict, "--overpredict");

	/* --push-res */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	io->limit_res = gtk_check_button_new_with_label("Limit prediction to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->limit_res),
	                   FALSE, FALSE, 0.0);
	io->push_res = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->push_res), 6);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->push_res),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("nm-1 above apparent resolution limit");
	gtk_label_set_markup(GTK_LABEL(label),
	                     "nm<sup>-1</sup> above apparent resolution limit");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	i_disable_if_not(io->limit_res, io->push_res);
	gtk_widget_set_tooltip_text(io->limit_res, "--push-res");

	/* --int-radii */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);

	label = gtk_label_new("Integration radii - inner:");
	gtk_widget_set_tooltip_text(label, "--int-radius=inner,middle,outer");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	io->ir_inn = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->ir_inn), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->ir_inn),
	                   FALSE, FALSE, 4.0);

	label = gtk_label_new("middle:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	io->ir_mid = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->ir_mid), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->ir_mid),
	                   FALSE, FALSE, 4.0);

	label = gtk_label_new("outer:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	io->ir_out = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->ir_out), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->ir_out),
	                   FALSE, FALSE, 4.0);

	/* --fix-profile-radius */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	io->fix_profile_radius_p = gtk_check_button_new_with_label("Fix the reflection radius to");
	gtk_box_pack_start(GTK_BOX(hbox),
	                   GTK_WIDGET(io->fix_profile_radius_p),
	                   FALSE, FALSE, 0.0);
	io->fix_profile_radius = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->fix_profile_radius), 6);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->fix_profile_radius),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("nm-1");
	gtk_label_set_markup(GTK_LABEL(label), "nm<sup>-1</sup>");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	i_disable_if_not(io->fix_profile_radius_p, io->fix_profile_radius);
	gtk_widget_set_tooltip_text(io->fix_profile_radius,
	                            "--fix-profile-radius");

	/* --fix-divergence */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("Beam divergence (full angle) for spot prediction:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label), FALSE, FALSE, 0.0);
	io->fix_divergence = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(io->fix_divergence), 6);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(io->fix_divergence),
	                   FALSE, FALSE, 4.0);
	label = gtk_label_new("mrad");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(hbox, "--fix-divergence");
	gtk_widget_show_all(box);

	return box;
}


static GtkWidget *advanced_parameters(CrystFELIndexingOpts *io)
{
	GtkWidget *box;
	GtkWidget *expander;
	GtkWidget *label;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	expander = gtk_expander_new("PinkIndexer");
	gtk_box_pack_start(GTK_BOX(box), expander, FALSE, FALSE, 8);
	label = gtk_label_new("Advanced options for this indexing method are "
	                      "currently not available through the GUI.");
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_container_add(GTK_CONTAINER(expander), label);

	expander = gtk_expander_new("XGandalf");
	gtk_box_pack_start(GTK_BOX(box), expander, FALSE, FALSE, 8);
	label = gtk_label_new("Advanced options for this indexing method are "
	                      "currently not available through the GUI.");
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_container_add(GTK_CONTAINER(expander), label);

	expander = gtk_expander_new("Felix");
	gtk_box_pack_start(GTK_BOX(box), expander, FALSE, FALSE, 8);
	label = gtk_label_new("Advanced options for this indexing method are "
	                      "currently not available through the GUI.");
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_container_add(GTK_CONTAINER(expander), label);

	expander = gtk_expander_new("TakeTwo");
	gtk_box_pack_start(GTK_BOX(box), expander, FALSE, FALSE, 8);
	label = gtk_label_new("Advanced options for this indexing method are "
	                      "currently not available through the GUI.");
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_container_add(GTK_CONTAINER(expander), label);

	expander = gtk_expander_new("FromFile");
	gtk_box_pack_start(GTK_BOX(box), expander, FALSE, FALSE, 8);
	label = gtk_label_new("Advanced options for this indexing method are "
	                      "currently not available through the GUI.");
	gtk_label_set_line_wrap(GTK_LABEL(label), TRUE);
	gtk_container_add(GTK_CONTAINER(expander), label);

	return box;
}


static void add_metadata_item(GtkListStore *model, const char *string)
{
	GtkTreeIter iter;
	gtk_list_store_append(model, &iter);
	gtk_list_store_set(model, &iter, 0, strdup(string), -1);
}


static gboolean add_metadata_sig(GtkWidget *button, GtkListStore *model)
{
	add_metadata_item(model, "/instrument/something");
	return FALSE;
}


static gboolean edit_metadata_sig(GtkCellRendererText *cell,
                                  const gchar *path_str,
                                  const gchar *new_text,
                                  GtkListStore *model)
{
	GtkTreeIter iter;
	gchar *old_text;
	GtkTreePath *path = gtk_tree_path_new_from_string(path_str);
	gtk_tree_model_get_iter(GTK_TREE_MODEL(model), &iter, path);
	gtk_tree_path_free(path);
	gtk_tree_model_get(GTK_TREE_MODEL(model), &iter, 0, &old_text, -1);
        g_free(old_text);
	gtk_list_store_set(model, &iter,
	                   0, g_strdup(new_text),
	                   -1);
	return FALSE;
}


static gboolean remove_metadata_sig(GtkWidget *button,
                                    GtkTreeView *treeview)
{
	GtkTreeSelection *sel;
	GtkTreeIter iter;
	GtkTreeModel *model;

	model = gtk_tree_view_get_model(treeview);
	sel = gtk_tree_view_get_selection(treeview);
	if ( gtk_tree_selection_get_selected(sel, NULL, &iter) != 0 ) {
		gchar *old_text;
		gtk_tree_model_get(GTK_TREE_MODEL(model), &iter, 0, &old_text, -1);
		g_free(old_text);
		gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
	}

	return FALSE;
}


static GtkWidget *stream_parameters(CrystFELIndexingOpts *io)
{
	GtkWidget *box;
	GtkWidget *treeview;
	GtkCellRenderer *renderer;
	GtkWidget *button;
	GtkWidget *hbox;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	/* --no-non-hits-in-stream */
	io->exclude_nonhits = gtk_check_button_new_with_label("Exclude skipped frames from stream");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->exclude_nonhits),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->exclude_nonhits,
	                            "--no-non-hits-in-stream");

	/* --no-peaks-in-stream */
	io->no_peaks_in_stream = gtk_check_button_new_with_label("Exclude peak search results from stream");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->no_peaks_in_stream),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->no_peaks_in_stream,
	                            "--no-peaks-in-stream");

	/* --no-refls-in-stream */
	io->no_refls_in_stream = gtk_check_button_new_with_label("Exclude integrated intensities from stream");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(io->no_refls_in_stream),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(io->no_refls_in_stream,
	                            "--no-refls-in-stream");

	io->copy_metadata_store = gtk_list_store_new(1, G_TYPE_STRING);

	treeview = gtk_tree_view_new_with_model(GTK_TREE_MODEL(io->copy_metadata_store));

	renderer = gtk_cell_renderer_text_new();
	g_object_set(renderer, "editable", TRUE, NULL);
	g_signal_connect(G_OBJECT(renderer), "edited",
	                 G_CALLBACK(edit_metadata_sig),
	                 io->copy_metadata_store);
	gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW(treeview),
	                                            -1,
	                                            "Metadata to copy to stream",
	                                            renderer,
	                                            "text", 0,
	                                            NULL);

	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(treeview),
	                   FALSE, FALSE, 4.0);
	gtk_widget_set_tooltip_text(box, "--copy-header / --copy-hdf5-field");

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
	gtk_container_set_border_width(GTK_CONTAINER(hbox), 8);
	button = gtk_button_new_from_icon_name("list-add", GTK_ICON_SIZE_BUTTON);
	gtk_button_set_label(GTK_BUTTON(button), "Add item");
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(add_metadata_sig), io->copy_metadata_store);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);
	button = gtk_button_new_from_icon_name("list-remove", GTK_ICON_SIZE_BUTTON);
	gtk_button_set_label(GTK_BUTTON(button), "Remove item");
	g_signal_connect(G_OBJECT(button), "clicked",
	                 G_CALLBACK(remove_metadata_sig), treeview);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(button),
	                   FALSE, FALSE, 4.0);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 4.0);

	gtk_widget_show_all(box);

	return box;
}


GtkWidget *crystfel_indexing_opts_new()
{
	CrystFELIndexingOpts *io;

	io = g_object_new(CRYSTFEL_TYPE_INDEXING_OPTS, NULL);

	io->cell_file = NULL;
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(io), GTK_POS_LEFT);

	gtk_notebook_append_page(GTK_NOTEBOOK(io),
	                         indexing_parameters(io),
	                         gtk_label_new("Indexing"));

	gtk_notebook_append_page(GTK_NOTEBOOK(io),
	                         integration_parameters(io),
	                         gtk_label_new("Integration"));

	gtk_notebook_append_page(GTK_NOTEBOOK(io),
	                         advanced_parameters(io),
	                         gtk_label_new("Advanced indexing"));

	io->stream_params = stream_parameters(io);
	gtk_notebook_append_page(GTK_NOTEBOOK(io), io->stream_params,
	                         gtk_label_new("Stream contents"));

	return GTK_WIDGET(io);
}


char *crystfel_indexing_opts_get_cell_file(CrystFELIndexingOpts *opts)
{
	return safe_strdup(opts->cell_file);
}


/* NULL means "automatic".
 * "none" means "no indexing" */
char *crystfel_indexing_opts_get_indexing_method_string(CrystFELIndexingOpts *opts)
{
	GtkTreePath *path;
	GtkTreeIter iter;
	char indm_str[1024];
	int first = 1;

	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->auto_indm)) ) {
		return NULL;
	}

	path = gtk_tree_path_new_from_string("0");
	gtk_tree_model_get_iter(GTK_TREE_MODEL(opts->indm_store),
	                        &iter, path);
	gtk_tree_path_free(path);

	indm_str[0] = '\0';
	do {
		gboolean enabled, prior_cell, prior_latt;
		gchar *name;

		gtk_tree_model_get(GTK_TREE_MODEL(opts->indm_store),
		                   &iter,
		                   0, &enabled,
		                   2, &prior_cell,
		                   3, &prior_latt,
		                   4, &name,
		                   -1);

		if ( enabled ) {
			if ( !first ) {
				strcat(indm_str, ",");
			}
			first = 0;
			strcat(indm_str, name);
			if ( prior_cell ) {
				strcat(indm_str, "-cell");
			}
			if ( prior_latt ) {
				strcat(indm_str, "-latt");
			}
		}

	} while ( gtk_tree_model_iter_next(GTK_TREE_MODEL(opts->indm_store),
	                                   &iter) );

	if ( indm_str[0] == '\0' ) {
		strcpy(indm_str, "none");
	}
	return strdup(indm_str);
}


int crystfel_indexing_opts_get_multi_lattice(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->multi));
}


int crystfel_indexing_opts_get_refine(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->refine));
}


int crystfel_indexing_opts_get_retry(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->retry));
}


int crystfel_indexing_opts_get_peak_check(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->check_peaks));
}


int crystfel_indexing_opts_get_cell_check(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->check_cell));
}



/* Values in 'tols' are in frac (not %) and rad */
void crystfel_indexing_opts_get_tolerances(CrystFELIndexingOpts *opts,
                                           float *tols)
{
	int i;
	for ( i=0; i<3; i++ ) {
		float tol;
		char *rval;
		const gchar *text = gtk_entry_get_text(GTK_ENTRY(opts->tols[i]));
		errno = 0;
		tol = strtod(text, &rval);
		if ( *rval != '\0' ) {
			printf("Invalid tolerance '%s'\n", text);
		} else {
			tols[i] = tol / 100.0;
		}
	}
	for ( i=3; i<6; i++ ) {
		float tol;
		char *rval;
		const gchar *text = gtk_entry_get_text(GTK_ENTRY(opts->tols[i]));
		errno = 0;
		tol = strtod(text, &rval);
		if ( *rval != '\0' ) {
			printf("Invalid tolerance '%s'\n", text);
		} else {
			tols[i] = deg2rad(tol);
		}
	}
}


int crystfel_indexing_opts_get_min_peaks(CrystFELIndexingOpts *opts)
{
	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->enable_hitfind)) ) {
		const gchar *text;
		int fewer_peaks;
		char *rval;
		text = gtk_entry_get_text(GTK_ENTRY(opts->ignore_fewer_peaks));
		errno = 0;
		fewer_peaks = strtod(text, &rval);
		if ( *rval != '\0' ) {
			printf("Invalid value for minimum number of peaks (%s)\n",
			      rval);
			return 0;
		}
		/* Subtract one because the dialog box says to skip
		 * frames with "FEWER THAN" this number */
		return fewer_peaks - 1;
	} else {
		return 0;
	}
}


char *crystfel_indexing_opts_get_integration_method_string(CrystFELIndexingOpts *opts)
{
	const gchar *id;
	char method[64];

	id = gtk_combo_box_get_active_id(GTK_COMBO_BOX(opts->integration_combo));
	if ( id == NULL ) return strdup("none");

	strcpy(method, id);
	if ( gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->centering)) ) {
		strcat(method, "-cen");
	}

	return strdup(method);
}


int crystfel_indexing_opts_get_overpredict(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->overpredict));
}


float crystfel_indexing_opts_get_push_res(CrystFELIndexingOpts *opts)
{
	if ( !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->limit_res)) ) {
		return INFINITY;
	} else {
		const gchar *text;
		float push_res;
		char *rval;
		text = gtk_entry_get_text(GTK_ENTRY(opts->push_res));
		errno = 0;
		push_res = strtof(text, &rval);
		if ( *rval != '\0' ) {
			printf("Invalid value for push-res (%s)\n",
			      rval);
			return INFINITY;
		}
		return push_res;
	}
}


void crystfel_indexing_opts_get_integration_radii(CrystFELIndexingOpts *opts,
                                                  float *ir_inn,
                                                  float *ir_mid,
                                                  float *ir_out)
{
	*ir_inn = get_float(opts->ir_inn);
	*ir_mid = get_float(opts->ir_mid);
	*ir_out = get_float(opts->ir_out);
}


int crystfel_indexing_opts_get_exclude_blanks(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->exclude_nonhits));
}


int crystfel_indexing_opts_get_exclude_peaks(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->no_peaks_in_stream));
}


int crystfel_indexing_opts_get_exclude_reflections(CrystFELIndexingOpts *opts)
{
	return gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(opts->no_refls_in_stream));
}


char **crystfel_indexing_opts_get_metadata_to_copy(CrystFELIndexingOpts *opts,
                                                   int *pn)
{
	GtkTreeIter iter;
	gboolean r;
	int n, i;
	char **arr;

	n = gtk_tree_model_iter_n_children(GTK_TREE_MODEL(opts->copy_metadata_store),
	                                   NULL);

	arr = malloc(n*sizeof(char *));
	if ( arr == NULL ) return NULL;

	r = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(opts->copy_metadata_store),
	                                  &iter);
	if ( r == FALSE ) return NULL;

	i = 0;
	do {
		gchar *header;
		gtk_tree_model_get(GTK_TREE_MODEL(opts->copy_metadata_store),
		                   &iter, 0, &header, -1);
		if ( i == n ) return NULL;
		arr[i++] = strdup(header);
	} while ( gtk_tree_model_iter_next(GTK_TREE_MODEL(opts->copy_metadata_store),
	                                   &iter) != FALSE );

	*pn = n;
	return arr;
}


double crystfel_indexing_opts_get_fixed_profile_radius(CrystFELIndexingOpts *opts,
                                                       int *active)
{
	*active = get_bool(opts->fix_profile_radius_p);
	return get_float(opts->fix_profile_radius)*1e9;
}


double crystfel_indexing_opts_get_fixed_divergence(CrystFELIndexingOpts *opts)
{
	return get_float(opts->fix_divergence)/1e3;
}


/********************** Setters *************************/


void crystfel_indexing_opts_set_show_stream_opts(CrystFELIndexingOpts *opts,
                                                 int val)
{
	opts->show_stream_opts = val;
	if ( val ) {
		gtk_widget_show_all(opts->stream_params);
	} else {
		gtk_widget_hide(opts->stream_params);
	}
}


void crystfel_indexing_opts_set_cell_file(CrystFELIndexingOpts *opts,
                                          const char *cell_file)
{
	if ( cell_file != NULL ) {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(opts->cell_chooser),
		                              cell_file);
		opts->cell_file = strdup(cell_file);
	} else {
		gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(opts->cell_chooser),
		                              "(null)");
		opts->cell_file = NULL;
	}
}


static const char *integration_method_id(IntegrationMethod meth)
{
	switch ( meth ) {
		case INTEGRATION_NONE : return "none";
		case INTEGRATION_RINGS : return "rings";
		case INTEGRATION_PROF2D : return "prof2d";
		default : return "none";
	}
}


void crystfel_indexing_opts_set_indexing_method_string(CrystFELIndexingOpts *opts,
                                                       const char *indm_str)
{
	GtkTreePath *path;
	GtkTreeIter iter;
	IndexingMethod *methods;
	int i, n;

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->auto_indm),
	                             (indm_str == NULL));
	if ( indm_str == NULL ) return;

	methods = parse_indexing_methods(indm_str, &n);
	if ( methods == NULL ) {
		ERROR("Failed to parse '%s'\n", indm_str);
		return;
	}

	path = gtk_tree_path_new_from_string("0");
	gtk_tree_model_get_iter(GTK_TREE_MODEL(opts->indm_store),
	                        &iter, path);
	gtk_tree_path_free(path);

	do {
		gchar *name;
		IndexingMethod this_method = 0;

		gtk_tree_model_get(GTK_TREE_MODEL(opts->indm_store),
		                   &iter,
		                   4, &name,
		                   -1);

		for ( i=0; i<n; i++ ) {
			char *str = base_indexer_str(methods[i]);
			if ( strcmp(str, name) == 0 ) {
				this_method = methods[i];
				break;
			}
		}

		gtk_list_store_set(opts->indm_store, &iter,
		                   0, (this_method != 0),
		                   2, (this_method & INDEXING_USE_CELL_PARAMETERS),
		                   3, (this_method & INDEXING_USE_LATTICE_TYPE),
		                   -1);

	} while ( gtk_tree_model_iter_next(GTK_TREE_MODEL(opts->indm_store),
	                                   &iter) );

	free(methods);
}


void crystfel_indexing_opts_set_multi_lattice(CrystFELIndexingOpts *opts,
                                              int multi)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->multi),
	                             multi);
}


void crystfel_indexing_opts_set_refine(CrystFELIndexingOpts *opts,
                                       int refine)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->refine),
	                             refine);
}


void crystfel_indexing_opts_set_retry(CrystFELIndexingOpts *opts,
                                      int retry)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->retry),
	                             retry);
}


void crystfel_indexing_opts_set_peak_check(CrystFELIndexingOpts *opts,
                                           int peak_check)

{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->check_peaks),
	                             peak_check);
}


void crystfel_indexing_opts_set_cell_check(CrystFELIndexingOpts *opts,
                                           int cell_check)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->check_cell),
	                             cell_check);
}


/* Values in 'tols' are in frac (not %) and rad */
void crystfel_indexing_opts_set_tolerances(CrystFELIndexingOpts *opts,
                                           float *tols)
{
	int i;
	for ( i=0; i<3; i++ ) {
		char tmp[64];
		snprintf(tmp, 63, "%f", tols[i]*100.0);
		gtk_entry_set_text(GTK_ENTRY(opts->tols[i]), tmp);
	}
	for ( i=3; i<6; i++ ) {
		char tmp[64];
		snprintf(tmp, 63, "%f", rad2deg(tols[i]));
		gtk_entry_set_text(GTK_ENTRY(opts->tols[i]), tmp);
	}
}


void crystfel_indexing_opts_set_min_peaks(CrystFELIndexingOpts *opts,
                                          int min_peaks)
{
	char tmp[64];

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->enable_hitfind),
	                             (min_peaks > 0));

	/* Plus one because dialog says skip when "fewer than" X peaks */
	snprintf(tmp, 63, "%i", min_peaks+1);
	gtk_entry_set_text(GTK_ENTRY(opts->ignore_fewer_peaks), tmp);
}


void crystfel_indexing_opts_set_integration_method_string(CrystFELIndexingOpts *opts,
                                                          const char *integr_str)
{
	IntegrationMethod meth;
	int err;

	meth = integration_method(integr_str, &err);
	if ( !err ) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->centering),
		                             meth & INTEGRATION_CENTER);
		gtk_combo_box_set_active_id(GTK_COMBO_BOX(opts->integration_combo),
		                            integration_method_id(meth & INTEGRATION_METHOD_MASK));
	}
}


void crystfel_indexing_opts_set_overpredict(CrystFELIndexingOpts *opts,
                                            int overpredict)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->overpredict),
	                             overpredict);
}


void crystfel_indexing_opts_set_push_res(CrystFELIndexingOpts *opts,
                                         float push_res)
{
	char tmp[64];

	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->limit_res),
	                             !isinf(push_res));

	snprintf(tmp, 63, "%f", push_res);
	gtk_entry_set_text(GTK_ENTRY(opts->push_res), tmp);
}


void crystfel_indexing_opts_set_integration_radii(CrystFELIndexingOpts *opts,
                                                  float ir_inn,
                                                  float ir_mid,
                                                  float ir_out)
{
	char tmp[64];

	snprintf(tmp, 63, "%.1f", ir_inn);
	gtk_entry_set_text(GTK_ENTRY(opts->ir_inn), tmp);

	snprintf(tmp, 63, "%.1f", ir_mid);
	gtk_entry_set_text(GTK_ENTRY(opts->ir_mid), tmp);

	snprintf(tmp, 63, "%.1f", ir_out);
	gtk_entry_set_text(GTK_ENTRY(opts->ir_out), tmp);
}


void crystfel_indexing_opts_set_metadata_to_copy(CrystFELIndexingOpts *opts,
                                                 char *const *headers,
                                                 int n)
{
	int i;
	gtk_list_store_clear(opts->copy_metadata_store);
	if ( headers == NULL ) return;
	for ( i=0; i<n; i++ ) {
		add_metadata_item(opts->copy_metadata_store,
		                  headers[i]);
	}
}


void crystfel_indexing_opts_set_exclude_blanks(CrystFELIndexingOpts *opts,
                                               int flag)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->exclude_nonhits),
	                             flag);
}


void crystfel_indexing_opts_set_exclude_peaks(CrystFELIndexingOpts *opts,
                                              int flag)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->no_peaks_in_stream),
	                             flag);
}


void crystfel_indexing_opts_set_exclude_reflections(CrystFELIndexingOpts *opts,
                                                    int flag)
{
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(opts->no_refls_in_stream),
	                             flag);
}


void crystfel_indexing_opts_set_fixed_profile_radius(CrystFELIndexingOpts *opts,
                                                     int active,
                                                     double val)
{
	char tmp[64];
	set_active(opts->fix_profile_radius_p, active);
	snprintf(tmp, 63, "%.3f", val/1e9);
	gtk_entry_set_text(GTK_ENTRY(opts->fix_profile_radius), tmp);
}


void crystfel_indexing_opts_set_fixed_divergence(CrystFELIndexingOpts *opts,
                                                 double val)
{
	char tmp[64];
	snprintf(tmp, 63, "%.3f", val*1e3);
	gtk_entry_set_text(GTK_ENTRY(opts->fix_divergence), tmp);
}
