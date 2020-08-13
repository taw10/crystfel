/*
 * crystfelindexingopts.h
 *
 * A GTK widget to set indexing options
 *
 * Copyright © 2020 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *  2020 Thomas White <taw@physics.org>
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

#include "crystfelindexingopts.h"


G_DEFINE_TYPE(CrystFELIndexingOpts,
              crystfel_indexing_opts,
              GTK_TYPE_NOTEBOOK)


static void crystfel_indexing_opts_class_init(CrystFELIndexingOptsClass *klass)
{
}


static void crystfel_indexing_opts_init(CrystFELIndexingOpts *io)
{
}


static void add_method(GtkListStore *store, const char *name)
{
	GtkTreeIter iter;
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter,
	                   0, FALSE,
	                   1, name,
	                   2, FALSE,
	                   3, FALSE,
	                   -1);
}


static void add_tol(GtkGrid *grid, const char *spec_t,
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
}


static GtkWidget *make_tolerances(CrystFELIndexingOpts *io)
{
	GtkWidget *grid;

	grid = gtk_grid_new();
	gtk_grid_set_row_spacing(GTK_GRID(grid), 4);
	gtk_grid_set_column_spacing(GTK_GRID(grid), 4);
	gtk_container_set_border_width(GTK_CONTAINER(grid), 6);

	add_tol(GTK_GRID(grid), "a", "%", 0, 0);
	add_tol(GTK_GRID(grid), "b", "%", 4, 0);
	add_tol(GTK_GRID(grid), "c", "%", 8, 0);
	add_tol(GTK_GRID(grid), "α", "°", 0, 1);
	add_tol(GTK_GRID(grid), "β", "°", 4, 1);
	add_tol(GTK_GRID(grid), "ɣ", "°", 8, 1);

	return grid;
}


static GtkWidget *make_indexing_methods(CrystFELIndexingOpts *io)
{
	GtkWidget *treeview;
	GtkListStore *store;
	GtkCellRenderer *renderer;
	GtkTreeViewColumn *column;

	store = gtk_list_store_new(4,
	                           G_TYPE_BOOLEAN,  /* Algo on */
	                           G_TYPE_STRING,   /* Algo name */
	                           G_TYPE_BOOLEAN,  /* Prior cell */
	                           G_TYPE_BOOLEAN); /* Prior latt */

	add_method(store, "DirAx");
	add_method(store, "MOSFLM");
	add_method(store, "XDS");
	add_method(store, "XGANDALF");
	add_method(store, "PinkIndexer");
	add_method(store, "TakeTwo");
	add_method(store, "ASDF");
	add_method(store, "Felix");

	treeview = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));

	renderer = gtk_cell_renderer_toggle_new();
	column = gtk_tree_view_column_new_with_attributes(NULL,
	                                                  renderer,
	                                                  "active", 0,
	                                                  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(treeview), column);

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

	renderer = gtk_cell_renderer_toggle_new();
	column = gtk_tree_view_column_new_with_attributes("Prior lattice type",
	                                                  renderer,
	                                                  "active", 3,
	                                                  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(treeview), column);

	return treeview;
}


static GtkWidget *indexing_parameters(CrystFELIndexingOpts *io)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *entry;
	GtkWidget *check;
	GtkWidget *filechooser;
	GtkWidget *expander;
	GtkWidget *frame;
	GtkWidget *indexing_methods;
	GtkWidget *tolerances;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 4);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	/* Use unit cell / Cell file chooser */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	check = gtk_check_button_new_with_label("Use unit cell");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);
	filechooser = gtk_file_chooser_button_new("Unit cell file",
	                                          GTK_FILE_CHOOSER_ACTION_OPEN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(filechooser),
	                   FALSE, FALSE, 0);

	/* Indexing method selector */
	check = gtk_check_button_new_with_label("Automatically choose the indexing methods");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);
	expander = gtk_expander_new("Select indexing methods and prior information");
	frame = gtk_frame_new(NULL);
	gtk_frame_set_shadow_type(GTK_FRAME(frame), GTK_SHADOW_IN);
	gtk_container_add(GTK_CONTAINER(expander),
	                  GTK_WIDGET(frame));
	indexing_methods = make_indexing_methods(io);
	gtk_container_add(GTK_CONTAINER(frame),
	                  GTK_WIDGET(indexing_methods));
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(expander),
	                   FALSE, FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(frame), 6);

	/* --multi */
	check = gtk_check_button_new_with_label("Attempt to find multiple lattices per frame");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);

	/* --no-refine (NB inverse) */
	check = gtk_check_button_new_with_label("Refine the indexing solution");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);

	/* --no-retry (NB inverse) */
	check = gtk_check_button_new_with_label("Retry indexing if unsuccessful");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);

	/* --no-check-peaks (NB inverse) */
	check = gtk_check_button_new_with_label("Check indexing solutions match peaks");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);

	/* --no-check-cell (NB inverse) and --tolerance */
	check = gtk_check_button_new_with_label("Check indexing solutions against reference cell");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);
	expander = gtk_expander_new("Unit cell tolerances");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(expander),
	                   FALSE, FALSE, 0);
	tolerances = make_tolerances(io);
	gtk_container_add(GTK_CONTAINER(expander), tolerances);

	/* --min-peaks (NB add one) */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	check = gtk_check_button_new_with_label("Skip frames with fewer than");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("peaks");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	return box;
}


static GtkWidget *integration_parameters(CrystFELIndexingOpts *io)
{
	GtkWidget *box;
	GtkWidget *combo;
	GtkWidget *check;
	GtkWidget *label;
	GtkWidget *entry;
	GtkWidget *hbox;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	/* --integration=method */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Integration method:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(combo),
	                   FALSE, FALSE, 0);
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "none",
	                "No integration (only spot prediction)");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "rings",
	                "Ring summation");
	gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(combo), "prof2d",
	                "Two dimensional profile fitting");

	/* -cen */
	check = gtk_check_button_new_with_label("Center integration boxes on observed reflections");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);

	/* --overpredict */
	check = gtk_check_button_new_with_label("Over-predict reflections (for post-refinement)");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);

	/* --push-res */
	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	check = gtk_check_button_new_with_label("Limit prediction to");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(check),
	                   FALSE, FALSE, 0);
	entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 4);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(entry),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("nm-1 above apparent resolution limit");
	gtk_label_set_markup(GTK_LABEL(label),
	                     "nm<sup>-1</sup> above apparent resolution limit");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	/* FIXME: fix-bandwidth, divergence, profile-radius */
	/* FIXME: --int-radius */

	return box;
}


GtkWidget *crystfel_indexing_opts_new()
{
	CrystFELIndexingOpts *io;

	io = g_object_new(CRYSTFEL_TYPE_INDEXING_OPTS, NULL);

	gtk_notebook_append_page(GTK_NOTEBOOK(io),
	                         indexing_parameters(io),
	                         gtk_label_new("Indexing"));

	gtk_notebook_append_page(GTK_NOTEBOOK(io),
	                         integration_parameters(io),
	                         gtk_label_new("Integration"));

	gtk_widget_show_all(GTK_WIDGET(io));
	return GTK_WIDGET(io);
}
