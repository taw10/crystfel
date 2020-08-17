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
#include <cell.h>
#include <cell-utils.h>
#include <integration.h>

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
	proj->indexing_opts = NULL;
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

	if ( proj->indexing_opts != NULL ) return FALSE;

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


static void set_indexing_opts(struct crystfelproject *proj,
                              CrystFELIndexingOpts *opts)
{
	/* Indexing */
	crystfel_indexing_opts_set_cell_file(opts, proj->indexing_params.cell_file);
	crystfel_indexing_opts_set_indexing_method_string(opts, proj->indexing_params.indexing_methods);
	crystfel_indexing_opts_set_multi_lattice(opts, proj->indexing_params.multi);
	crystfel_indexing_opts_set_refine(opts, !proj->indexing_params.no_refine);
	crystfel_indexing_opts_set_retry(opts, !proj->indexing_params.no_retry);
	crystfel_indexing_opts_set_peak_check(opts, !proj->indexing_params.no_peak_check);
	crystfel_indexing_opts_set_cell_check(opts, !proj->indexing_params.no_cell_check);
	crystfel_indexing_opts_set_tolerances(opts, proj->indexing_params.tols);
	crystfel_indexing_opts_set_min_peaks(opts, proj->indexing_params.min_peaks);

	/* Integration */
	crystfel_indexing_opts_set_integration_method_string(opts, proj->indexing_params.integration_method);
	crystfel_indexing_opts_set_overpredict(opts, proj->indexing_params.overpredict);
	crystfel_indexing_opts_set_push_res(opts, proj->indexing_params.push_res);
}


static void get_indexing_opts(struct crystfelproject *proj,
                              CrystFELIndexingOpts *opts)
{
	/* Indexing */
	proj->indexing_params.cell_file = crystfel_indexing_opts_get_cell_file(opts);
	proj->indexing_params.indexing_methods = crystfel_indexing_opts_get_indexing_method_string(opts);
	proj->indexing_params.multi = crystfel_indexing_opts_get_multi_lattice(opts);
	proj->indexing_params.no_refine = !crystfel_indexing_opts_get_refine(opts);
	proj->indexing_params.no_retry = !crystfel_indexing_opts_get_retry(opts);
	proj->indexing_params.no_peak_check = !crystfel_indexing_opts_get_peak_check(opts);
	proj->indexing_params.no_cell_check = !crystfel_indexing_opts_get_cell_check(opts);
	proj->indexing_params.min_peaks = crystfel_indexing_opts_get_min_peaks(opts);

	/* Integration */
	proj->indexing_params.integration_method = crystfel_indexing_opts_get_integration_method_string(opts);
	proj->indexing_params.overpredict = crystfel_indexing_opts_get_overpredict(opts);
	proj->indexing_params.push_res = crystfel_indexing_opts_get_push_res(opts);
}


static IndexingFlags indexing_flags(struct index_params *params)
{
	IndexingFlags fl = 0;

	if ( !params->no_retry ) fl |= INDEXING_RETRY;
	if ( params->multi ) fl |= INDEXING_MULTI;
	if ( !params->no_refine ) fl |= INDEXING_REFINE;
	if ( !params->no_peak_check ) fl |= INDEXING_CHECK_PEAKS;
	if ( !params->no_cell_check ) fl |= INDEXING_CHECK_CELL;

	return fl;
}


static void run_indexing_once(struct crystfelproject *proj)
{
	IndexingPrivate *ipriv;
	UnitCell *cell;
	IntegrationMethod int_method;
	int i;
	int err;

	if ( proj->indexing_params.cell_file != NULL ) {
		cell = load_cell_from_file(proj->indexing_params.cell_file);
	} else {
		cell = NULL;
	}

	ipriv = setup_indexing(proj->indexing_params.indexing_methods,
	                       cell,
	                       proj->dtempl,
	                       proj->indexing_params.tols,
	                       indexing_flags(&proj->indexing_params),
	                       NULL, NULL, NULL, NULL);

	index_pattern(proj->cur_image, ipriv);

	for ( i=0; i<proj->cur_image->n_crystals; i++ ) {
		crystal_set_profile_radius(proj->cur_image->crystals[i], 0.02e9);
		crystal_set_mosaicity(proj->cur_image->crystals[i], 0.0);
	}

	err = 0;
	int_method = integration_method(proj->indexing_params.integration_method,
	                                &err);

	integrate_all_5(proj->cur_image, int_method, PMODEL_XSPHERE,
	                proj->indexing_params.push_res,
	                3, 4, 5,  /* FIXME */
	                INTDIAG_NONE, 0, 0, 0, NULL,
	                proj->indexing_params.overpredict);

	cleanup_indexing(ipriv);
}


static void index_one_response_sig(GtkWidget *dialog, gint resp,
                                   struct crystfelproject *proj)
{
	if ( resp == GTK_RESPONSE_OK ) {
		get_indexing_opts(proj,
		                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));
		run_indexing_once(proj);
	}

	gtk_widget_destroy(dialog);
	proj->indexing_opts = NULL;
}


gint index_one_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;

	if ( proj->indexing_opts != NULL ) return FALSE;

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

	proj->indexing_opts = crystfel_indexing_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox),
	                   GTK_WIDGET(proj->indexing_opts),
	                   FALSE, FALSE, 8.0);
	set_indexing_opts(proj,
	                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
}
