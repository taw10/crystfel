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
#include <predict-refine.h>

#include "crystfelimageview.h"
#include "crystfelindexingopts.h"
#include "gui_project.h"
#include "crystfel_gui.h"
#include "gui_peaksearch.h"

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


static int run_indexing_all(struct crystfelproject *proj,
                            int backend_idx, const char *job_title,
                            const char *job_notes)
{
	struct crystfel_backend *be;
	void *job_priv;

	be = &proj->backends[backend_idx];
	job_priv = be->run_indexing(job_title, job_notes, proj,
	                            be->indexing_opts_priv);

	if ( job_priv != NULL ) {
		add_running_task(proj, "Indexing all frames",
		                 be, job_priv);
		return 0;
	} else {
		return 1;
	}
}


struct new_job_params {
	struct crystfelproject *proj;
	GtkWidget *indexing_backend_combo;
	GtkWidget *indexing_backend_opts_widget;
	GtkWidget *indexing_backend_opts_box;
	GtkWidget *job_title_entry;
	GtkWidget *job_notes_text;
};


static char *get_all_text(GtkTextView *view)
{
	GtkTextBuffer *buf;
	GtkTextIter start, end;

	buf = gtk_text_view_get_buffer(view);

	gtk_text_buffer_get_start_iter(buf, &start);
	gtk_text_buffer_get_end_iter(buf, &end);

	return gtk_text_buffer_get_text(buf, &start, &end, FALSE);
}


static void index_all_response_sig(GtkWidget *dialog, gint resp,
                                   struct new_job_params *njp)
{
	if ( resp == GTK_RESPONSE_OK ) {

		int backend_idx;
		const char *job_title;
		char *job_notes;

		get_indexing_opts(njp->proj,
		                  CRYSTFEL_INDEXING_OPTS(njp->proj->indexing_opts));

		backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(njp->indexing_backend_combo));
		if ( backend_idx < 0 ) return;

		job_title = gtk_entry_get_text(GTK_ENTRY(njp->job_title_entry));
		job_notes = get_all_text(GTK_TEXT_VIEW(njp->job_notes_text));

		if ( job_title[0] == '\0' ) {
			ERROR("You must provide a job name.\n");
			return;
		}

		free(njp->proj->indexing_new_job_title);
		njp->proj->indexing_new_job_title = strdup(job_title);

		if ( run_indexing_all(njp->proj, backend_idx,
		                      job_title, job_notes) == 0 )
		{
			gtk_widget_destroy(dialog);
			njp->proj->indexing_opts = NULL;
		}

		free(job_notes);

	} else {
		gtk_widget_destroy(dialog);
		njp->proj->indexing_opts = NULL;
	}
}


static void indexing_backend_changed_sig(GtkWidget *combo,
                                         struct new_job_params *njp)
{
	int backend_idx;
	struct crystfel_backend *be;

	backend_idx = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
	if ( backend_idx < 0 ) return;
	njp->proj->indexing_backend_selected = backend_idx;

	be = &njp->proj->backends[backend_idx];

	if ( njp->indexing_backend_opts_widget != NULL ) {
		gtk_widget_destroy(njp->indexing_backend_opts_widget);
	}

	njp->indexing_backend_opts_widget = be->make_indexing_parameters_widget(be->indexing_opts_priv);

	gtk_box_pack_start(GTK_BOX(njp->indexing_backend_opts_box),
	                   GTK_WIDGET(njp->indexing_backend_opts_widget),
	                   FALSE, FALSE, 0);
	gtk_widget_show_all(njp->indexing_backend_opts_widget);
}


static GtkWidget *make_job_opts(struct crystfelproject *proj,
                                struct new_job_params *njp)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	GtkWidget *scroll;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Job name:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	njp->job_title_entry = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->job_title_entry),
	                   TRUE, TRUE, 2.0);
	if ( proj->indexing_new_job_title != NULL ) {
		gtk_entry_set_text(GTK_ENTRY(njp->job_title_entry),
		                   proj->indexing_new_job_title);
	}

	label = gtk_label_new("This name will be used for a working subfolder");
	gtk_label_set_markup(GTK_LABEL(label),
	                     "<i>This name will be used for a working subfolder</i>");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	gtk_entry_set_placeholder_text(GTK_ENTRY(njp->job_title_entry),
	                               "indexing-trial-1");

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   TRUE, TRUE, 0);
	label = gtk_label_new("Notes:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);
	njp->job_notes_text = gtk_text_view_new();
	scroll = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_add(GTK_CONTAINER(scroll), njp->job_notes_text);
	gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroll),
	                                    GTK_SHADOW_ETCHED_IN);
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(scroll),
	                   TRUE, TRUE, 2.0);

	label = gtk_label_new("The notes above will be placed in the job's folder as 'notes.txt'");
	gtk_label_set_markup(GTK_LABEL(label),
	                     "<i>The notes above will be placed in the job's folder as 'notes.txt'</i>");
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	return box;
}


static GtkWidget *make_backend_opts(struct new_job_params *njp)
{
	GtkWidget *box;
	GtkWidget *hbox;
	GtkWidget *label;
	int i;

	box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 8);
	gtk_container_set_border_width(GTK_CONTAINER(box), 8);

	hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 8);
	gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(hbox),
	                   FALSE, FALSE, 0);
	label = gtk_label_new("Batch system:");
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(label),
	                   FALSE, FALSE, 0);

	njp->indexing_backend_combo = gtk_combo_box_text_new();
	gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(njp->indexing_backend_combo),
	                   FALSE, FALSE, 0);

	for ( i=0; i<njp->proj->n_backends; i++ ) {
		gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(njp->indexing_backend_combo),
		                          njp->proj->backends[i].name,
		                          njp->proj->backends[i].friendly_name);
	}

	njp->indexing_backend_opts_box = gtk_box_new(GTK_ORIENTATION_VERTICAL,
	                                              0);
	gtk_box_pack_start(GTK_BOX(box),
	                   GTK_WIDGET(njp->indexing_backend_opts_box),
	                   FALSE, FALSE, 0);
	njp->indexing_backend_opts_widget = NULL;

	/* njp->indexing_backend_opts{_box} must exist before the following */
	g_signal_connect(G_OBJECT(njp->indexing_backend_combo), "changed",
	                 G_CALLBACK(indexing_backend_changed_sig), njp);
	gtk_combo_box_set_active(GTK_COMBO_BOX(njp->indexing_backend_combo),
	                         njp->proj->indexing_backend_selected);

	return box;
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


static void free_new_job_params(gpointer njp, GClosure *closure)
{
	free(njp);
}


gint index_all_sig(GtkWidget *widget, struct crystfelproject *proj)
{
	GtkWidget *dialog;
	GtkWidget *content_area;
	GtkWidget *vbox;
	GtkWidget *backend_page;
	GtkWidget *job_page;
	struct new_job_params *njp;

	if ( proj->indexing_opts != NULL ) return FALSE;

	njp = malloc(sizeof(struct new_job_params));
	if ( njp == NULL ) return FALSE;

	njp->proj = proj;

	dialog = gtk_dialog_new_with_buttons("Index all frames",
	                                     GTK_WINDOW(proj->window),
	                                     GTK_DIALOG_DESTROY_WITH_PARENT,
	                                     "Cancel", GTK_RESPONSE_CANCEL,
	                                     "Run", GTK_RESPONSE_OK,
	                                     NULL);

	g_signal_connect_data(G_OBJECT(dialog), "response",
	                      G_CALLBACK(index_all_response_sig),
	                      njp, free_new_job_params, 0);

	vbox = gtk_vbox_new(FALSE, 0.0);
	content_area = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
	gtk_container_add(GTK_CONTAINER(content_area), vbox);
	gtk_container_set_border_width(GTK_CONTAINER(content_area), 8);

	proj->indexing_opts = crystfel_indexing_opts_new();
	gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(proj->indexing_opts),
	                   FALSE, FALSE, 8.0);
	set_indexing_opts(proj,
	                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));

	job_page = make_job_opts(proj, njp);
	gtk_notebook_prepend_page(GTK_NOTEBOOK(proj->indexing_opts),
	                          job_page,
	                          gtk_label_new("Job name"));

	backend_page = make_backend_opts(njp);
	gtk_notebook_append_page(GTK_NOTEBOOK(proj->indexing_opts),
	                         backend_page,
	                         gtk_label_new("Cluster/batch system"));

	gtk_dialog_set_default_response(GTK_DIALOG(dialog),
	                                GTK_RESPONSE_OK);
	gtk_widget_show_all(dialog);

	return FALSE;
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
	char *methods;
	int i;
	int err;
	TakeTwoOptions *taketwoopts;
	XGandalfOptions *xgandalf_opts;
	PinkIndexerOptions *pinkIndexer_opts;
	FelixOptions *felix_opts;

	if ( proj->indexing_params.cell_file != NULL ) {
		cell = load_cell_from_file(proj->indexing_params.cell_file);
	} else {
		cell = NULL;
	}

	update_peaks(proj);

	if ( proj->indexing_params.indexing_methods == NULL ) {
		methods = detect_indexing_methods(cell);
		STATUS("Auto-detected indexng methods: %s\n",
		       methods);
	} else {
		methods = strdup(proj->indexing_params.indexing_methods);
	}

	/* Get default options for the indexing methods.
	 * The GUI current does not allow them to be changed */
	default_method_options(&taketwoopts,
	                       &xgandalf_opts,
	                       &pinkIndexer_opts,
	                       &felix_opts);

	ipriv = setup_indexing(methods, cell, proj->dtempl,
	                       proj->indexing_params.tols,
	                       indexing_flags(&proj->indexing_params),
	                       taketwoopts, xgandalf_opts,
	                       pinkIndexer_opts, felix_opts);
	free(methods);

	index_pattern(proj->cur_image, ipriv);

	for ( i=0; i<proj->cur_image->n_crystals; i++ ) {
		crystal_set_profile_radius(proj->cur_image->crystals[i], 0.02e9);
		crystal_set_mosaicity(proj->cur_image->crystals[i], 0.0);
		if ( refine_radius(proj->cur_image->crystals[i],
		                   proj->cur_image) )
		{
			ERROR("WARNING: Radius determination failed\n");
		}
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

	STATUS("Number of crystals: %i\n",
	       proj->cur_image->n_crystals);
	for ( i=0; i<proj->cur_image->n_crystals; i++ ) {
		double a, b, c, al, be, ga;
		cell_get_parameters(crystal_get_cell(proj->cur_image->crystals[i]),
		                    &a, &b, &c, &al, &be, &ga);
		STATUS(" %2i: %.2f %.2f %.2f A,  %.2f %.2f %.2f deg\n",
		       i, a*1e10, b*1e10, c*1e10,
		       rad2deg(al), rad2deg(be), rad2deg(ga));
	}

}


static void index_one_response_sig(GtkWidget *dialog, gint resp,
                                   struct crystfelproject *proj)
{
	GtkWidget *w;

	if ( resp == GTK_RESPONSE_OK ) {
		gtk_combo_box_set_active(GTK_COMBO_BOX(proj->results_combo),
		                         0);
		get_indexing_opts(proj,
		                  CRYSTFEL_INDEXING_OPTS(proj->indexing_opts));
		run_indexing_once(proj);
	}

	gtk_widget_destroy(dialog);
	proj->indexing_opts = NULL;
	w =  gtk_ui_manager_get_widget(proj->ui, "/ui/mainwindow/view/refls");
	gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w), 1);
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


static void add_arg(char **args, int pos, const char *label)
{
	args[pos] = strdup(label);
}


static void add_arg_float(char **args, int pos, const char *label,
                          float val)
{
	char *str = malloc(64);
	if ( str == NULL ) return;
	snprintf(str, 63, "--%s=%f", label, val);
	args[pos] = str;
}


static void add_arg_int(char **args, int pos, const char *label,
                        int val)
{
	char *str = malloc(64);
	if ( str == NULL ) return;
	snprintf(str, 63, "--%s=%i", label, val);
	args[pos] = str;
}


GFile *get_crystfel_path_gfile()
{
	GFile *self;
	GFileInfo *self_info;
	const char *self_target;
	GFile *tar;
	GFile *parent_dir;
	GError *error = NULL;

	self = g_file_new_for_path("/proc/self/exe");
	self_info = g_file_query_info(self, "standard",
	                              G_FILE_QUERY_INFO_NONE,
	                              NULL, &error);
	if ( self_info == NULL ) return NULL;

	self_target = g_file_info_get_symlink_target(self_info);
	if ( self_target == NULL ) return NULL;

	tar = g_file_new_for_path(self_target);
	if ( tar == NULL ) return NULL;

	parent_dir = g_file_get_parent(tar);
	if ( parent_dir == NULL ) return NULL;

	g_object_unref(self);
	g_object_unref(self_info);
	g_object_unref(tar);

	return parent_dir;
}


char *get_crystfel_path_str()
{
	char *path;
	GFile *crystfel_path = get_crystfel_path_gfile();
	if ( crystfel_path == NULL ) return NULL;
	path = g_file_get_path(crystfel_path);
	g_object_unref(crystfel_path);
	return path;
}


static char *get_indexamajig_exe()
{
	GFile *crystfel_path;
	char *indexamajig_path;
	GFile *indexamajig;

	crystfel_path = get_crystfel_path_gfile();
	if ( crystfel_path == NULL ) return NULL;

	indexamajig = g_file_get_child(crystfel_path, "indexamajig");
	if ( indexamajig == NULL ) return NULL;

	indexamajig_path = g_file_get_path(indexamajig);
	g_object_unref(indexamajig);
	g_object_unref(crystfel_path);

	return indexamajig_path;
}


char **indexamajig_command_line(const char *geom_filename,
                                const char *n_thread_str,
                                const char *files_list,
                                const char *stream_filename,
                                struct peak_params *peak_search_params,
                                struct index_params *indexing_params)
{
	char **args;
	char tols[2048];
	char *indexamajig_path;
	int n_args = 0;

	args = malloc(64*sizeof(char *));
	if ( args == NULL ) return NULL;

	indexamajig_path = get_indexamajig_exe();
	if ( indexamajig_path == NULL ) {
		ERROR("Couldn't determine indexamajig path. "
		      "This is OK provided the executable path is set "
		      "correctly.\n");
		indexamajig_path = strdup("indexamajig");
	}

	/* The basics */
	add_arg(args, n_args++, indexamajig_path);
	add_arg(args, n_args++, "-i");
	add_arg(args, n_args++, files_list);
	if ( geom_filename != NULL ) {
		add_arg(args, n_args++, "-g");
		add_arg(args, n_args++, geom_filename);
	}
	add_arg(args, n_args++, "-o");
	add_arg(args, n_args++, stream_filename);
	add_arg(args, n_args++, "-j");
	add_arg(args, n_args++, n_thread_str);

	/* Peak search */
	add_arg(args, n_args++, "--peaks");
	add_arg(args, n_args++, str_peaksearch(peak_search_params->method));

	if ( peak_search_params->method == PEAK_ZAEF ) {
		add_arg_float(args, n_args++, "threshold",
		              peak_search_params->threshold);
		add_arg_float(args, n_args++, "min-squared-gradient",
		              peak_search_params->min_sq_gradient);
		add_arg_float(args, n_args++, "min-snr",
		              peak_search_params->min_snr);

	} else if ( peak_search_params->method == PEAK_PEAKFINDER8 ) {
		add_arg_float(args, n_args++, "threshold",
		              peak_search_params->threshold);
		add_arg_float(args, n_args++, "min-snr",
		              peak_search_params->min_snr);
		add_arg_int(args, n_args++, "min-pix-count",
		            peak_search_params->min_pix_count);
		add_arg_int(args, n_args++, "max-pix-count",
		            peak_search_params->max_pix_count);
		add_arg_int(args, n_args++, "local-bg-radius",
		            peak_search_params->local_bg_radius);
		add_arg_int(args, n_args++, "min-res",
		            peak_search_params->min_res);
		add_arg_int(args, n_args++, "max-res",
		            peak_search_params->max_res);
	}

	if ( indexing_params->min_peaks > 0 ) {
		add_arg_int(args, n_args++, "min-peaks",
		            indexing_params->min_peaks);
	}

	/* Indexing */
	if ( indexing_params->indexing_methods != NULL ) {
		add_arg(args, n_args++, "--indexing");
		add_arg(args, n_args++, indexing_params->indexing_methods);
	}
	if ( indexing_params->cell_file != NULL ) {
		add_arg(args, n_args++, "-p");
		add_arg(args, n_args++, indexing_params->cell_file);
	}
	snprintf(tols, 2048, "--tolerance=%f,%f,%f,%f,%f,%f",
	         rad2deg(indexing_params->tols[0]),
	         rad2deg(indexing_params->tols[1]),
	         rad2deg(indexing_params->tols[2]),
	         rad2deg(indexing_params->tols[3]),
	         rad2deg(indexing_params->tols[4]),
	         rad2deg(indexing_params->tols[5]));
	add_arg(args, n_args++, tols);
	if ( indexing_params->multi ) add_arg(args, n_args++, "--multi");
	if ( indexing_params->no_refine ) add_arg(args, n_args++, "--no-refine");
	if ( indexing_params->no_retry ) add_arg(args, n_args++, "--no-retry");
	if ( indexing_params->no_peak_check ) add_arg(args, n_args++, "--no-peak-check");
	if ( indexing_params->no_cell_check ) add_arg(args, n_args++, "--no-cell-check");

	/* Integration */
	add_arg(args, n_args++, "--integration");
	add_arg(args, n_args++, indexing_params->integration_method);
	if ( indexing_params->overpredict ) args[n_args++] = "--overpredict";
	if ( !isinf(indexing_params->push_res) ) {
		add_arg_float(args, n_args++, "push-res",
		              indexing_params->push_res);
	}

	args[n_args] = NULL;
	return args;
}
